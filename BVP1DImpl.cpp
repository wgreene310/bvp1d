// Copyright (C) 2016 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <limits>

using std::cout;
using std::endl;

#include "BVP1DImpl.h"
#include "BVP1dException.h"
#include "cubicInterp.h"
#include "BVP1dOptions.h"
#include "GaussLobattoIntRule.h"
#include "SunVector.h"

#define USE_LAPACK 1

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_dense.h>
#if USE_LAPACK
#include <kinsol/kinsol_lapack.h>
#endif

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


template<class T1, class T2>
void BVP1DImpl::calcPhi(const T1 &y, T2 &phi)
{
  const int numNodes = mesh.size();
  Eigen::Map<const Eigen::MatrixXd> yMat(y.data(), numDepVars, numNodes);
  if (numParams) {
    parameters = y.segment(numDepVars*numNodes, numParams);
    //cout << "p=" << parameters.transpose() << endl;
  }
  int phiOff = gVec.rows();
  bvp.bcFunc(yMat.col(0), yMat.col(numNodes - 1), parameters, gVec);
  phi.head(phiOff) = gVec;
  bvp.odeFunc(mesh[0], yMat.col(0), parameters, fim1);
  fRHS.col(0) = fim1;
  for (int i = 1; i < numNodes; i++) {
    double hi = mesh[i] - mesh[i - 1];
    const auto &yi = yMat.col(i);
    bvp.odeFunc(mesh[i], yMat.col(i), parameters, fi);
    const auto &yim1 = yMat.col(i - 1);
    double xm2 = mesh[i - 1] + hi / 2.;
    yim2 = (yim1 + yi) / 2. - hi / 8.*(fi - fim1);
    bvp.odeFunc(xm2, yim2, parameters, fim2);
    phi.segment(phiOff, numDepVars) = yi - yim1 - hi / 6.*(fim1 + 4 * fim2 + fi);
    fim1 = fi;
    fRHS.col(i) = fi;
    phiOff += numDepVars;
  }
#if 0
  cout << "phi\n" << phi << endl;
#endif
}

template<class T>
void BVP1DImpl::calcError(const T &u)
{
  GaussLobattoIntRule intRule(5);
  const int numIntPts = 3;
  GaussLobattoIntRule::IRule intPts[numIntPts];
  for (int i = 0; i < numIntPts; i++)
    intPts[i] = intRule.getPoint(i + 1);

  const int numNodes = mesh.size();
  Eigen::Map<const Eigen::MatrixXd> y(u.data(), numDepVars, numNodes);

  if (numParams) 
    parameters = u.segment(numDepVars*numNodes, numParams);
  int numEl = mesh.size() - 1;
  residualError.resize(numEl);
  residualError.setZero();
  RealVector fim2X(numDepVars);
  double absOvRel = options.getAbsTol() / options.getRelTol();
  double maxErr = 0;
  for (int e = 0; e < numEl; e++) {
    double elemErr = 0;
    double h = mesh(e + 1) - mesh(e);
    double jac = h / 2;
    for (int i = 0; i < numIntPts; i++) {
      double s = intPts[i].xi;
      cubicInterp(y.col(e), fRHS.col(e), y.col(e + 1), fRHS.col(e + 1),
        h, s, yim2, fim2);
#if 0
      cout << "e=" << e << "i=" << i <<  "yim2=" << yim2.transpose() <<
        " fim2=" << fim2.transpose() << endl;
#endif
      double xm2 = mesh[e] + h / 2.*(s+1);
      bvp.odeFunc(xm2, yim2, parameters, fim2X);
#if 0
      cout << "e=" << e << "i=" << i <<  "yim2=" << yim2.transpose() <<
        " fim2X=" << fim2X.transpose() << endl;
#endif
#if 0
      cout << "xm2=" << xm2 << " fim2=" << fim2.transpose() << 
        " fim2X=" << fim2X.transpose() << endl;
      double err = (fim2 - fim2X).cwiseAbs().maxCoeff();
      printf("Element: %2d, pt=%d, err=%12.3e\n", e, i, err);
#endif
      auto denom = fim2X.array().abs().max(absOvRel);
      //cout << "denom=" << denom.transpose() << endl;
      auto errI = ((fim2 - fim2X).array()/denom).matrix();
      double errNorm = errI.dot(errI);
      //elemErr += (fim2 - fim2X).array().square().matrix()*intPts[i].wt*jac;
      elemErr += errNorm*intPts[i].wt*jac;
    }
    residualError(e) = sqrt(elemErr);
  }
}

namespace {

  void print(SparseMat &a) {
    auto ad = a.toDense();
    for (int i = 0; i < a.rows(); i++) {
      for (int j = 0; j < a.cols(); j++)
        printf("%12.4e,", ad(i, j));
      printf("\n");
    }
  }

  void prtLocVec(N_Vector v, const char *name) {
    printf("Vector: %s, loc(v)=%p, loc(v.data)=%p\n", name, v, NV_DATA_S(v));
  }

  void check_flag(void *flagvalue, const char *funcname, int opt)
  {
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
      char msg[256];
      sprintf(msg,
        "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      throw BVP1dException("bvp1d:sundials_mem_alloc", msg);
    }

    /* Check if flag < 0 */
    else if (opt == 1) {
      errflag = (int *)flagvalue;
      if (*errflag < 0) {
        char msg[256];
        sprintf(msg,
          "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
        throw BVP1dException("bvp1d:sundials_error", msg);
      }
    }
  }

  int funcKinsol(N_Vector y, N_Vector phi, void *user_data) {

    BVP1DImpl *bvp = (BVP1DImpl*) user_data;
    int neq = NV_LENGTH_S(y);
    MapVec yVec(NV_DATA_S(y), neq);
    MapVec phiVec(NV_DATA_S(phi), neq);
    bvp->calcPhi(yVec, phiVec);
    return 0;
  }

}


BVP1DImpl::BVP1DImpl(BVPDefn &bvp, RealVector &initMesh, RealMatrix &yInit,
  RealVector &parameters, BVP1dOptions &options) :
  bvp(bvp), initMesh(initMesh), initSolution(yInit), parameters(parameters),
  options(options)
{
  kmem = 0;
  mesh = initMesh;
  currentInitSoln = initSolution;
  if (yInit.cols() != mesh.size())
    throw BVP1dException("bvp1d:solinit_x_y_inconsistent",
    "The number of columns in y must equal the length of the x array.");
  numDepVars = yInit.rows();
  numParams = parameters.size();
  gVec.resize(numDepVars + numParams);
  fi.resize(numDepVars);
  fim1.resize(numDepVars);
  fim2.resize(numDepVars);
  yim2.resize(numDepVars);
#if 0
  cout << "YInit\n" << yInit << endl;
  cout << "pInit=" << parameters << endl;
#endif
}


BVP1DImpl::~BVP1DImpl()
{
  if (kmem) KINFree(&kmem);
}

int BVP1DImpl::solve(Eigen::MatrixXd &solMat, RealMatrix &yPrime, 
  RealVector &paramVec)
{
  int nMax = options.getNMax();
  if (nMax < 0) {
    nMax = (int)std::floor(1000. / numDepVars);
  }
  const double o3 = 1. / 3.;
  double relTol = options.getRelTol();
  double maxErr = 0;
  while (true) {
    const int numNodes = mesh.size();
    int err = solveFixedMesh(solMat, yPrime, paramVec);
    if (err) return err;
    maxErr = residualError.maxCoeff();
    if (maxErr <= relTol)
      return 0;
    // refine mesh
    RealVector newMesh;
    RealMatrix newInitSoln;
    refineMesh(solMat, newMesh, newInitSoln);
    //cout << "newMesh\n" << newMesh << endl;
    //cout << "newInitSoln\n" << newInitSoln.transpose() << endl;
    if (newMesh.size() > nMax)
      break;
    mesh = newMesh;
    currentInitSoln = newInitSoln;
  }

  printf("Maximum residual error of %11.3e exceeds RelTol value of %11.3e.\n",
    maxErr, relTol);

  return 0;
}

void BVP1DImpl::refineMesh(const RealMatrix &sol, RealVector &newMesh,
  RealMatrix &newInitSoln)
{
  const double o3 = 1. / 3.;
  double relTol = options.getRelTol();
  const int numNodes = mesh.size();
  int numEl = numNodes - 1;
  const int maxNewNodes = numNodes + 2 * numEl;
  newMesh.resize(maxNewNodes);
  newInitSoln.resize(numDepVars, maxNewNodes);
  int numNewN = 1;
  newMesh[0] = mesh[0];
  newInitSoln.col(0) = sol.col(0);
  for (int e = 0; e < numEl; e++) {
    double errE = residualError[e];
    double h = mesh[e + 1] - mesh[e];
    if (errE > 100 * relTol) {
      // add two new points in this interval
      double h3 = h / 3;
      newMesh[numNewN] = mesh[e] + h3;
      newMesh[numNewN + 1] = mesh[e] + 2 * h3;
      Eigen::Ref<Eigen::VectorXd> cN = newInitSoln.col(numNewN);
      cubicInterp(sol.col(e), fRHS.col(e),
        sol.col(e + 1), fRHS.col(e + 1),
        h, -o3, cN, fim2);
      Eigen::Ref<Eigen::VectorXd> cN1 = newInitSoln.col(numNewN+1);
      cubicInterp(sol.col(e), fRHS.col(e),
        sol.col(e + 1), fRHS.col(e + 1),
        h, o3, cN1, fim2);
      numNewN += 2;
    }
    else if (errE > relTol) {
      // add one new point in this interval
      double h2 = h / 2;
      newMesh[numNewN] = mesh[e] + h2;
      Eigen::Ref<Eigen::VectorXd> cN = newInitSoln.col(numNewN);
      cubicInterp(sol.col(e), fRHS.col(e),
        sol.col(e + 1), fRHS.col(e + 1), h, 0, cN, fim2);
      numNewN++;
    }
    newMesh[numNewN] = mesh[e + 1];
    newInitSoln.col(numNewN) = sol.col(e + 1);
    numNewN++;
  }
  newMesh.conservativeResize(numNewN);
  newInitSoln.conservativeResize(numDepVars, numNewN);
  //printf("Num nodes: old mesh=%d, new mesh=%d\n", numNodes, numNewN);
}

int BVP1DImpl::solveFixedMesh(Eigen::MatrixXd &solMat, RealMatrix &yPrime, 
  RealVector &paramVec)
{
  const int numNodes = mesh.size();
  const int numYEqns = numDepVars*numNodes;
  const int neq = numYEqns + numParams;
  solMat.resize(numDepVars, numNodes);
  fRHS.resize(numDepVars, numNodes);

  SunVector u(neq);
  double *uData = &u[0];
  MapMat uMat(uData, numDepVars, numNodes);
  uMat = currentInitSoln;
  //cout << "yInit\n" << uMat << endl;
  if (numParams) {
    std::copy_n(parameters.data(), numParams, &uData[numYEqns]);
    paramVec.resize(numParams);
  }

  kmem = KINCreate();
  check_flag((void *) kmem, "KINCreate", 0);
  int flag = KINInit(kmem, funcKinsol, u());
  check_flag(&flag, "KINInit", 1);
  int ier = KINSetUserData(kmem, this);
  check_flag(&ier, "KINSetUserData", 1);
  /* Specify stopping tolerance based on residual */
  double fnormtol = options.getAbsTol();
  double scsteptol = fnormtol;
  flag = KINSetFuncNormTol(kmem, fnormtol);
  check_flag(&flag, "KINSetFuncNormTol", 1);
  flag = KINSetScaledStepTol(kmem, scsteptol);
  check_flag(&flag, "KINSetScaledStepTol", 1);
#if USE_LAPACK
  flag = KINLapackDense(kmem, neq);
  check_flag(&flag, "KINLapackDense", 1);
#else
  flag = KINDense(kmem, neq);
  check_flag(&flag, "KINDense", 1);
#endif

  SunVector scale(neq);
  scale.setConstant(1);

  const double mxnewtstep = 1e5;
  flag = KINSetMaxNewtonStep(kmem, mxnewtstep);
  check_flag(&flag, "KINSetMaxNewtonStep", 1);

  /* Call main solver */
  int strat = KIN_LINESEARCH;
  flag = KINSol(kmem,           /* KINSol memory block */
    u(),         /* initial guess on input; solution vector */
    strat,     /* global strategy choice */
    scale(),          /* scaling vector, for the variable cc */
    scale());         /* scaling vector for function values fval */
  check_flag(&flag, "KINSol", 1);

  //cout << uMat << endl;
  MapVec uVec(uData, neq);
  calcError(uVec);

  solMat = uMat;
  yPrime = fRHS;
  if (numParams)
    std::copy_n(&uData[numYEqns], numParams, paramVec.data());

  KINFree(&kmem);
  kmem = 0;

  return 0;
}

RealVector BVP1DImpl::linspace(double start, double end, int n)
{
  RealVector v(n);
  int nm1 = n - 1;
  double dx = (end - start) / nm1;
  double x = start;
  for (int i = 0; i < n; i++) {
    v[i] = x;
    x += dx;
  }
  return v;
}
