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

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_klu.h>

#include <FDJacobian.h>

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

template<class T>
void BVP1DImpl::calcPhi(const T &y, T &phi)
{
  Eigen::Map<const Eigen::MatrixXd> yMat(y.data(), numDepVars, numNodes);
  MapMat phiMat(phi.data(), numDepVars, numNodes);
  bvp.bcFunc(yMat.col(0), yMat.col(numNodes - 1), parameters, gVec);
  phiMat.col(0) = gVec;
  bvp.odeFunc(mesh[0], yMat.col(0), parameters, fim1);
  for (int i = 1; i < numNodes; i++) {
    double hi = mesh[i] - mesh[i - 1];
    auto &yi = yMat.col(i);
    bvp.odeFunc(mesh[i], yMat.col(i), parameters, fi);
    auto &yim1 = yMat.col(i - 1);
    double xm2 = mesh[i - 1] + hi / 2.;
    yim2 = (yim1 + yi) / 2. - hi / 8.*(fi - fim1);
    bvp.odeFunc(mesh[i], yim2, parameters, fim2);
    phiMat.col(i) = yi - yim1 - hi / 6.*(fim1 + 4 * fim2 + fi);
    fim1 = fi;
  }
#if 0
  cout << "phi\n" << phiMat << endl;
#endif
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
    printf("Vector: %s, loc(v)=%d, loc(v.data)=%d\n", name, v, NV_DATA_S(v));
  }

  int funcBathe(N_Vector u, N_Vector f, void *user_data) {

    BVP1DImpl *bvp = (BVP1DImpl*) user_data;

    int neq = NV_LENGTH_S(u);
    MapVec uVec(NV_DATA_S(u), neq);
    
    prtLocVec(u, "u_res");

    MapVec fVec(NV_DATA_S(f), neq);

    //cout << bvp->J << endl;
    fVec = bvp->J*uVec - bvp->rhs;
    //cout << "J*u=" << (bvp->J*uVec).transpose() << endl;
    cout << "u=" << uVec.transpose() << endl;
    cout << "f=" << fVec.transpose() << endl;

    return 0;
  }

  int jacBathe(N_Vector u, N_Vector f, SlsMat Jac, void *user_data,
    N_Vector tmp1, N_Vector tmp2) {
    prtLocVec(u, "u_jac");
    BVP1DImpl *bvp = (BVP1DImpl*) user_data;
    SparseMat &jacEig = bvp->J;
    int neq = NV_LENGTH_S(u);
    int nnz = jacEig.nonZeros();
    std::copy_n(jacEig.outerIndexPtr(), neq + 1, Jac->colptrs);
    std::copy_n(jacEig.innerIndexPtr(), nnz, Jac->rowvals);
    std::copy_n(jacEig.valuePtr(), nnz, Jac->data);
    Jac->NNZ = nnz;
    //PrintSparseMat(Jac);

    return 0;
  }

  static int check_flag(void *flagvalue, char *funcname, int opt)
  {
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
      fprintf(stderr,
        "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      return(1);
    }

    /* Check if flag < 0 */
    else if (opt == 1) {
      errflag = (int *)flagvalue;
      if (*errflag < 0) {
        fprintf(stderr,
          "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
        return(1);
      }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
      fprintf(stderr,
        "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      return(1);
    }

    return(0);
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


BVP1DImpl::BVP1DImpl(BVPDefn &bvp, RealVector &mesh, RealMatrix &yInit,
  RealVector &parameters) :
  bvp(bvp), mesh(mesh), yInit(yInit), parameters(parameters)
{
  if (yInit.cols() != mesh.size())
    throw BVP1dException("bvp1d:solinit_x_y_inconsistent",
    "The number of columns in y must equal the length of the x array.");
  numDepVars = yInit.rows();
  numNodes = mesh.rows();
  numParams = parameters.size();
  gVec.resize(numDepVars + numParams);
  fi.resize(numDepVars);
  fim1.resize(numDepVars);
  fim2.resize(numDepVars);
  yim2.resize(numDepVars);
}


BVP1DImpl::~BVP1DImpl()
{
}

int BVP1DImpl::batheTest() {
  typedef Eigen::Triplet<double> Triplet;
  Eigen::Triplet<double> triplets[] = { Triplet(0, 0, 2), Triplet(0, 1, -2), Triplet(0, 4, -1),
    Triplet(1, 1, 3), Triplet(1, 2, -2),
    Triplet(2, 2, 5), Triplet(2, 3, -3), Triplet(3, 3, 10),
    Triplet(3, 4, 4), Triplet(4, 4, 10) };
  const int n = 5;
  
  //printf("nnz=%d\n", nnz);
  std::vector<Triplet> eigTriplets;
  for (int i = 0; i<sizeof(triplets) / sizeof(triplets[0]); i++) {
    const Triplet &t = triplets[i];
    int r = t.row();
    int c = t.col();
    eigTriplets.push_back(t);
    if (r != c)
      eigTriplets.push_back(Triplet(c, r, t.value()));
  }

  J.resize(n, n);
  J.setFromTriplets(eigTriplets.begin(), eigTriplets.end());
  const int nnz = J.nonZeros();

  //cout << J << endl;
  //print(J);

  rhs.resize(n);
  rhs << 0, 1, 0, 0, 0;
  Eigen::VectorXd uExact(n);
  uExact << 636., 619., 292., 74., 34.;

  //cout << J*uExact << endl;

  N_Vector u = N_VNew_Serial(n);
  prtLocVec(u, "u_main");
  MapVec uVec(NV_DATA_S(u), n);
  uVec.setConstant(1);
  
  void *kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);
  int flag = KINInit(kmem, funcBathe, u);
  if (check_flag(&flag, "KINInit", 1)) return(1);
  int ier = KINSetUserData(kmem, this);
  if (check_flag(&ier, "KINSetUserData", 1)) return(1);
  /* Specify stopping tolerance based on residual */
  double fnormtol = 1e-8;
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
  flag = KINKLU(kmem, n, nnz);
  if (check_flag(&flag, "KINKLU", 1)) return(1);
  flag = KINSlsSetSparseJacFn(kmem, jacBathe);
  if (check_flag(&flag, "KINSlsSetSparseJacFn", 1)) return(1);

  N_Vector scale = N_VNew_Serial(n);
  N_VConst_Serial(1, scale);

  /* Call main solver */
  int strat = KIN_PICARD;
  strat = KIN_NONE;
  strat = KIN_LINESEARCH;
  flag = KINSol(kmem,           /* KINSol memory block */
    u,         /* initial guess on input; solution vector */
    strat,     /* global strategy choice */
    scale,          /* scaling vector, for the variable cc */
    scale);         /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1)) return(1);

  prtLocVec(u, "u_main");

  cout << uVec << endl;

  return 0;
}

int BVP1DImpl::solve(Eigen::MatrixXd &solMat)
{
  const int neq = numDepVars*numNodes + numParams;
  solMat.resize(numDepVars, numNodes);

  N_Vector u = N_VNew_Serial(neq);
  MapMat uMat(NV_DATA_S(u), numDepVars, numNodes); // FIXME!!
  //uMat.setConstant(1);
  uMat = yInit;
  //cout << "yInit\n" << uMat << endl;

  void *kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);
  int flag = KINInit(kmem, funcKinsol, u);
  if (check_flag(&flag, "KINInit", 1)) return(1);
  int ier = KINSetUserData(kmem, this);
  if (check_flag(&ier, "KINSetUserData", 1)) return(1);
  /* Specify stopping tolerance based on residual */
  double fnormtol = 1e-5, scsteptol=fnormtol;
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
  flag = KINSetScaledStepTol(kmem, scsteptol);
  if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);
  flag = KINDense(kmem, neq);
  if (check_flag(&flag, "KINDense", 1)) return(1);

  N_Vector scale = N_VNew_Serial(neq);
  N_VConst_Serial(1, scale);

  /* Call main solver */
  int strat = KIN_PICARD;
  strat = KIN_NONE;
  strat = KIN_LINESEARCH;
  flag = KINSol(kmem,           /* KINSol memory block */
    u,         /* initial guess on input; solution vector */
    strat,     /* global strategy choice */
    scale,          /* scaling vector, for the variable cc */
    scale);         /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1)) return(1);

  //cout << uMat << endl;

  solMat = uMat;

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
