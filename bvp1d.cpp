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

#include <stdio.h>
#include <stdexcept>
#include <iostream>
#include <limits>

using std::cout;
using std::endl;

#include <mex.h>

#include "BVP1DImpl.h"
#include "MexInterface.h"

#define FUNC_NAME "bvp1d"

namespace {

}

class MexBVP : public BVP1DImpl::BVPDefn {
public:
  MexBVP(RealMatrix &yInit, MexInterface &mexInt, 
    const mxArray *odeFuncPtr, const mxArray *bcFuncPtr) :
    yInit(yInit), mexInt(mexInt), odeFuncPtr(odeFuncPtr), bcFuncPtr(bcFuncPtr) {
    numDepVars = yInit.rows();
    mxX = mxCreateDoubleScalar(0);
    mxY = mxCreateDoubleMatrix(numDepVars, 1, mxREAL);
    mxY0 = mxCreateDoubleMatrix(numDepVars, 1, mxREAL);
    mxYn = mxCreateDoubleMatrix(numDepVars, 1, mxREAL);
  }
  ~MexBVP() {
    mxDestroyArray(mxX);
    mxDestroyArray(mxY);
    mxDestroyArray(mxY0);
    mxDestroyArray(mxYn);
  }
  virtual void odeFunc(double x, const RealVector &y,
    const RealVector &p, RealVector &fVec) {
#if 0
    const double c = .2, f = .3;
    fVec[0] = y[1] / c;
    fVec[1] = -f;
#else
    *mxGetPr(mxX) = x;
    std::copy_n(y.data(), numDepVars, mxGetPr(mxY));
    const mxArray *inArgs[] = { odeFuncPtr, mxX, mxY };
    RealVector *outArgs[] = { &fVec };
    mexInt.callMatlab(inArgs, 3, outArgs, 1);
#endif
  }
  virtual void bcFunc(const RealVector &y0, const RealVector &yn,
    const RealVector &p, RealVector &g) {
#if 1
    std::copy_n(y0.data(), numDepVars, mxGetPr(mxY0));
    std::copy_n(yn.data(), numDepVars, mxGetPr(mxYn));
    const mxArray *inArgs[] = { bcFuncPtr, mxY0, mxYn };
    RealVector *outArgs[] = { &g };
    mexInt.callMatlab(inArgs, 3, outArgs, 1);
#else
    g << y0[0], yn[1];
#endif
  }
private:
  RealMatrix &yInit;
  MexInterface &mexInt;
  int numDepVars;
  mxArray *mxX, *mxY, *mxY0, *mxYn;
  const mxArray *odeFuncPtr, *bcFuncPtr;
};

void mexFunction(int nlhs, mxArray*
  plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 3 && nrhs != 4)
    mexErrMsgIdAndTxt("bvp1d:nrhs",
    FUNC_NAME " requires three or four input arguments.");

  for (int i = 0; i < 1; i++) {
    if (!mxIsFunctionHandle(prhs[i])) {
      char msg[80];
      sprintf(msg, "Argument %d is not a function handle.", i + 1);
      mexErrMsgIdAndTxt("bvp1d:arg_not_func", msg);
    }
  }

  const mxArray *solinit = prhs[2];
  if (!mxIsStruct(solinit))
    mexErrMsgIdAndTxt("bvp1d:arg3_not_struct", 
    "Argument three must be a struct.");
  const char *reqFieldNames[] = { "x", "y" };
  const int numFields = sizeof(reqFieldNames) / sizeof(reqFieldNames[0]);
  mxArray *mxFlds[numFields];
  for (int i = 0; i < numFields; i++) {
    const char *fn = reqFieldNames[i];
    mxFlds[i] = mxGetField(solinit, 0, fn);
    if (!mxFlds[i]) {
      char msg[80];
      sprintf(msg, "solinit is missing the required %s array.", fn);
      mexErrMsgIdAndTxt("bvp1d:solinit_field_missing", msg);
    }
  }

  RealVector parameters;
  mxArray *mxParams = mxGetField(solinit, 0, "parameters");

  mxArray *sol = mxCreateStructMatrix(1, 1, numFields, &reqFieldNames[0]);
  MexInterface mexInt;
  if (mxParams)
    parameters = mexInt.fromMxArrayVec(mxParams);
  try {
    RealVector mesh = mexInt.fromMxArrayVec(mxFlds[0]);
    //cout << mesh << endl;
    RealMatrix yInit = mexInt.fromMxArray(mxFlds[1]);
    if (yInit.cols() != mesh.size())
      mexErrMsgIdAndTxt("bvp1d:solinit_x_y_inconsistent", 
      "The number of columns in y must equal the length of the x array.");
    MexBVP bvpDef(yInit, mexInt, prhs[0], prhs[1]);
    BVP1DImpl bvp(bvpDef, mesh, yInit, parameters);
    Eigen::MatrixXd y;
    int err = bvp.solve(y);
    if (err)
      mexErrMsgIdAndTxt("bvp1d:solve_failure",
      "Unable to solve BVP.");
    mxSetField(sol, 0, "x", mexInt.toMxArray(mesh));
    mxSetField(sol, 0, "y", mexInt.toMxArray(y));
  }
  catch (const std::exception &ex) {
    mexErrMsgIdAndTxt("bvp1d:exception", ex.what());
  }
  catch (...) {
    mexErrMsgIdAndTxt("bvp1d:internal_err", "Internal error in bvp1d.\n");
  }
#if 1
  // for now just return the input mesh
  //mxSetField(sol, 0, "x", mxX);
#endif
  plhs[0] = sol;
}