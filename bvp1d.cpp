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

#include <boost/algorithm/string.hpp>

#include <mex.h>

#include "BVP1DImpl.h"
#include "MexInterface.h"
#include "BVP1dOptions.h"
#include "BVP1dException.h"

#define FUNC_NAME "bvp1d"

namespace {
  BVP1dOptions getOptions(const mxArray *opts) {
    if (!mxIsStruct(opts))
      mexErrMsgIdAndTxt("pde1d:optins_type",
      "The options argument to " FUNC_NAME " must be a struct.");
    BVP1dOptions bvpOpts;
    int n = mxGetNumberOfFields(opts);
    for (int i = 0; i < n; i++) {
      const char *ni = mxGetFieldNameByNumber(opts, i);
      mxArray *val = mxGetFieldByNumber(opts, 0, i);
      if (mxIsEmpty(val)) continue;
      if (boost::iequals(ni, "abstol"))
        bvpOpts.setAbsTol(mxGetScalar(val));
      else if (boost::iequals(ni, "reltol"))
        bvpOpts.setRelTol(mxGetScalar(val));
      else if (boost::iequals(ni, "nmax"))
        bvpOpts.setNMax((int) mxGetScalar(val));
      else if (boost::iequals(ni, "stats")) {
        const int buflen = 1024;
        char buf[buflen];
        mxGetString(val, buf, buflen);
        bool doStats;
        if (boost::iequals(buf, "on"))
          doStats = true;
        else if (boost::iequals(buf, "off"))
          doStats = false;
        else
          mexErrMsgIdAndTxt("bvp1d:invalidStats",
          "The value of the \"Stats\" option must be either \"On\" or \"Off\".");
        bvpOpts.setPrintStats(doStats);
      }
      else {
        char msg[1024];
        sprintf(msg, "The options argument contains the field \"%s\".\n"
          "This is not a currently-supported option and will be ignored.",
          ni);
        mexWarnMsgIdAndTxt("bvp1d:unknown_option", msg);
      }
    }
    return bvpOpts;
  }

  void complexArgError(const char *fldName) {
    char msg[256];
    sprintf(msg, "The \"%s\" structure field is complex. " FUNC_NAME
      " does not currently support complex equations.", fldName);
    throw BVP1dException("bvp1d:complex", msg);
  }
}

class MexBVP : public BVP1DImpl::BVPDefn {
public:
  MexBVP(RealMatrix &yInit, RealVector &param, MexInterface &mexInt, 
    const mxArray *odeFuncPtr, const mxArray *bcFuncPtr) :
    yInit(yInit), parameters(param), mexInt(mexInt), 
    odeFuncPtr(odeFuncPtr), bcFuncPtr(bcFuncPtr) {
    numDepVars = yInit.rows();
    mxX = mxCreateDoubleScalar(0);
    mxY = mxCreateDoubleMatrix(numDepVars, 1, mxREAL);
    mxY0 = mxCreateDoubleMatrix(numDepVars, 1, mxREAL);
    mxYn = mxCreateDoubleMatrix(numDepVars, 1, mxREAL);
    numParameters = parameters.size();
    if (numParameters)
      mxP = mxCreateDoubleMatrix(numParameters, 1, mxREAL);
    else
      mxP = 0;
  }
  ~MexBVP() {
    mxDestroyArray(mxX);
    mxDestroyArray(mxY);
    mxDestroyArray(mxY0);
    mxDestroyArray(mxYn);
    if (mxP)
      mxDestroyArray(mxP);
  }
  virtual void odeFunc(double x, const RealVector &y,
    const RealVector &p, RealVector &fVec) {
    *mxGetPr(mxX) = x;
    std::copy_n(y.data(), numDepVars, mxGetPr(mxY));
    int numInArgs = 3;
    if (numParameters) {
      std::copy_n(p.data(), numParameters, mxGetPr(mxP));
      numInArgs++;
    }
    const mxArray *inArgs[] = { odeFuncPtr, mxX, mxY, mxP };
    RealVector *outArgs[] = { &fVec };
    mexInt.callMatlab(inArgs, numInArgs, outArgs, 1);
  }
  virtual void bcFunc(const RealVector &y0, const RealVector &yn,
    const RealVector &p, RealVector &g) {
    std::copy_n(y0.data(), numDepVars, mxGetPr(mxY0));
    std::copy_n(yn.data(), numDepVars, mxGetPr(mxYn));
    int numInArgs = 3;
    if (numParameters) {
      std::copy_n(p.data(), numParameters, mxGetPr(mxP));
      numInArgs++;
    }
    const mxArray *inArgs[] = { bcFuncPtr, mxY0, mxYn, mxP };
    RealVector *outArgs[] = { &g };
    mexInt.callMatlab(inArgs, numInArgs, outArgs, 1);
  }
private:
  RealMatrix &yInit;
  RealVector &parameters;
  MexInterface &mexInt;
  int numDepVars, numParameters;
  mxArray *mxX, *mxY, *mxY0, *mxYn, *mxP;
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
  const char *fieldNames[] = { "x", "y", "solver", "yp", "parameters" };
  const int numReqFields = 2;
    mxArray *mxFlds[numReqFields];
  for (int i = 0; i < numReqFields; i++) {
    const char *fn = fieldNames[i];
    mxFlds[i] = mxGetField(solinit, 0, fn);
    if (!mxFlds[i]) {
      char msg[80];
      sprintf(msg, "solinit is missing the required %s array.", fn);
      mexErrMsgIdAndTxt("bvp1d:solinit_field_missing", msg);
    }
    if (mxIsComplex(mxFlds[i]))
      complexArgError(fn);
  }

  int numFields = numReqFields + 2;
  RealVector parameters;
  mxArray *mxParams = mxGetField(solinit, 0, "parameters");
  MexInterface mexInt;
  if (mxParams) {
    parameters = mexInt.fromMxArrayVec(mxParams);
    if (mxIsComplex(mxParams))
      complexArgError("parameters");
    numFields++;
  }

  mxArray *sol = mxCreateStructMatrix(1, 1, numFields, &fieldNames[0]);

  try {
    RealVector mesh = mexInt.fromMxArrayVec(mxFlds[0]);
    //cout << mesh << endl;
    RealMatrix yInit = mexInt.fromMxArray(mxFlds[1]);
    if (yInit.cols() != mesh.size())
      mexErrMsgIdAndTxt("bvp1d:solinit_x_y_inconsistent", 
      "The number of columns in y must equal the length of the x array.");
    MexBVP bvpDef(yInit, parameters, mexInt, prhs[0], prhs[1]);
    BVP1dOptions opts;
    if (nrhs == 4)
      opts = getOptions(prhs[3]);
    BVP1DImpl bvp(bvpDef, mesh, yInit, parameters, opts);
    RealMatrix y, yPrime;
    RealVector p;
    int err = bvp.solve(y, yPrime, p);
    if (err)
      mexErrMsgIdAndTxt("bvp1d:solve_failure",
      "Unable to solve BVP.");
    mxSetField(sol, 0, "solver", mxCreateString("bvp4c"));
    mxSetField(sol, 0, "x", mexInt.toMxArray(bvp.getMesh().transpose()));
    mxSetField(sol, 0, "y", mexInt.toMxArray(y));
    mxSetField(sol, 0, "yp", mexInt.toMxArray(yPrime));
    mxSetField(sol, 0, "error", mexInt.toMxArray(bvp.getError()));
    if (mxParams)
      mxSetField(sol, 0, "parameters", mexInt.toMxArray(p));
  }
  catch (const BVP1dException &ex) {
    mexErrMsgIdAndTxt(ex.getId(), ex.what());
  }
  catch (const std::exception &ex) {
    mexErrMsgIdAndTxt("bvp1d:exception", ex.what());
  }
  catch (...) {
    mexErrMsgIdAndTxt("bvp1d:internal_err", "Internal error in bvp1d.\n");
  }
  plhs[0] = sol;
}