#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include <mex.h>

class MexInterface {
public:
  typedef Eigen::VectorXd RealVector;
  MexInterface(int maxOutArgs=1);
  ~MexInterface();
  void print(const mxArray *a, const char *name);
  std::string getFuncNameFromHandle(const mxArray *fh);
  void callMatlab(const mxArray *inArgs[], int nargin,
    RealVector *outArgs[], int nargout);
  template<class T>
  mxArray *toMxArray(const T &a)
  {
    mxArray *ma = mxCreateDoubleMatrix(a.rows(), a.cols(), mxREAL);
    double *dest = mxGetPr(ma);
    const double *src = a.data();
    std::copy_n(src, a.cols()*a.rows(), dest);
    return ma;
  }
  Eigen::MatrixXd fromMxArray(const mxArray *a);
  Eigen::VectorXd fromMxArrayVec(const mxArray *a);
private:
  const int maxOutArgs;
  std::vector<mxArray*> matOutArgs;
};

