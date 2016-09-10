#pragma once

#include <nvector/nvector_serial.h>
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_klu.h>

#include <Eigen/SparseCore>

struct SunSparseMat : public _SlsMat {
  SunSparseMat(Eigen::SparseMatrix<double> &eMat) {
    M = eMat.rows();
    N = eMat.cols();
    NNZ = eMat.nonZeros();
    data = eMat.valuePtr();
    rowvals = eMat.innerIndexPtr();
    colptrs = eMat.outerIndexPtr();
  }
};

class FiniteDiffJacobian {
public:
  typedef Eigen::SparseMatrix<double> SparseMat;
  FiniteDiffJacobian(SparseMat &jacPattern);
  ~FiniteDiffJacobian();
  void calcJacobian(N_Vector uu, N_Vector r,
    KINSysFn rf, void *userData, SlsMat Jac, int *numFuncEvals=0);
private:
  int neq, nnz;
  Eigen::VectorXi indrow, jpntr, ngrp;
  int maxgrp, mingrp;
  double sqrtEps;
};


