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

#pragma once

#include <nvector/nvector_serial.h>
#include <kinsol/kinsol.h>
//#include <kinsol/kinsol_klu.h>

#include <Eigen/SparseCore>

#if 0
struct SunSparseMat : public _SlsMat {
  SunSparseMat(Eigen::SparseMatrix<double> &eMat) {
    M = static_cast<int>(eMat.rows());
    N = static_cast<int>(eMat.cols());
    NNZ = static_cast<int>(eMat.nonZeros());
    data = eMat.valuePtr();
    rowvals = eMat.innerIndexPtr();
    colptrs = eMat.outerIndexPtr();
  }
};
#endif

class FiniteDiffJacobian {
public:
  typedef Eigen::SparseMatrix<double> SparseMat;
  typedef Eigen::Map<SparseMat> SparseMap;
  FiniteDiffJacobian(SparseMat &jacPattern);
  ~FiniteDiffJacobian();
  void calcJacobian(N_Vector uu, N_Vector r,
    KINSysFn rf, void *userData, SparseMap Jac, int *numFuncEvals=0);
private:
  void copyIndices(int *outerInd, int *innerInd);
  int neq, nnz;
  Eigen::VectorXi indrow, jpntr, ngrp;
  int maxgrp, mingrp;
};


