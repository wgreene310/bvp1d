

#include "FiniteDiffJacobian.h"
#include "FDJacobian.h"
#include "SunVector.h"

FiniteDiffJacobian::FiniteDiffJacobian(SparseMat &jacPattern)
{
  neq = jacPattern.rows();
  nnz = jacPattern.nonZeros();
  indrow.resize(nnz);
  jpntr.resize(nnz);
  ngrp.resize(neq);
  Eigen::VectorXi indcol(nnz);
  int ii = 0;
  for (int k = 0; k < jacPattern.outerSize(); ++k)
    for (SparseMat::InnerIterator it(jacPattern, k); it; ++it) {
    indrow[ii] = it.row() + 1;   // row index
    indcol[ii] = it.col() + 1;
    ii++;
    }
  const int liwa = 6 * neq;
  Eigen::VectorXi ipntr(neq + 1), iwa(liwa);
  int info;
  dsm_(&neq, &neq, &nnz, indrow.data(), indcol.data(), ngrp.data(),
    &maxgrp, &mingrp, &info, ipntr.data(), jpntr.data(), iwa.data(), &liwa);
  //printf("info=%d, maxgrp=%d, mingrp=%d\n", info, maxgrp, mingrp);

  sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
}


FiniteDiffJacobian::~FiniteDiffJacobian()
{
}


void FiniteDiffJacobian::calcJacobian(N_Vector uu,  N_Vector r,
  KINSysFn rFunc, void *userData, SlsMat jac, int *numFuncEvals)
{
  typedef Eigen::Map<Eigen::VectorXd> MapVec;
  MapVec u(NV_DATA_S(uu), neq);
  MapVec resvec(NV_DATA_S(r), neq);
  MapVec jacData(jac->data, nnz);

  rFunc(uu, r, userData);
  Eigen::VectorXd u0 = u;
  Eigen::VectorXd r0 = resvec;

  Eigen::VectorXd d(neq), fjacd(neq), fjac(nnz);
  const int col = 1;

  for (int numgrp = 1; numgrp <= maxgrp; numgrp++) {
    for (int j = 0; j < neq; j++) {
      d[j] = 0;
      if (ngrp[j] == numgrp) {
        d[j] = sqrtEps*std::max(u0[j], 1.0);
      }
      u[j] = u0[j] + d[j];
    }
    rFunc( uu, r, userData);
    fjacd = resvec - r0;
    fdjs_(&neq, &neq, &col, indrow.data(), jpntr.data(), ngrp.data(),
      &numgrp, d.data(), fjacd.data(), fjac.data());
  }
  jacData = fjac;
  if (numFuncEvals)
    *numFuncEvals = maxgrp;
 

  // copy the row and column pointers
  indrow.array() -= 1;
  jpntr.array() -= 1;
  std::copy_n(jpntr.data(), neq + 1, jac->colptrs);
  std::copy_n(indrow.data(), nnz, jac->rowvals);
  indrow.array() += 1;
  jpntr.array() += 1;
}
