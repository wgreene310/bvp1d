#include <stdio.h>

#include <kinsol/kinsol.h>
#include <kinsol/kinsol_dense.h>
#if USE_KLU
#include <kinsol/kinsol_klu.h>
#endif

#include "BVPSolverStats.h"


BVPSolverStats::BVPSolverStats(bool recordStats) :
recordStats(recordStats)
{
  numFuncEvals = numSolveIterations = numJacEval = 0;
  numFuncEvalJac = 0;
  numPtsCurrentMesh = 0;
  numMeshUpdates = 0;
}

void BVPSolverStats::update(void *kinsolMemory, int numMeshPts) {
  if (!recordStats) return;
  numMeshUpdates++;
  numPtsCurrentMesh = numMeshPts;
  long nfevals;
  KINGetNumFuncEvals(kinsolMemory, &nfevals);
  numFuncEvals += nfevals;
  long nSolvIters;
  KINGetNumNonlinSolvIters(kinsolMemory, &nSolvIters);
  numSolveIterations += nSolvIters;
  long numJac, numFuncJac = 0;
#if USE_KLU
  KINSlsGetNumJacEvals(kinsolMemory, &numJac);
#else
  KINDlsGetNumJacEvals(kinsolMemory, &numJac);
  KINDlsGetNumFuncEvals(kinsolMemory, &numFuncJac);
#endif
  numJacEval += numJac;
  numFuncEvalJac += numFuncJac;
}

void BVPSolverStats::incrementJacobianFuncCalls(int numCalls) {
  numFuncEvalJac += numCalls;
}

void BVPSolverStats::print() const
{
  if (!recordStats) return;
  printf("Number of mesh refinements = %d.\n", numMeshUpdates);
  printf("Number of points in final mesh = %d.\n", numPtsCurrentMesh);
  printf("Number of residual function evaluations = %d\n", numFuncEvals);
  printf("Number of nonlinear solver iterations = %d\n", numSolveIterations);
  printf("%d Jacobian evaluations requiring %d residual function evaluations.\n",
    numJacEval, numFuncEvalJac);
}