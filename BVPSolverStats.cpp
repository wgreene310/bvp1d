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
  printf("Number of residual function evaluations = %ld\n", numFuncEvals);
  printf("Number of nonlinear solver iterations = %ld\n", numSolveIterations);
  printf("%ld Jacobian evaluations requiring %ld residual function evaluations.\n",
    numJacEval, numFuncEvalJac);
}