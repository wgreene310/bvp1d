#pragma once

class BVPSolverStats
{
public:
  BVPSolverStats(bool recordStats);
  void update(void *kinsolMemory, int numMeshPts);
  void incrementJacobianFuncCalls(int numCalls);
  void print() const;
private:
  bool recordStats;
  long numFuncEvals, numSolveIterations, numJacEval;
  long numFuncEvalJac;
  int numPtsCurrentMesh;
  int numMeshUpdates;
};

