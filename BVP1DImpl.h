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

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::VectorXd RealVector;
typedef Eigen::MatrixXd RealMatrix;

class BVP1dOptions;

class BVP1DImpl {
public:
  class BVPDefn {
  public:
    virtual void odeFunc(double x, const RealVector &y,
      const RealVector &p, RealVector &f) = 0;
    virtual void bcFunc(const RealVector &y0, const RealVector &yn,
      const RealVector &p, RealVector &g) = 0;
  };
  BVP1DImpl(BVPDefn &bvp, RealVector &mesh, RealMatrix &yInit,
    RealVector &parameters, BVP1dOptions &options);
  ~BVP1DImpl();
  template<class T1, class T2>
  void calcPhi(const T1 &y, T2 &phi);
  template<class T>
  void calcError(const T &u, RealVector &err);
  int batheTest();
  int solve(RealMatrix &solMat, RealVector &paramVec);
  static RealVector linspace(double start, double end, int n);
  SparseMat J;
  RealVector rhs;
private:
  BVPDefn &bvp;
  RealVector &mesh;
  RealMatrix &yInit;
  RealVector parameters;
  const BVP1dOptions &options;
  int numNodes, numDepVars, numParams;
  RealVector gVec, fi, fim1, fim2, yim2;
};

