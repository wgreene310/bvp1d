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

class BVP1DImpl {
public:
  class BVPDefn {
  public:
    virtual void odeFunc(double x, const RealVector &y,
      RealVector &f) = 0;
    virtual void bcFunc(const RealVector &y0, const RealVector &yn,
      RealVector &g) = 0;
  };
  BVP1DImpl(BVPDefn &bvp, RealVector &mesh, RealMatrix &yInit);
  ~BVP1DImpl();
  //Eigen::MatrixXd solve();
  template<class T>
  void calcPhi(const T &y, T &phi);
  int batheTest();
  int solve(Eigen::MatrixXd &solMat);
  static RealVector linspace(double start, double end, int n);
  SparseMat J;
  RealVector rhs;
private:
  BVPDefn &bvp;
  RealVector &mesh;
  RealMatrix &yInit;
  //RealVector phi;
  int numNodes, numDepVars;
  RealVector fi, fim1, fim2, yim2;
};

