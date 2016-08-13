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

