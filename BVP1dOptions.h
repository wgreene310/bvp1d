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

#ifndef BVP1dOptions_h
#define BVP1dOptions_h

class BVP1dOptions
{
public:
  BVP1dOptions(double absTol = 1e-6, double relTol=1e-3) :
    absTol(absTol), relTol(relTol) {
    vectorizedFuncs = stats = false;
    nMax = -1;
  }
  double getAbsTol() const { return absTol;  }
  void setAbsTol(double tol) { absTol = tol; }
  double getRelTol() const { return relTol; }
  void setRelTol(double tol) { relTol = tol; }
  bool isVectorized() const { return vectorizedFuncs; }
  void setVectorized(bool isVec) { vectorizedFuncs = isVec; }
  bool printStats() const { return stats; }
  void setPrintStats(bool sts) { stats = sts; }
  int getNMax() const { return nMax; }
  void setNMax(int nmax) { nMax = nmax; }
private:
  double absTol, relTol;
  bool vectorizedFuncs;
  bool stats;
  int nMax;
};

#endif

