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

#ifndef cubicInterp_h
#define cubicInterp_h

#include <Eigen/Core>

template<class T1, class T2, class T3, class T4, class T5, class T6>
void cubicInterp(const T1 &f1, const T2 &fp1, const T3 &f2, const T4 &fp2,
  double L, double s, T5 &fs, T6 &dfDx) {
  double sp1 = s + 1, sm1 = s - 1;;
  double ps2 = sp1*sp1;
  double ms2 = sm1*sm1;
  Eigen::Vector4d H, dHds;
  H << ms2*(2 + s) / 4, ms2*sp1 / 4, ps2*(2 - s) / 4, ps2*sm1 / 4;
  dHds << (3 * sm1*sp1) / 4., (sm1*(1 + 3 * s)) / 4., (-3 * sm1*sp1) / 4., (sp1*(3 * s - 1)) / 4.;

  double L2 = L / 2.;
  fs = H[0] * f1 + H[1] * fp1*L2 + H[2] * f2 + H[3] * fp2*L2;
  dfDx =2/L*(dHds[0] * f1 + dHds[1] * fp1*L2 + dHds[2] * f2 + dHds[3] * fp2*L2);
}

#endif