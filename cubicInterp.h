
#ifndef cubicInterp_h
#define cubicInterp_h

#include <Eigen/Core>

template<class T, class T1>
void cubicInterp(const T &f1, const T &fp1, const T &f2, const T &fp2,
  double L, double s, T1 &fs, T1 &fps) {
  double sp1 = s + 1, sm1 = s - 1;;
  double ps2 = sp1*sp1;
  double ms2 = sm1*sm1;
  Eigen::Vector4d H, f, dHds;
  H << ms2*(2 + s) / 4, ms2*sp1 / 4, ps2*(2 - s) / 4, ps2*sm1 / 4;
  f << f1, fp1*L / 2, f2, fp2*L / 2;
  fs = H.transpose()*f;
  dHds << (3 * sm1*sp1) / 4, (sm1*(1 + 3 * s)) / 4, (-3 * sm1*sp1) / 4, (sp1*(3 * s - 1)) / 4;
  fps = dHds.transpose()*f;
}

#endif