#include "wlp_fdtd2d.h"


double jz(double t){

  const double f = FREQ;
  double sig = SRC_SIGMA;
  double t0 = T0_time;
  double jz_value = 0.0; // [A/m^2]
  
  //ガウシアンパルスx正弦波
  jz_value = AMP * std::exp( - std::pow((t - t0)/sig, 2.0) / 2.0 ) * std::sin(2.0*M_PI*f*(t-t0));

  return jz_value;
}
