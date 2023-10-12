#include "wlp_fdtd2d.h"
#include <gauss_quadrature.h>

double function(double t, void *p){
  const double s = TIME_SCALE; 
  int n = *((int*)p);
  return jz(t) * wlp( n, s*t );
}

void calc_coef_jz(double *coef_jz, const int p_max){

  const double s = TIME_SCALE; 
  /* 積分する関数 */
  for(int n = 0; n <= p_max; n++){
    double ini_ans = AndoLab::gauss_quadrature( function, (void*)(&n), T_min_jz, T_max_jz );
    coef_jz[n] = s * AndoLab::converged_gauss_quadrature( function, (void*)(&n), T_min_jz, T_max_jz, ini_ans, 1.0e-8 );
  }
 
//  std::ifstream ifs;
//  std::string FREQ_TYPE;
//  if(FREQ > 1.0) FREQ_TYPE = "high"; else FREQ_TYPE = "low";
//  if(FREQ_TYPE == "high")ifs.open("../coef_jz/coef_jz_high_s2.dat");
//  if(FREQ_TYPE == "low")ifs.open("../coef_jz/coef_jz_low_s0.2.dat");
//  if(!ifs) {std::cout << "No coef_jz file." << std::endl; std::exit(0);}
//  for(int i = 0; i <= p_max; i++) ifs >> coef_jz[i];
//  ifs.close();
//  std::cout << "coef_jz[100] = " << coef_jz[100] << std::endl;

  return;
}  
