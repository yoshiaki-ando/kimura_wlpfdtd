#include "wlp_fdtd2d.h"


void display_progress(std::string process, int prog, int prog_max){

  const int prog_step = 20;
  int pct = 100 * prog/prog_max;
  static int pct_flag = -1;
  
  if( pct % prog_step == 0 && pct != pct_flag ){
    pct_flag = pct;
    if(pct == 0)
      std::cout << process << "... "<< std::endl;
    std::cout << pct << " %. " << "(" << prog << "/" << prog_max << ")"<< std::endl;
  }
  
  return;
}
