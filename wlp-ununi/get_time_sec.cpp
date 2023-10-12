#include "wlp_fdtd2d.h"

double get_time_sec(void){
  nanoseconds t = steady_clock::now().time_since_epoch();
  return ( t.count()*1.0e-9 );
}
