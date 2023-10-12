#include "wlp_fdtd2d.h"

void calc_d(double *d, int N2,  int N3, int N4,  int N5, int N7,
            double d_max, double d_min){

    for(int i = 0; i < N7; i++){
        d[i] = 0.0;
    }
    
    for(int i=0; i < N2; i++){
        d[i] = d_max;
    }
    for(int i=N2; i < N3; i++){
        d[i] = d[i-1]/r;
    }
    for(int i=N3; i < N4; i++){
        d[i] = d_min;
    }
    for(int i=N4; i < N5; i++){
        d[i] = d[i-1]*r;
    }
    for(int i=N5; i < N7; i++){
        d[i] = d_max;
    }

    return;
}

void calc_d_h(double *d_h, double *d, int N2,  int N3, int N4,  int N5, int N7){
   
  for(int i = 0; i < N7; i++){
      d_h[i] = 0.0;
  }

  for(int i=1; i < N2; i++){
  	d_h[i] = ((d[i]+d[i-1])/2.0);
  }
  d_h[N2] = (d[N2-1]/2.0 + d[N2]*(1.0-beta));
  for(int i=N2+1; i < N3; i++){
  	d_h[i] = (d[i-1]*beta + d[i]*(1.0-beta));
  }
  d_h[N3] = (d[N3-1]*beta + d[N3]/2.0);
  for(int i=N3+1; i < N4; i++){
  	d_h[i] = ((d[i]+d[i-1])/2.0);
  }
  d_h[N4] = (d[N4-1]/2.0 + d[N4]*beta);
  for(int i=N4+1; i < N5; i++){
  	d_h[i] = (d[i-1]*(1.0-beta)+d[i]*beta);
  }
  d_h[N5] = (d[N5-1]*(1.0-beta) + d[N5]/2.0);
  for(int i=N5+1; i < N7; i++){
  	d_h[i] = ((d[i]+d[i-1])/2.0);
  }

    return;
}