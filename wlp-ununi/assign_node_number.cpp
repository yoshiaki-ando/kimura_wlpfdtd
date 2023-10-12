#include "wlp_fdtd2d.h"

void assign_node_number(double ***node, int *total_unknowns){

  const int l = L_PML;
  int node_number = 0;

  
  // Boundary (not update)
  for( int i = 0; i <= NX7; i++ ){
    for(int n = 0; n < 3; n++){
      node[i][0][n] = -2;
      node[i][NY7][n] = -2;
    }
  }
  for( int j = 0; j <= NY7; j++ ){
    for(int n = 0; n < 3; n++){
      node[0][j][n] = -2;
      node[NX7][j][n] = -2;
    }
  }

  for( int i = 1; i < NX7; i++ ){
    for( int j = 1; j < NY7; j++ ){
      if( (i <= l) || (NX7-l <= i) ||
  	  (j <= l) || (NY7-l <= j) ){ // PML
        node[i][j][EZ] = node_number++;
        node[i][j][EZX] = node_number++;
        node[i][j][EZY] = node_number++;
      } else { // Analysis Region
        node[i][j][EZ] = node_number++;
      }
    }
  }
  
  (*total_unknowns) = node_number;

  
  return;
}
