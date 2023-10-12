#include "wlp_fdtd2d.h"

/***************************************************************************************
 *  Analysis Region, Eq.(20):         nd(i,j)(EZ) >= 0 && nd(i,j)(EZX(EZY)) == -1
 *  PML (Free Space), Eq.(59),(60):   nd(i,j)(EZ) == -1 && nd(i,j)(EZX(EZY)) >= 0
 *  PML (Ground), Eq.(88),(89),(90):  nd(i,j)(EZ) >= 0 && nd(i,j)(EZX(EZY)) >= 0
 *  Ez = Ezx + Ezy (except PML with conductive media)
 ***************************************************************************************/

void compose_coef_matrix(double *dx, double *dy, double *dx_h, double *dy_h,
					  int total_unknowns, int *n, css **S, csn **N, 
					  double *Dx, double *Dy, double *Cx, double *Cy,
                      double *const_hx1, double *const_hx2,
                      double *const_hy1, double *const_hy2,
                      double **conduct_media,
                      double *conduct_pmlx, double *conduct_pmly,
                      double ***node)
{ 

  constexpr double s { TIME_SCALE };
  int row = 0;
  double coef = 0.0;
  double const_d0 = 0.0;

  cs *T = cs_spalloc(0,0,1,1,1);

  for( int i = 1; i < NX7; i++ ){
    for( int j = 1; j < NY7; j++ ){
      const_d0 = 2.0*conduct_media[i][j] / EPS0/s;

      if( node[i][j][EZ] >= 0 && node[i][j][EZX] == -1 ){
			/*** Eq.(20) ******************************************/
				row = node[i][j][EZ];

				coef = (1.0 + const_d0 +
					Dx[i]*Cx[i] + Dx[i]*Cx[i-1] +
					Dy[j]*Cy[j] + Dy[j]*Cy[j-1] );
				  cs_entry(T, row, row, coef);

				coef = - Dx[i]*Cx[i];
				if( node[i+1][j][EZ] >= 0 ){
					cs_entry(T, row, node[i+1][j][EZ], coef );
				} //else if( node[i+1][j][EZX] >= 0 ){
					//cs_entry(T, row, node[i+1][j][EZX], coef );
					//cs_entry(T, row, node[i+1][j][EZY], coef );
				//}
				coef = - Dx[i]*Cx[i-1];
				if( node[i-1][j][EZ] >= 0 ){
					cs_entry(T, row, node[i-1][j][EZ], coef );
				} //else if( node[i-1][j][EZX] >= 0 ){
					//cs_entry(T, row, node[i-1][j][EZX], coef );
					//cs_entry(T, row, node[i-1][j][EZY], coef );
				//}

				coef = - Dy[j]*Cy[j];
				if( node[i][j+1][EZ] >= 0 ){
					cs_entry(T, row, node[i][j+1][EZ], coef );
				} //else if( node[i][j+1][EZX] >= 0 ){
					//cs_entry(T, row, node[i][j+1][EZX], coef );
					//cs_entry(T, row, node[i][j+1][EZY], coef );
				//}
				coef = -Dy[j]*Cy[j-1];
				if( node[i][j-1][EZ] >= 0 ){
					cs_entry(T, row, node[i][j-1][EZ], coef );
				} //else if( node[i][j-1][EZX] >= 0 ){
					//cs_entry(T, row, node[i][j-1][EZX], coef );
				  //cs_entry(T, row, node[i][j-1][EZY], coef );
				//}
			/*** Eq.(20) (up to here) *******************************/
      } else {
			/*** Eq.(88) ********************************************/
				row = node[i][j][EZX];

				coef = EPS0 * dx_h[i] * (s/2.0+conduct_pmlx[i]);//dx変更
				cs_entry(T, row, row, coef);

				coef = const_hy1[i] + const_hy1[i-1];
				cs_entry(T, row, node[i][j][EZ], coef);

				coef = - const_hy1[i];
				if( node[i+1][j][EZ] >= 0 ){
					cs_entry(T, row, node[i+1][j][EZ], coef);
				} //else if( node[i+1][j][EZX] >= 0 ){
	  			//cs_entry(T, row, node[i+1][j][EZX], coef);
	  			//cs_entry(T, row, node[i+1][j][EZY], coef);
				//}

				coef = - const_hy1[i-1];
				if( node[i-1][j][EZ] >= 0 ){
	  			cs_entry(T, row, node[i-1][j][EZ], coef);
				} //else if( node[i-1][j][EZX] >= 0 ){
	  			//cs_entry(T, row, node[i-1][j][EZX], coef);
	  			//cs_entry(T, row, node[i-1][j][EZY], coef);
				//}
			/*** Eq.(88) (up to here) *******************************/

			/*** Eq.(89) ********************************************/
				row = node[i][j][EZY];

				coef = EPS0 * dy_h[j] * (s/2.0+conduct_pmly[j]);//dy変更
				cs_entry(T, row, row, coef);

				coef = const_hx1[j] + const_hx1[j-1];
				cs_entry(T, row, node[i][j][EZ], coef);

				coef = - const_hx1[j];
				if( node[i][j+1][EZ] >= 0 ){
					cs_entry(T, row, node[i][j+1][EZ], coef);
				} //else if( node[i][j+1][EZX] >= 0 ){
					//cs_entry(T, row, node[i][j+1][EZX], coef);
					//cs_entry(T, row, node[i][j+1][EZY], coef);
				//}

				coef = - const_hx1[j-1];
				if( node[i][j-1][EZ] >= 0 ){
					cs_entry(T, row, node[i][j-1][EZ], coef);
				} //else if( node[i][j-1][EZX] >= 0 ){
					//cs_entry(T, row, node[i][j-1][EZX], coef);
					//cs_entry(T, row, node[i][j-1][EZY], coef);
				//}
			/*** Eq.(89) (up to here) *******************************/

			/*** Eq.(90) ********************************************/
				row = node[i][j][EZ];

				coef = s/2.0 + conduct_media[i][j]/EPS0;
				cs_entry(T, row, row, coef);

				coef = - s/2.0;
				cs_entry(T, row, node[i][j][EZX], coef);
				cs_entry(T, row, node[i][j][EZY], coef);
			/*** Eq.(90) (up to here) *******************************/
      }
    }
  }

  //LU decomp 
  cs *A = cs_triplet(T);
  cs_spfree(T);
  *n = A-> n;
  *S = cs_sqr(A, 1, 0);
  *N = cs_lu(A, (*S), 1e-14);

  return;
}
