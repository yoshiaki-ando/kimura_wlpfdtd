#include "wlp_fdtd2d.h"

/***************************************************************************************
 *  Analysis Region, Eq.(20):         nd[i][j][EZ] >= 0 && nd[i][j](EZX[EZY] == -1
 *  PML (Free Space), Eq.(59),(60):   nd[i][j][EZ] == -1 && nd[i][j](EZX[EZY] >= 0
 *  PML (Ground), Eq.(88),(89),(90):  nd[i][j][EZ] >= 0 && nd[i][j](EZX[EZY] >= 0
 *  Ez = Ezx + Ezy (except PML with conductive media)
 ***************************************************************************************/


void compose_right_side(double *dx, double *dy, double *dx_h, double *dy_h,
                       double *b,
                       double *coef_jz,
                       double **sum_coef_ez,
                       double **sum_coef_ezx, double **sum_coef_ezy,
                       double **sum_coef_hx, double **sum_coef_hy,
                       double *Dx, double *Dy,
                       double *const_hx1, double *const_hx2,
                       double *const_hy1, double *const_hy2,
                       double ***node, int p)
{
  constexpr int i_jz { int(X_SOURCE/DX + 0.5 + NX3) };
  constexpr int j_jz { int(Y_SOURCE/DY + 0.5 + NY3) };
  constexpr double s { TIME_SCALE };
  int row = 0;

  for( int i = 1; i < NX7; i++ ){
    for( int j = 1; j < NY7; j++ ){

      if( node[i][j][EZ] >= 0 && node[i][j][EZX] == -1 ){
				// Eq.(20)
				row = node[i][j][EZ];
				b[row] = ( - 2.0 * Dx[i] * ( sum_coef_hy[i][j]-sum_coef_hy[i-1][j] )
								   + 2.0 * Dy[j] * ( sum_coef_hx[i][j]-sum_coef_hx[i][j-1] )
								   - 2.0 * sum_coef_ez[i][j] ); // coef_jz is considered outside the loop.
      } else {
         
				// Eq.(88)
				row = node[i][j][EZX];
				b[row] = ( - EPS0* dx_h[i] *s * sum_coef_ezx[i][j]//dx変更
								- const_hy2[i] * sum_coef_hy[i][j]
								+ const_hy2[i-1] * sum_coef_hy[i-1][j] );

				// Eq.(89)
				row = node[i][j][EZY];
				b[row] = ( - EPS0* dy_h[j] *s * sum_coef_ezy[i][j]//dy変更
								+ const_hx2[j] * sum_coef_hx[i][j]
								- const_hx2[j-1] * sum_coef_hx[i][j-1] );
				// Eq.(90)
				row = node[i][j][EZ];
				b[row] = s * ( - sum_coef_ez[i][j]
							+ sum_coef_ezx[i][j]
							+ sum_coef_ezy[i][j] );
      }
    }
  }

  row = node[i_jz][j_jz][EZ];
  b[row] += - 2.0/EPS0/s * coef_jz[p];

  return;
}
