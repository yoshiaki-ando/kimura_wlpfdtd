#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <ctime>
#include <vector>
#include <complex>

#include "../csparse.h"

constexpr double C0 (2.99792458e8);
constexpr double MU0 (4.0e-7 * M_PI);
constexpr double EPS0 (1.0 / C0 / C0 / MU0);
constexpr double Z0 (std::sqrt(MU0/EPS0));
constexpr double Y0 (1.0/Z0);

constexpr double AMP (1.0);
constexpr double FREQ (0.05);
constexpr double PERIOD (1.0/FREQ);

constexpr double DT { 1.0/100 }; //1/fs[s^-1]
constexpr double DF { 0.0002 }; //基本周波数(逆数がRT)

constexpr int NT { (int)std::round(double(1.0/DF/DT)) };
constexpr double RT { NT*DT };

//波源 
constexpr double SRC_SIGMA { 0.45*PERIOD };
constexpr double T0_time { 7.0*SRC_SIGMA };

/* Jzの係数導出用 */
constexpr double T_max_jz {T0_time + 6.0*SRC_SIGMA};
constexpr double T_min_jz {T0_time - 6.0*SRC_SIGMA};

#define EZ 0
#define EZX 1
#define EZY 2
#define HX 3
#define HY 4
#define TOTAL 5

#define LOC_DATA ( std::string("data/"))
#define LOC_PLOT ( std::string("plot/"))

/****************Biak****************/
constexpr double RX {150.0e3}; 
constexpr double RY {100.0e3}; 
constexpr double DX {RX/1500};
constexpr double DY {RY/1000};
constexpr int NX { int(RX/DX) };
constexpr int NY { int(RY/DY) };

constexpr double X_BOUNDARY {0.5*RX};
constexpr double Y_SURFACE {0.6*RY};
constexpr double Distance_btw_Src_Obs { 105.0e3 };
constexpr double X_SOURCE {(RX/2.0) - Distance_btw_Src_Obs/2.0};
constexpr double Y_SOURCE {Y_SURFACE - 20.0e3}; //ISC
constexpr double X_OBSERVE {(RX/2.0) + Distance_btw_Src_Obs/2.0};
constexpr double Y_OBSERVE {Y_SURFACE};
/*****************************************/

// 0.001-0.1Hz
constexpr int P_ORDER (1000);
constexpr double TIME_SCALE ( 0.2 );

// 0.1-10Hz
// constexpr int P_ORDER (5000);
// constexpr double TIME_SCALE ( 2.0 );

constexpr double TIME_MARGIN (0.0);

constexpr int L_PML (12); 
constexpr double M_PML (4.0);
constexpr double R_PML (1.0e-6);

/***************追加**********************/
constexpr double Rx_total { 5.e9 }; // 全解析領域の全長

// 増加率から不均一セル数を求める
constexpr double dx_max { 1.e8 };  // 最大セルサイズ(x方向)
constexpr double r_target { 1.6 }; // 不均一セルの増加率
constexpr int NMX { int(std::ceil(std::log(dx_max/DX)/std::log(r_target))-1) }; // 不均一セルの個数(片側, x方向)
constexpr double r { std::exp( std::log(dx_max/DX)/(NMX+1) ) }; // 不均一セルの増加比
constexpr int NMY { int(std::log(dx_max/DY)/std::log(r)-1) }; // 不均一セルの個数(片側, y方向)
constexpr double dy_max { std::pow(r,NMY+1)*DY };  // 最大セルサイズ(y方向)

constexpr double LX1 { (DX*r - dx_max)/(1-r) };  // 不均一セルの全長 (x方向)
constexpr double LY1 { (DY*r - std::pow(r,NMY+1)*DY)/(1-r) };  // 不均一セルの全長 (y方向)

constexpr int NCX { int(std::ceil( (Rx_total-2*LX1)/2/dx_max - L_PML )) };  // 不均一とPML間の緩衝セル数(切り上げ)
constexpr int NCY { int(std::ceil( (Rx_total-2*LY1)/2/dy_max - L_PML )) };  // 不均一とPML間の緩衝セル数(切り上げ)
constexpr double beta { 1.0/(1.0+std::sqrt(r)) };

constexpr int NX1 { L_PML };       //PMLのセル数
constexpr int NX2 { NX1 + NCX };    //左端から左側不均一セルの始まりまでのセル数
constexpr int NX3 { NX2 + NMX };    //左端から均一セルのの始まりまでのセル数
constexpr int NX4 { NX3 + NX };    //左端から均一セルのルの終わりまでのセル数 
constexpr int NX5 { NX4 + NMX };    //左端から右側不均一セルの終わりまでのセル数
constexpr int NX6 { NX5 + NCX };    //左端から右側PMLの始まりまでのセル数
constexpr int NX7 { NX6 + NX1};    //X方向の全セル数

constexpr int NY1 { L_PML };       
constexpr int NY2 { NY1 + NCY };       
constexpr int NY3 { NY2 + NMY };      
constexpr int NY4 { NY3 + NY };     
constexpr int NY5 { NY4 + NMY };       
constexpr int NY6 { NY5 + NCY };     
constexpr int NY7 { NY6 + NY1};   

constexpr double LX2 { LX1 + (NX2)*(dx_max) };  //左端から解析領域セルの左端までの距離
constexpr double LY2 { LY1 + (NY2)*(dy_max) };  //下端から解析領域セルの下端までの距離
/*******************************************/


using namespace std::chrono;
std::string get_date_str(void);
double get_time_sec(void);
void display_progress(std::string process, int prog, int prog_max);

double *allocate_memory_d1d(const int N_of_Array, const double Initial_value);
double **allocate_memory_d2d(const int M, const int N, const double Initial_value);
void free_memory_d2d(double **mat);
double ***allocate_memory_d3d(int ir, int jr, int kr, double Initial_value);
void free_memory_d3d(double ***mat);
void assign_node_number(double ***node, int *total_unknowns);

void compose_coef_matrix(double *dx, double *dy, double *dx_h, double *dy_h,
                      int total_unknowns, int *n, css **S, csn **N, 
					            double *Dx, double *Dy, double *Cx, double *Cy,
                      double *const_hx1, double *const_hx2,
                      double *const_hy1, double *const_hy2,
                      double **conduct_media,
                      double *conduct_pmlx, double *conduct_pmly,
                      double ***node); 
void calc_coef_jz(double *coef, const int pmax);
void calc_conductivity(double *x_pos, double *ypos,
                       double **conduct_media,
		               double *conduct_pmlx_ez, double *conduct_pmly_ez,
		               double *conduct_pmlx_h, double *conduct_pmly_h
                       );
void calc_const_param(double *Dx, double *Dy, double *Cx, double *Cy,
                      double *dx, double *dy, double *dx_h, double *dy_h);
void calc_const_pml(double *dx, double *dy,
                    double *const_hx1, double *const_hx2,
		                double *const_hy1, double *const_hy2,
                    double *conduct_pmlx_ez, double *conduct_pmly_ez,
                    double *conduct_pmlx_h, double *conduct_pmly_h
                    );
void calc_wlp(double **wlp, const int p_max, const int n_max, const double Dt);

double jz(double t);
void compose_right_side(double *dx, double *dy, double *dx_h, double *dy_h,
                       double *b,
                       double *coef_jz,
                       double **sum_coef_ez,
                       double **sum_coef_ezx, double **sum_coef_ezy,
                       double **sum_coef_hx, double **sum_coef_hy,
                       double *Dx, double *Dy,
                       double *const_hx1, double *const_hx2,
                       double *const_hy1, double *const_hy2,
                       double ***node, int p);
void update_e(double **coef_ez,
              double **sum_coef_ez,
              double **sum_coef_ezx, double **sum_coef_ezy,
              double *x,
              double ***node,
              int p);
void update_h(double *dx, double *dy,
              double **coef_hx, double **coef_hy,
              double **sum_coef_hx, double **sum_coef_hy,
              double **coef_ez,
              double *Cx, double *Cy,
              double *const_hx1, double *const_hx2,
              double *const_hy1, double *const_hy2,
              double ***node,
              int p);
void expand_field_p(double coef_field, double *field, 
                    double *wlp,
                    int p);
void expand_field_2d_p(double **coef_field, double ***field,
                           int nx, int ny, int nt,
                           double *wlp,
                           int p);
void output_field_1d(double* field, int forN, double delta, 
                     std::string filename);
void output_field_2d(double ***field, int nx, int ny, int nt, double delta, int n_thin_out,
                     std::string field_type, std::string directory_of_data);
void output_conductivity(double *x_pos, double *y_pos,
                         double **condct_media, 
                         std::string directory_of_data);
void make_profile(double time_pre, double time_fdtd, double time_post, 
                  double time_total, std::string file_name);
void fourier_transform(double *field, std::complex <double> *FTdata);
void fast_fourier_transform(double *field, std::complex <double> *FTdata, int nt = NT, double dt = DT);
void calc_transfer_function(std::complex<double> *A, 
                            std::complex<double> *B, 
                            double *T);
void output_spectrum(std::complex <double> *FTdata, 
                    int forN, double delta, 
                    std::string file_name);

void calc_d(double *d, int N2,  int N3, int N4,  int N5, int N7,
            double d_max, double d_min);
void calc_d_h(double *d_h, double *d, int N2,  int N3, int N4,  int N5, int N7);
void set_coordinate(double *dx, const int nx, double *x, const double x_orig,
                    double *dy, const int ny, double *y, const double y_orig);

double weighting(const double x, const double w, const double t);
double wlp(const int p, const double t);
void calc_wlp_for_output(const int p, double **w, double **v, double *ret, const double dt, const int max_n);
