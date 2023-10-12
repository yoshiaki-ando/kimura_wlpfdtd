#include "params.h"

const double m_eps = 1e-10;

class coordinate_sph_orth{
  public:
    coordinate_sph_orth operator+ (coordinate_sph_orth a);
    coordinate_sph_orth operator- (coordinate_sph_orth a);
    coordinate_sph_orth(void);
    coordinate_sph_orth(double given_lat, double given_lon);
    coordinate_sph_orth(double a, double b, double c);
    coordinate_sph_orth(double given_lat, double given_lon, double a, double b, double c);
    void sph_coordinate(void);
    void orth_coordinate(void);
    double abs(void);
    void print_params();
    void vec_line_on_sph(coordinate_sph_orth a, coordinate_sph_orth b, double t);
    double lat; /* Latitude [rad] */
    double lon; /* Longitude [rad] */
    /* 半径1の球面上の直交座標, 基準は球の中心 */
    double x;
    double y;
    double z;
}; 

double inner_product(coordinate_sph_orth a, coordinate_sph_orth b);
double arg(coordinate_sph_orth a, coordinate_sph_orth b);
double dist(coordinate_sph_orth a, coordinate_sph_orth b);
double calc_quadratic_equation(double a, double b, double c, const int i, const bool fl);
void divide_param(coordinate_sph_orth a, coordinate_sph_orth b, double* t, const int nt);
std::vector<double> cross_product(std::vector<double> v1, std::vector<double> v2);
std::vector<double> calc_coef_plane(coordinate_sph_orth v1, coordinate_sph_orth v2, coordinate_sph_orth v3, std::vector<double> sig);
double linear_interpolaton(std::vector<double> coef, coordinate_sph_orth v);
