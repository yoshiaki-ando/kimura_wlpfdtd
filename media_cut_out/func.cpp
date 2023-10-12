#include "cut_out.h"

coordinate_sph_orth coordinate_sph_orth::operator+ (coordinate_sph_orth a){
  coordinate_sph_orth b(x + a.x, y + a.y, z + a.z);
  return b;
}
coordinate_sph_orth coordinate_sph_orth::operator- (coordinate_sph_orth a){
  coordinate_sph_orth b(x - a.x, y - a.y, z - a.z);
  return b;
}

coordinate_sph_orth::coordinate_sph_orth(void){
  x = 0.0; y = 0.0; z = 0.0; lat = 0.0; lon = 0.0;
}
coordinate_sph_orth::coordinate_sph_orth(double given_lat, double given_lon){
  lat = given_lat; 
  lon = given_lon;
  x = std::cos(lat)*std::cos(lon);
  y = std::cos(lat)*std::sin(lon);
  z = std::sin(lat);
}
coordinate_sph_orth::coordinate_sph_orth(double a, double b, double c){
  x = a; y = b; z = c;
}

coordinate_sph_orth::coordinate_sph_orth(double given_lat, double given_lon, double a, double b, double c){
  lat = given_lat; lon = given_lon;
  x = a; y = b; z = c;
}

void coordinate_sph_orth::sph_coordinate(void){
  lat = std::asin(z);
  if(y >= 0)
    lon = std::acos(x/std::cos(lat)); // 0~pi
  else if(y < 0)
    lon = -std::acos(x/std::cos(lat)); // pi~2pi
}
void coordinate_sph_orth::orth_coordinate(void){
  x = std::cos(lat)*std::cos(lon);
  y = std::cos(lat)*std::sin(lon);
  z = std::sin(lat);
}

double coordinate_sph_orth::abs(void) {
  return std::sqrt(x*x + y*y + z*z);
}

void coordinate_sph_orth::print_params(void){
 std::cout << "(x, y, z) = (" << x << ", " << y << ", " << z << "), (lat, lon) = (" << lat/M_PI*180.0 << ", " << lon/M_PI*180.0 << ")[°]" << std::endl;
}

double inner_product(coordinate_sph_orth a, coordinate_sph_orth b){
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

double arg(coordinate_sph_orth a, coordinate_sph_orth b){
  return std::acos(inner_product(a,b)/a.abs()/b.abs());
}

double dist(coordinate_sph_orth a, coordinate_sph_orth b){
  const double r = R_EARTH;
  return r*acos(std::sin(a.lat)*std::sin(b.lat) + std::cos(a.lat)*std::cos(b.lat)*std::cos(b.lon-a.lon));
}

double calc_quadratic_equation(double a, double b, double c, const int i, const bool fl) {
  double d = (-b-std::sqrt(b*b-4.0*a*c)) / (2.0*a); 
  double e = (-b+std::sqrt(b*b-4.0*a*c)) / (2.0*a);
  //std::cout << i << ": "<< d << " " << e << std::endl;
  if(!fl) return d;
  else return e;
}

void divide_param(coordinate_sph_orth a, coordinate_sph_orth b, double* t, const int nt){
  // r(t) = r_s + t(r_e - r_s)の直線上を球面上で等分するparam:tを計算
  
  const double dangle = arg(a, b)/(dist(a,b)/dr);
  // 0 < t < 1の外側の範囲は最大で(pi-arg(r_s, r_e)/2.0)までしか計算できない
  if( (Rr-dist(a,b))/2./dr*dangle > (M_PI-arg(a,b))/2.0 ) {
    std::cout << "[Error] You can not calculate outside of arg(a,b) +-(pi-arg(a,b))/2[rad]." << std::endl;
    std::cout << "Please let cut length(Rr) smaller than now (but Rr > dist(a,b))" << std::endl; 
    std::exit(0);
  }

  bool fl=false; // aとbの偏角が90度を超えた: true
  coordinate_sph_orth c = b-a;
  double A = inner_product(a, c);
  for(int i = 0; i <= nt; i++){
    double angle = (i - (Rr-dist(a,b))/2./dr )* dangle;
    std::cout << angle << std::endl;
    double d = A*A-std::pow(c.abs()*a.abs()*std::cos(angle),2);
    double e = 2.0*A*( std::pow(a.abs(),2) * (1.0-std::pow(std::cos(angle),2)) );
    double f = std::pow(a.abs(),4)*(1.0-std::pow(std::cos(angle), 2));
    if(angle > M_PI/2. || angle < 0.0) {
      fl = true;
    }else{
      fl = false;
    }
    t[i] = calc_quadratic_equation(d, e, f, i, fl); 
  }
  return;
}

void coordinate_sph_orth::vec_line_on_sph(coordinate_sph_orth a, coordinate_sph_orth b, double t) {
  // r(t) = r_s + t(r_e - r_s)
  x = a.x + t*(b.x-a.x);
  y = a.y + t*(b.y-a.y);
  z = a.z + t*(b.z-a.z);

  double magnitude = std::sqrt(x*x+y*y+z*z);
  x/=magnitude;
  y/=magnitude;
  z/=magnitude;

  return;
}

std::vector<double> cross_product(std::vector<double> v1, std::vector<double> v2) {
  std::vector<double> ret(3);
  double a, b, c;
  a = v1[1]*v2[2] - v2[1]*v1[2];
  b = v1[2]*v2[0] - v2[2]*v1[0];
  c = v1[0]*v2[1] - v2[0]*v1[1];
  ret = {a, b, c};
  return ret;
}

std::vector<double> calc_coef_plane(coordinate_sph_orth v1, coordinate_sph_orth v2, coordinate_sph_orth v3, std::vector<double> sig) {
  std::vector<double> coef(3);
  std::vector<double> r1, r2(3);
  r1 = {v2.lon-v1.lon, v2.lat-v1.lat, sig[1]-sig[0]}; 
  r2 = {v3.lon-v1.lon, v3.lat-v1.lat, sig[2]-sig[0]}; 
  double d;
  coef = cross_product(r1, r2);
  d = -(coef[0]*v1.lon + coef[1]*v1.lat + coef[2]*sig[0]);
  coef.resize(4);
  coef[3] = d;
  return coef;
}

double linear_interpolaton(std::vector<double> coef, coordinate_sph_orth v) {
  return -(coef[0]*v.lon + coef[1]*v.lat + coef[3])/coef[2];
}
