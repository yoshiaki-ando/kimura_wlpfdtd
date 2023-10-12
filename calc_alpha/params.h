#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

const double C0 (2.99792458e8);
const double MU0 (4.0e-7 * M_PI);
const double EPS0 (1.0 / C0 / C0 / MU0);

const std::string file_ilp_pre ("input_data/ilp_pre.dat");
const std::string file_ilp_near ("input_data/ilp_near.dat");

// Spitak
/*
const std::string file_dir ("input_data/spitak/ft_bx_div_Iz.dat");
constexpr double first_freq (0.1);
constexpr double last_freq (1.0);
constexpr double obs_value (0.2e-9);
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "spitak";
const std::string eq_mag = "6.75";
const std::string eq_type = "near";
*/

// Biak
/*
const std::string file_dir ("input_data/biak/ft_by_div_Iz.dat");
// constexpr double first_freq (0.005);
// constexpr double last_freq (0.03);
constexpr double first_freq (0.005);
constexpr double last_freq (0.005);
constexpr double obs_value (0.16e-9); // 垂直
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "biak";
const std::string eq_mag = "8.19";
const std::string eq_type = "pre";
*/

// Alum Rock
/*
const std::string file_dir ("input_data/alumrock/ft_bx_div_Iz.dat");
constexpr double first_freq (0.001);
constexpr double last_freq (12);
constexpr double obs_value (1.4e-9); // 水平
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "alumrock";
const std::string eq_mag = "5.56";
const std::string eq_type = "near";
*/

// Guam 
/*
const std::string file_dir ("input_data/guam/ft_by_div_Iz.dat");
constexpr double first_freq (0.01);
constexpr double last_freq (0.01);
constexpr double obs_value (0.1e-9); // By
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "guam";
const std::string eq_mag = "7.75";
const std::string eq_type = "pre";
*/

// laquila
/*
const std::string file_dir ("input_data/laquila/ft_by_div_Iz.dat");
constexpr double first_freq (0.01);
constexpr double last_freq (0.015);
constexpr double obs_value (40.e-12); // 0.01~0.015Hz: 35~40pT(By)
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "laquila";
const std::string eq_mag = "6.32";
const std::string eq_type = "pre";
*/

// kagoshima
/*
const std::string file_dir ("input_data/kagoshima/ft_by_div_Iz.dat");
constexpr double first_freq (0.007);
constexpr double last_freq (0.013);
constexpr double obs_value (1.e-12); // 0.01+-0.003Hz: 1pT(By)
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "kagoshima";
const std::string eq_mag = "6.1";
const std::string eq_type = "pre";
*/

// Izu swarm
/*
const std::string file_dir ("input_data/izu/ft_ez_div_Iz.dat");
constexpr double first_freq (0.01);
constexpr double last_freq (0.01);
constexpr double obs_value (2.e-6); // 0.01Hz: 1mV/km = 1uV/m
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "izu";
const std::string eq_mag = "6.11";
const std::string eq_type = "pre";
*/

// iwate_nairiku_hokubu
/*
const std::string file_dir ("input_data/iwate_nairiku_hokubu/ft_by_div_Iz.dat");
constexpr double first_freq (0.007);
constexpr double last_freq (0.013);
constexpr double obs_value (0.3e-12); // By 0.00030nT/√Hz
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "iwate_nairiku_hokubu";
const std::string eq_mag = "5.82";
const std::string eq_type = "pre";
*/

// wenchuan
/*
const std::string file_dir ("input_data/wenchuan/ft_ez_div_Iz.dat");
constexpr double first_freq (0.1);
constexpr double last_freq (10.0);
constexpr double obs_value (0.11e-3); // Ez: 1.3mV/m
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "wenchuan";
const std::string eq_mag = "7.91";
const std::string eq_type = "near";
*/

// jammu

const std::string file_dir ("input_data/jammu_kashmir/ft_by_div_Iz.dat");
constexpr double first_freq { 0.01 };
constexpr double last_freq { 0.1 };
constexpr double obs_value { 2.0e-9 }; // By: 2nT
constexpr bool fl_per_band = true; //帯域割フラグ true:する, false:しない
const std::string eq_name = "jammu";
const std::string eq_mag = "5.18";
const std::string eq_type = "near";
