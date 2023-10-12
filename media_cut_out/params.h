constexpr double R_EARTH = 6378.137e3;
#define deg_to_rad(deg) (((deg)/360.0)*2.0*M_PI)
#define rad_to_deg(rad) (((rad)/2.0/M_PI)*360.0)

/*
// Spitak EQ(1988)
const std::string eq_name {"spitak"};
constexpr double dr = 1e2; // 導電率の取得間隔
constexpr double Rr = 150e3; // 切り取る全長
constexpr double lat_epi = 41.0; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 44.20; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 42.1;
constexpr double lon_obs = 44.68;
*/

// Biak EQ(1996)
/*
const std::string eq_name {"biak"};
constexpr double dr = 1e2;
constexpr double Rr = 150e3;
constexpr double lat_epi = -0.27; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 136.54; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = -1.08;
constexpr double lon_obs = 136.05;
*/

// Guam EQ(1993)
/*
const std::string eq_name {"guam"};
constexpr double dr = 1e2;
constexpr double Rr = 100e3;
constexpr double lat_epi = 13.0; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 144.75; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 13.58;
constexpr double lon_obs = 144.87;
*/

// Alum Rock EQ(2007)
/*
const std::string eq_name {"alumrock"};
constexpr double dr = 1e2;
constexpr double Rr = 100e3;
constexpr double lat_epi = 37.432; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = -121.776; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 37.415;
constexpr double lon_obs = -121.780;
*/
/*
// kagoshima EQ(2007)
const std::string eq_name {"kagoshima"};
constexpr double dr = 1e2;
constexpr double Rr = 150e3;
constexpr double lat_epi = 32.0; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 130.3; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 31.48;
constexpr double lon_obs = 130.72;
*/

/*
// L'Aquila EQ(2009)
const std::string eq_name {"laquila"};
constexpr double dr = 1e2;
constexpr double Rr = 100e3;
constexpr double lat_epi = 42.368; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 13.353; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 42.23;
constexpr double lon_obs = 13.19;
*/

// Izu swarm EQ(2000)

const std::string eq_name {"izu"};
constexpr double dr = 1e2;
constexpr double Rr = 100e3;
constexpr double lat_epi = 34.229 ; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 139.272; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 34.421;
constexpr double lon_obs = 139.283;


// Iwate nairiku hokubu EQ(1998)
/*
const std::string eq_name {"iwate_nairiku_hokubu"};
constexpr double dr = 1e2;
constexpr double Rr = 100e3;
constexpr double lat_epi = 39.8; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 140.864; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 39.88;
constexpr double lon_obs = 140.94;
*/

// Wenchuan
/*
const std::string eq_name {"wenchuan"};
constexpr double dr = 5e3;
constexpr double Rr = 2000e3;
constexpr double lat_epi = 31.02; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 103.37; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 39.3;
constexpr double lon_obs = 115.9;
*/

// Jammu kashmir
/*
const std::string eq_name {"jammu_kashmir"};
constexpr double dr = 5e3;
constexpr double Rr = 1500e3;
constexpr double lat_epi = 35.32; // 90 = 北緯90度 -90 = 南緯90度
constexpr double lon_epi = 74.83; // -180 = 西経180度 180° = 東経180度
constexpr double lat_obs = 27.20;
constexpr double lon_obs = 78.0;
*/

constexpr double lat_epi_rad = deg_to_rad(lat_epi);
constexpr double lon_epi_rad = deg_to_rad(lon_epi);
constexpr double lat_obs_rad = deg_to_rad(lat_obs);
constexpr double lon_obs_rad = deg_to_rad(lon_obs);
