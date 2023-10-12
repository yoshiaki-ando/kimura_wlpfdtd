#include "wlp_fdtd2d.h"


std::string get_date_str(void){
  std::time_t t = system_clock::to_time_t( system_clock::now() );
  std::tm *date = std::localtime(&t);
    
  std::ostringstream oss_date;
  oss_date
  << std::setw(2) << std::setfill('0') << date->tm_year - 100
	   << std::setw(2) << std::setfill('0') << date->tm_mon + 1
	   << std::setw(2) << std::setfill('0') << date->tm_mday
	   << "_"
	   << std::setw(2) << std::setfill('0') << date->tm_hour
	   << std::setw(2) << std::setfill('0') << date->tm_min;
    
  return oss_date.str();
}
