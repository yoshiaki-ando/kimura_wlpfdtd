#include "wlp_fdtd2d.h"


void make_profile(double time_pre, double time_fdtd, double time_post, double time_total, std::string file_name){
  
  std::ofstream ofs_prof( file_name.c_str() );
  ofs_prof
  << "Computational Region (Rx X Ry) = " << RX/1.0e3 << " km X " << RY/1.0e3 << " km" << std::endl
  << "Cell Size (Dx X Dy) = " << DX/1.0e3 << " km X " << DY/1.0e3 << " km" << std::endl
  << "Number of Cells (Nx, Ny) = (" << NX << ", " << NY << ")" << std::endl
  << "Analysis Time (RT) = " << RT << " s" << std::endl
  << "Time Step (for output) (DT) = " << DT << " s" << std::endl
  << "Number of Time Steps (NT) = " << NT << std::endl 
  << std::endl
  << "Total Computational Region (Rtx X Rty) = " << (LX2*2+RX)/1e3 << " km X " << (LY2*2+RY)/1e3 << " km" << std::endl
  << "Number of Buffer Cells (NCX, NCY) = (" << NCX << ", "<<  NCY << ")" << std::endl
  << "Number of Ununiform Cells (NMX, NMY) = (" << NMX << ", " << NMY << ")" << std::endl
  << "Increase Rate of Cell Size = " << r << std::endl
  << "Position Rate of H in Ununiform Cell(<=0.5) = " << beta << std::endl
  << "Maximum Cell Size: dx_max = " << dx_max/1.0e3 << ", dy_max = " << dy_max/1.0e3 << " [km]" << std::endl
  << "Total Length of Ununiform Cells: " << "LX1 = " << LX1/1.0e3 << ", " << "LY1 = " << LY1/1.0e3 << " [km]" << std::endl
  << "Length of Extended Area Including PML(One side): " << "LX2 = " << LX2/1.0e3 << ", " << "LY2 = " << LY2/1.0e3 << " [km]" << std::endl
  << std::endl
  << "Origin Point (x0, y0) = (" << X_BOUNDARY/1.0e3 << " km, " << Y_SURFACE/1.0e3 << " km)" << std::endl
  << "Observatory Point (x_o-x0, y_o-y0) = (" << (X_OBSERVE - X_BOUNDARY)/1.0e3 << " km, " << (Y_OBSERVE - Y_SURFACE)/1.0e3 << " km)" << std::endl 
  << "Source Point (x_s-x_0, y_s-y0) = (" << (X_SOURCE - X_BOUNDARY)/1.0e3 << " km, " << (Y_SOURCE - Y_SURFACE)/1.0e3 << " km)" << std::endl 
  << std::endl
  << "P_ORDER = " << P_ORDER << std::endl
  << "TIME_SCALE = " << TIME_SCALE << std::endl
  << std::endl
  << "PML_LAYER: " << L_PML << std::endl
  << "*Jz*" << std::endl
  << "Interval of integration(T_min_jz ~ T_max_jz) = " << T_min_jz << " s ~ " << T_max_jz << " s"<< std::endl
  << std::endl
  << "time_pre = " << time_pre << " s" << std::endl
  << "time_fdtd = " << time_fdtd << " s" << std::endl
  << "time_post = "<< time_post << " s" << std::endl
  << "time_total = " << time_total << " s" << std::endl;
  ofs_prof.close();
  
  return;
}
