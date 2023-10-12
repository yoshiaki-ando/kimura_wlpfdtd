#include "wlp_fdtd2d.h"

void output_conductivity(double *x_pos, double *y_pos,
                         double **condct_media,
                         std::string directory_of_data){

  std::ofstream ofs_condct(directory_of_data + "/conductivity.dat");
  for( int i = 0; i <= NX7; i++ ){
    for( int j = 0; j <= NY7; j ++ ){
      ofs_condct << x_pos[i] << " " << y_pos[j] << " " << condct_media[i][j] << std::endl;
    }
    ofs_condct << std::endl;
  }
  ofs_condct.close();


  //output gnuplot script 
  std::string unit = "[m]";
  double unit_scale = 1.0;
  if( RX/1.0e3 >= 1.0 ){
    unit = "[km]";
    unit_scale = 1.0e3;
  }
  
  std::ofstream ofs_script(directory_of_data + "/plot_conductivity.gpt");
  ofs_script
  << "reset" << std::endl
  << "set term png size 1000, 600" << std::endl
  << "set output \""
  << directory_of_data + "/conductivity.png\"" << std::endl
  << std::endl
  << "unset key" << std::endl
  << "set tics font \"Arial, 16\"" << std::endl
  << "set size ratio 0.75" << std::endl
  << "set pm3d map" << std::endl
  << "set palette rgbformulae 30, 31, 32" << std::endl
  << "set grid" << std::endl
  << std::endl
  << "set xlabel \"x " << unit << "\"" << std::endl
  << "set xlabel font \"Arial, 16\"" << std::endl
  << "set ylabel \"y " << unit << "\"" << std::endl
  << "set ylabel font \"Arial, 16\"" << std::endl
  << "set ylabel offset -2, 0" << std::endl
  << "set cblabel \"Conductivity [S/m]\"" << std::endl
  << "set cblabel font \"Arial, 16\"" << std::endl
  << "set cblabel offset 2, 0" << std::endl
  << std::endl
  << "set xrange [" << -X_BOUNDARY/unit_scale << ":"
  << (RX-X_BOUNDARY)/unit_scale << "]" << std::endl
  << "set yrange [" << -Y_SURFACE/unit_scale << ":"
  << (RY-Y_SURFACE)/unit_scale << "]" << std::endl
  << "set cbrange [1.0e-5:10]" << std::endl
  << "set logscale cb" << std::endl
  << "set format cb \"%.1e\"" << std::endl
  << std::endl
  << "splot \"" + directory_of_data + "/conductivity.dat\" "
  << "using ($1/" << unit_scale << "):"
  << "($2/" << unit_scale << "):3" << std::endl
  << std::endl;
  
  ofs_script.close();
  return;
}
