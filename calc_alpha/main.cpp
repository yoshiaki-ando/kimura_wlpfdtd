#include "params.h"

int main(void){
  const int n_ilp = 200;

  // 伝達関数を入力
  double tmp;
  std::ifstream ifs(file_dir);
  std::string tmps; 
  std::vector<double> tran(2);
  std::vector<double> f(2);
  for(int n = 0; n <= 1; n++) {
    ifs >> f[n] >> tran[n];
  }
  const double df = f[1] - f[0];
  const int nf_first = std::round(first_freq/df);
  const int nf_last = std::round(last_freq/df);
  f.resize(nf_last+1);
  tran.resize(nf_last+1);

  for(int n = 2; n <= nf_last; n++) {
    ifs >> f[n] >> tran[n];
  }
  ifs.close();

  bool fl_same_freq = false;
  if(nf_first==nf_last) fl_same_freq=true;

  std::vector<std::vector<double>> current(n_ilp, std::vector<double>(nf_last+1, 0.0));
  std::vector<std::vector<double>> b_lp(n_ilp, std::vector<double>(nf_last+1, 0.0));

  std::string input_dir;
  if(eq_type == "pre") input_dir = file_ilp_pre;
  if(eq_type == "near") input_dir = file_ilp_near;
  std::ifstream ifs_input(input_dir);
  for(int i = 0; i < n_ilp; i++) {
    double a, b;
    ifs_input >> a >> b;
    for(int n = 0; n <= nf_last; n++) {
      double f = n*df;
      current[i][n] = std::pow(10, a*std::log10(f/(9.e-3))+b);
      // b_lp = I_lp(f)H(f)
      b_lp[i][n] = current[i][n] * tran[n];
    }
  }
  ifs_input.close();

  std::vector<double> int_field(n_ilp, 0.0);
  for(int i = 0; i < n_ilp; i++) {
    if(!fl_same_freq){
      for(int n = nf_first; n < nf_last; n++) {
        // ∫(b_lp^2)df
        int_field[i] += (std::pow(b_lp[i][n],2)+std::pow(b_lp[i][n+1],2))*df/2.0;
      }
      // (1/BW)∫(b_lp^2)df 
      if(fl_per_band) int_field[i] /= (last_freq - first_freq);
    }else{
      int_field[i] = b_lp[i][nf_first];
    }
  }

  std::vector<double> alpha(n_ilp, 0.0);
  for(int i = 0; i < n_ilp; i++) {
      if(!fl_same_freq) alpha[i] = obs_value/std::sqrt(int_field[i]);
      else alpha[i] = obs_value/int_field[i];
  }

  std::sort(alpha.begin(), alpha.end());
  // for(int i = 0; i < n_ilp; i++) {
  //   std::cout << alpha[i] << std::endl;
  // }

  std::ofstream ofs("output_data/" + eq_name + ".dat");
    ofs << "#obs_value = " << obs_value << std::endl;
    ofs << "#magnitude alpha\n#(log(alpha[0])+log(alpha[199])/2" << std::endl;
    ofs << eq_mag << " " << alpha[0] << std::endl;
    ofs << eq_mag << " " << alpha[199] << std::endl;
    ofs << (std::log10(alpha[0]) + std::log10(alpha[199]))/2.0 << std::endl;
  ofs.close();

  return 0;
}
