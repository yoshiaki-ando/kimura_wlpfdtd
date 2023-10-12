#include <fftw3.h>
#include "wlp_fdtd2d.h"

void fourier_transform(double *field, std::complex <double> *FTdata){
  std::complex <double> zj(0.0, 1.0);
  double fs = 1.0/DT; //サンプリングレート
  double Df = (fs / NT); //基本周波数
  for(int n = 0; n < NT+1; n++){
    double t = n*DT;
    for(int m = 0; m < NT/2; m++){
      double omega = 2.0 * M_PI * Df * m;
      FTdata[m] += field[n] * exp( - zj * omega * t) * DT;
    }
  }

  return;
}

void fast_fourier_transform(double *field, std::complex <double> *FTdata, int nt, double dt){
  
  const int REAL = 0;
  const int IMAG = 1;
  std::complex <double> zj(0.0, 1.0);

  fftw_complex *in, *out;//入力と出力を定義
  fftw_plan p;//高速フーリエ変換の実行形式変数の定義
  const int fftsize = nt;

  in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *fftsize);
  p = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


  for(int n = 0; n < fftsize; n++){//入力を入力
    in[n][REAL] = field[n];
    in[n][IMAG] = 0.0;
  }

  fftw_execute(p);//高速フーリエ変換を実行

  for(int n=0; n<fftsize/2; n++){
    FTdata[n] = (out[n][REAL] + zj * out[n][IMAG])*dt;
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  return;
}

void output_spectrum(std::complex <double> *FTdata, int forN, double delta, std::string file_name){
  std::ofstream ofs_FTobs( file_name );
  for( int n = 0; n <= forN; n++ ){
    ofs_FTobs << n * delta << " " << std::scientific << std::abs(FTdata[n]) << std::endl;
  }
  ofs_FTobs.close();

  return;
}

void calc_transfer_function(std::complex<double> *A, std::complex<double> *B, double *T){
  for(int n = 0; n <= NT/2; n++){
    T[n] = std::abs(A[n])/std::abs(B[n]);
  }

  return;
}
