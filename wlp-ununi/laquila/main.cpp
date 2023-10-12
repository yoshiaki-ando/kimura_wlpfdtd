#include "wlp_fdtd2d.h"
#include <experimental/filesystem>

int main(void){
  double start_pre = get_time_sec();

  constexpr double p_max = P_ORDER;
  constexpr int n_max = NT+1;

  const std::string str_date = get_date_str();
  std::string directory_of_data =  LOC_DATA+str_date;
  //create a directory
	std::experimental::filesystem::create_directory(directory_of_data);
  
  constexpr int i_obs = int(X_OBSERVE/DX + 0.5 + NX3);
  constexpr int j_obs = int(Y_OBSERVE/DY + 0.5 + NY3);
  constexpr int i_src = int(X_SOURCE/DX + 0.5 + NX3);
  constexpr int j_src = int(Y_SOURCE/DY + 0.5 + NY3);
  std::cout << "obs_point = (" << i_obs << ", " << j_obs << ")" << std::endl;

  // 時間波形
  double *hx = allocate_memory_d1d(n_max, 0.0);
  double **hx_up = allocate_memory_d2d(2, n_max, 0.0);
  double *hy = allocate_memory_d1d(n_max, 0.0);
  double *hy_left = allocate_memory_d1d(n_max, 0.0);
  double *ez = allocate_memory_d1d(n_max, 0.0);
  double *jz = allocate_memory_d1d(n_max, 0.0);

  // 展開係数和
  double **sum_coef_ez = allocate_memory_d2d(NX7+1, NY7+1, 0.0);
  double **sum_coef_ezx = allocate_memory_d2d(NX7+1, NY7+1, 0.0);
  double **sum_coef_ezy = allocate_memory_d2d(NX7+1, NY7+1, 0.0);
  double **sum_coef_hx = allocate_memory_d2d(NX7+1, NY7, 0.0);
  double **sum_coef_hy = allocate_memory_d2d(NX7, NY7+1, 0.0);  

  //dx,dyを計算
  double *dx = allocate_memory_d1d(NX7, 0.0);
  double *dy = allocate_memory_d1d(NY7, 0.0);
  calc_d(dx, NX2, NX3, NX4, NX5, NX7, dx_max, DX);
  calc_d(dy, NY2, NY3, NY4, NY5, NY7, dy_max, DY);
  double *dx_h = allocate_memory_d1d(NX7, 0.0);
  double *dy_h = allocate_memory_d1d(NY7, 0.0);
  calc_d_h(dx_h, dx, NX2, NX3, NX4, NX5, NX7);
  calc_d_h(dy_h, dy, NY2, NY3, NY4, NY5, NY7);

  //セルの位置を計算
  double *x_pos = allocate_memory_d1d(NX7+1, 0.0);
  double *y_pos = allocate_memory_d1d(NY7+1, 0.0);
  set_coordinate(dx, NX7+1, x_pos, -X_BOUNDARY, dy, NY7+1, y_pos, -Y_SURFACE); //dx, dy 追加

  std::cout << "(x_obs, y_obs) = "<< x_pos[i_obs]*1e-3 << "km, " << y_pos[j_obs]*1e-3 << "km" << std::endl;
  std::cout << "(x_src, y_src) = "<< x_pos[i_src]*1e-3 << "km, " << y_pos[j_src]*1e-3 << "km" << std::endl;

  //展開係数
  double **coef_ez = allocate_memory_d2d(NX7+1, NY7+1, 0.0);
  double **coef_hy = allocate_memory_d2d(NX7, NY7+1, 0.0);
  double **coef_hx = allocate_memory_d2d(NX7+1, NY7, 0.0);
  double *coef_jz   = allocate_memory_d1d(p_max+1, 0.0);
  calc_coef_jz(coef_jz, p_max);

  //ノード番号(右辺の一次元行列の番号を割り当てる)
  double ***node = allocate_memory_d3d(NX7+1, NY7+1, 3, -1.0);
  int total_unknowns = 0;
  assign_node_number(node, &total_unknowns);

  // モデルの導電率とその計算(pmlの導電率設定を含む)
  double **conduct_media = allocate_memory_d2d(NX7+1, NY7+1, 0.0);
  double *conduct_pmlx_ez = allocate_memory_d1d(NX7+1, 0.0);
  double *conduct_pmly_ez = allocate_memory_d1d(NY7+1, 0.0);
  double *conduct_pmlx_h = allocate_memory_d1d(NX7+1, 0.0);
  double *conduct_pmly_h = allocate_memory_d1d(NY7+1, 0.0);
  calc_conductivity(x_pos, y_pos, conduct_media, conduct_pmlx_ez, conduct_pmly_ez, conduct_pmlx_h, conduct_pmly_h);

  double *Dx = allocate_memory_d1d(NX7, 0.0);
  double *Dy = allocate_memory_d1d(NY7, 0.0);
  double *Cx = allocate_memory_d1d(NX7, 0.0);
  double *Cy = allocate_memory_d1d(NY7, 0.0);
  calc_const_param(Dx, Dy, Cx, Cy,
                   dx, dy, dx_h, dy_h);
  
  // 各定数とその計算(導電率設定)
  double *const_hx1 = allocate_memory_d1d(NY7, 0.0);
  double *const_hx2 = allocate_memory_d1d(NY7, 0.0);
  double *const_hy1 = allocate_memory_d1d(NX7, 0.0);
  double *const_hy2 = allocate_memory_d1d(NX7, 0.0);
  calc_const_pml(dx, dy,
                const_hx1, const_hx2,
		            const_hy1, const_hy2,
                conduct_pmlx_ez, conduct_pmly_ez,
                conduct_pmlx_h, conduct_pmly_h
                );

  // 左辺の係数行列(更新なし)を生成
  css *S;
  csn *N;
  int nmat;
  compose_coef_matrix(dx, dy, dx_h, dy_h,
                      total_unknowns, &nmat, &S, &N, 
                      Dx, Dy, Cx, Cy,
                      const_hx1, const_hx2,
                      const_hy1, const_hy2,
                      conduct_media,
                      conduct_pmlx_ez, conduct_pmly_ez,
                      node);

  // 一次元の展開係数行列
  double *b = allocate_memory_d1d(total_unknowns, 0.0);
  double *x = allocate_memory_d1d(total_unknowns, 0.0);

  // ラゲール多項式の1次前、2次前の時間波形
  // (漸化式で定義されるためp_max分のラゲール多項式の配列を用意する必要はない)
  double **w = allocate_memory_d2d(3, n_max, 0.0);
  double **v = allocate_memory_d2d(3, n_max, 0.0);
  double *wlp = allocate_memory_d1d(n_max, 0.0);
  
  double time_pre = get_time_sec() - start_pre;
  std::cout << "前処理... 完了. / " << " " << time_pre << " [s]" << std::endl;

  double start_fdtd = get_time_sec();
  
  for(int p=0; p<= p_max; p++){

    //右辺の一次元行列を生成(更新)
    compose_right_side(dx, dy, dx_h, dy_h,
                       b,
                       coef_jz,
                       sum_coef_ez,
                       sum_coef_ezx, sum_coef_ezy,
                       sum_coef_hx, sum_coef_hy,
                       Dx, Dy,
                       const_hx1, const_hx2,
                       const_hy1, const_hy2,
                       node, p);

    //LU solve
    cs_ipvec(nmat, N->Pinv, b, x);
    cs_lsolve(N->L, x);
    cs_usolve(N->U, x);
    cs_ipvec(nmat, S->Q, x, b);

    //電界、磁界の展開係数和を更新
    update_e(coef_ez,
             sum_coef_ez,
             sum_coef_ezx, sum_coef_ezy,
             b, // bには解Ez_pが代入されている
             node,
             p);
    update_h(dx, dy,
             coef_hx, coef_hy,
             sum_coef_hx, sum_coef_hy,
             coef_ez,
             Cx, Cy,
             const_hx1, const_hx2,
             const_hy1, const_hy2,
             node,
             p);

    calc_wlp_for_output(p, w, v, wlp, DT, NT);
    // 展開係数から時間波形に展開
    expand_field_p(coef_jz[p], jz, wlp, p);
    expand_field_p(coef_hx[i_obs][j_obs], hx, wlp, p);
    for(int i = 0; i < 2; i++)
      expand_field_p(coef_hx[i_obs][j_obs+(i+1)], hx_up[i], wlp, p);
    expand_field_p(coef_hy[i_obs][j_obs], hy, wlp, p);
    expand_field_p(coef_hy[i_obs-1][j_obs], hy_left, wlp, p);
    expand_field_p(coef_ez[i_obs][j_obs], ez, wlp, p);
    display_progress("展開係数導出", p, p_max);

  }

  double time_fdtd = get_time_sec() - start_fdtd;
  std::cout << "展開係数導出... 完了. / " << " " << time_fdtd << " [s]." << std::endl;

  // メモリ解放
  cs_sfree (S) ;
  cs_nfree (N) ;
  free_memory_d2d(sum_coef_ez);
  free_memory_d2d(sum_coef_ezx);
  free_memory_d2d(sum_coef_ezy);
  free_memory_d2d(sum_coef_hx);
  free_memory_d2d(sum_coef_hy);
  delete [] wlp;
  free_memory_d2d(w);
  free_memory_d2d(v);

  // 以降, 時間波形出力, フーリエ変換, スペクトル出力, 伝達関数計算などFDTD以外の処理
  double start_post = get_time_sec();

  output_field_1d(jz, n_max, DT, directory_of_data + "/jz.dat");
  output_field_1d(ez, n_max, DT, directory_of_data + "/ez.dat");
  std::complex <double> *FT_jz = new std::complex <double> [NT/2+1];
  fast_fourier_transform(jz, FT_jz);
  output_spectrum(FT_jz, NT/2, DF, directory_of_data + "/ft_jz.dat");
  
  std::complex <double> *FT_hx = new std::complex <double> [NT/2+1];
  std::complex <double> *FT_hx_up1 = new std::complex <double> [NT/2+1];
  std::complex <double> *FT_hx_up2 = new std::complex <double> [NT/2+1];
  std::complex <double> *FT_hy = new std::complex <double> [NT/2+1];
  std::complex <double> *FT_hy_left = new std::complex <double> [NT/2+1];
  std::complex <double> *FT_ez = new std::complex <double> [NT/2+1];
  fast_fourier_transform(hx, FT_hx);
  fast_fourier_transform(hx_up[0], FT_hx_up1);
  fast_fourier_transform(hx_up[1], FT_hx_up2);
  fast_fourier_transform(hy, FT_hy);
  fast_fourier_transform(hy_left, FT_hy_left);
  fast_fourier_transform(ez, FT_ez);

  double *hx_div_jz = allocate_memory_d1d(NT/2+1, 0.0);
  double *hx_up1_div_jz = allocate_memory_d1d(NT/2+1, 0.0);
  double *hx_up2_div_jz = allocate_memory_d1d(NT/2+1, 0.0);
  double *hy_div_jz = allocate_memory_d1d(NT/2+1, 0.0);
  double *hy_left_div_jz = allocate_memory_d1d(NT/2+1, 0.0);
  double *ez_div_jz = allocate_memory_d1d(NT/2+1, 0.0);
  calc_transfer_function(FT_hx, FT_jz, hx_div_jz);
  calc_transfer_function(FT_hx_up1, FT_jz, hx_up1_div_jz);
  calc_transfer_function(FT_hx_up2, FT_jz, hx_up2_div_jz);
  calc_transfer_function(FT_hy, FT_jz, hy_div_jz);
  calc_transfer_function(FT_hy_left, FT_jz, hy_left_div_jz);
  calc_transfer_function(FT_ez, FT_jz, ez_div_jz);

  std::ofstream ofs_bx(directory_of_data + "/ft_bx_div_Iz.dat");
  for(int i = 0; i < NT/2+1; i++){
    ofs_bx << i*DF << " " 
    << ( ( 3.0/8.0*hx_up2_div_jz[i] ) - (5.0/4.0*hx_up1_div_jz[i]) + (15.0/8.0*hx_div_jz[i]) )*MU0/DX/DY 
    << std::endl;
  }
  ofs_bx.close();
  std::ofstream ofs_by(directory_of_data + "/ft_by_div_Iz.dat");
  for(int i = 0; i < NT/2+1; i++){
    ofs_by << i*DF << " " 
    << ( ( hy_div_jz[i] + hy_left_div_jz[i] )/2 )*MU0/DX/DY << std::endl;
  }
  ofs_by.close();
  std::ofstream ofs_ez(directory_of_data + "/ft_ez_div_Iz.dat");
  for(int i = 0; i < NT/2+1; i++){
    ofs_ez << i*DF << " " << ez_div_jz[i] /DX/DY << std::endl;
  }
  ofs_ez.close();
  
  double time_post = get_time_sec() - start_post; 
  std::cout << "後処理... 完了. / " << " " << time_post << " [s]" << std::endl;
  double time_total = time_pre + time_fdtd + time_post;
  std::cout << "全計算 完了. / " << " " << time_total << " [s]." << std::endl;

  output_conductivity(x_pos, y_pos, conduct_media, directory_of_data);
  make_profile(time_pre, time_fdtd, time_post, time_total, directory_of_data + "/profile.txt"); 

  return 0;
}
