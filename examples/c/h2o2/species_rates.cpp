#include "species_rates.hpp"
#include <math.h>
#include <stdlib.h>




void ndot_reset(double *__restrict__ dphi)
{
  for (int i = 0; i <= 9; ++i)
  {
    int const j = 0;

    dphi[i] = 0.0;
  }
}
void wdot_reset(double *__restrict__ wdot)
{
  for (int i = 0; i <= 8; ++i)
  {
    int const j = 0;

    wdot[i] = 0.0;
  }
}
void get_concentrations(double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ conc)
{
  double T_val;
  double V_inv;
  double n;
  double n_sum;

  {
    int const j = 0;

    T_val = fmin(10000.0, fmax(100.0, phi[0]));
    V_inv = 1.0 / fmax(1e-30, phi[1]);
    n_sum = 0.0;
    conc[8] = P_arr[0] / (8314.4621 * T_val);
    for (int i = 0; i <= 7; ++i)
    {
      n = fmax(1e-300, phi[2 + i]);
      n_sum = n_sum + n;
      conc[i] = n * V_inv;
    }
    conc[8] = conc[8] + -1.0 * fmax(1e-300, n_sum) * V_inv;
  }
}
void a_only_simple(double const *__restrict__ phi, double *__restrict__ kf)
{
  {
    int const i = 3;

    {
      int const j = 0;

      kf[3] = simple_A[3];
    }
  }
}
void beta_int_simple(double const *__restrict__ phi, double *__restrict__ kf)
{
  double Tinv;
  double Tval;
  int b_end;
  int i_0;
  int i_1;
  double kf_temp;

  {
    int const j = 0;

    Tval = fmin(10000.0, fmax(100.0, phi[0]));
    Tinv = 1.0 / fmin(10000.0, fmax(100.0, phi[0]));
    for (int i = 0; i <= 3; ++i)
    {
      i_1 = simple_rtype_1_map[i];
      i_0 = simple_rtype_1_inds[i];
      b_end = simple_beta[i_0];
      kf_temp = simple_A[i_0];
      kf_temp = kf_temp * fast_powi(Tval, b_end);
      kf[i_1] = kf_temp;
    }
  }
}
void rateconst_fullsimple(double const *__restrict__ phi, double *__restrict__ kf)
{
  double Tinv;
  int i_0;
  int i_1;
  double logT;

  {
    int const j = 0;

    logT = log(fmin(10000.0, fmax(100.0, phi[0])));
    Tinv = 1.0 / fmin(10000.0, fmax(100.0, phi[0]));
    for (int i = 0; i <= 22; ++i)
    {
      i_1 = simple_rtype_2_map[i];
      i_0 = simple_rtype_2_inds[i];
      kf[i_1] = exp(fmin(690.775527898, simple_A[i_0] + logT * simple_beta[i_0] + -1.0 * simple_Ta[i_0] * Tinv));
    }
  }
}
void eval_thd_body_concs(double const *__restrict__ P_arr, double const *__restrict__ phi, double const *__restrict__ conc, double *__restrict__ thd_conc)
{
  int not_spec;
  int offset;
  int spec_end;
  int spec_ind;
  double thd_temp;

  for (int i = 0; i <= 5; ++i)
  {
    int const j = 0;

    offset = thd_offset[i];
    spec_end = thd_offset[1 + i];
    not_spec = thd_type[i] != 2;
    thd_temp = P_arr[0] * not_spec / (8314.4621 * phi[0]);
    for (int ispec = offset; ispec <= -1 + spec_end; ++ispec)
    {
      spec_ind = thd_spec[ispec];
      thd_temp = thd_temp + (thd_eff[ispec] + -1.0 * not_spec) * conc[spec_ind];
    }
    thd_conc[i] = thd_temp;
  }
}
void rateconst_fullfall(double const *__restrict__ phi, double *__restrict__ kf_fall)
{
  double Tinv;
  double logT;

  {
    int const j = 0;

    logT = log(fmin(10000.0, fmax(100.0, phi[0])));
    Tinv = 1.0 / fmin(10000.0, fmax(100.0, phi[0]));
    {
      int const i = 0;

      kf_fall[0] = exp(fmin(690.775527898, fall_A[0] + logT * fall_beta[0] + -1.0 * fall_Ta[0] * Tinv));
    }
  }
}
void red_pres(double const *__restrict__ phi, double const *__restrict__ thd_conc, double const *__restrict__ kf, double const *__restrict__ kf_fall, double *__restrict__ Pr)
{
  double k0;
  double kinf;

  {
    int const i = 20;

    {
      int const j = 0;

      if (!fall_type[0])
        kinf = kf[20];
      if (fall_type[0])
        kinf = kf_fall[0];
      if (!fall_type[0])
        k0 = kf_fall[0];
      if (fall_type[0])
        k0 = kf[20];
      Pr[0] = thd_conc[5] * k0 / kinf;
    }
  }
}
void fall_troe(double const *__restrict__ Pr, double const *__restrict__ phi, double *__restrict__ Fi, double *__restrict__ Fcent, double *__restrict__ Atroe, double *__restrict__ Btroe)
{
  double Atroe_temp;
  double Btroe_temp;
  double Fcent_temp;
  double T;
  double logFcent;
  double logPr;

  {
    int const j = 0;

    T = fmin(10000.0, fmax(100.0, phi[0]));
    {
      int const i = 0;

      logPr = log10(fmax(1e-300, Pr[0]));
      Fcent_temp = troe_a[0] * exp(fmin(690.775527898, -1.0 * T * troe_T1[0])) + (1.0 + -1.0 * troe_a[0]) * exp(fmin(690.775527898, -1.0 * T * troe_T3[0]));
      if (troe_T2[0] != 0.0)
        Fcent_temp = Fcent_temp + exp(fmin(690.775527898, -1.0 * troe_T2[0] / T));
      Fcent[0] = Fcent_temp;
      logFcent = log10(fmax(1e-300, Fcent_temp));
      Btroe_temp = -1.1762 * logFcent + -1.0 * 0.14 * logPr + 0.806;
      Btroe[0] = Btroe_temp;
      Atroe_temp = -0.67 * logFcent + logPr + -0.4;
      Atroe[0] = Atroe_temp;
      Fi[0] = pow(fmax(1e-300, Fcent_temp), 1.0 / (Atroe_temp * Atroe_temp / (Btroe_temp * Btroe_temp) + 1.0));
    }
  }
}
void rateconst_Kc(double const *__restrict__ b, double *__restrict__ Kc, double const *__restrict__ kf, double *__restrict__ kr)
{
  double B_sum;
  double Kc_temp;
  double P_sum;
  int P_sum_end;
  double P_val;
  int net_nu;
  int offset;
  int spec_end;
  int spec_ind;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    offset = net_reac_to_spec_offsets[i];
    if (!(nu_sum[i] > 0))
      P_val = 8314.4621 / 101325.0;
    if (nu_sum[i] > 0)
      P_val = 101325.0 / 8314.4621;
    P_sum_end = abs(nu_sum[i]);
    P_sum = fast_powi(P_val, P_sum_end);
    B_sum = 0.0;
    spec_end = net_reac_to_spec_offsets[1 + i];
    for (int ispec = offset; ispec <= -1 + spec_end; ++ispec)
    {
      net_nu = reac_to_spec_nu[2 * ispec] + -1 * reac_to_spec_nu[1 + 2 * ispec];
      spec_ind = rxn_to_spec[ispec];
      if (net_nu != 0)
        B_sum = B_sum + net_nu * b[spec_ind];
    }
    B_sum = exp(fmin(690.775527898, B_sum));
    Kc_temp = P_sum * B_sum;
    kr[i] = kf[i] / Kc_temp;
    Kc[i] = Kc_temp;
  }
}
void ci_thd(double const *__restrict__ thd_conc, double *__restrict__ pres_mod)
{
  for (int i = 0; i <= 4; ++i)
  {
    int const j = 0;

    pres_mod[i] = thd_conc[i];
  }
}
void ci_fall(double const *__restrict__ Fi, double const *__restrict__ Pr, double *__restrict__ pres_mod)
{
  double ci_temp;

  {
    int const i = 0;

    {
      int const j = 0;

      ci_temp = Fi[0] / (1.0 + Pr[0]);
      if (!fall_type[0])
        ci_temp = ci_temp * Pr[0];
      pres_mod[5] = ci_temp;
    }
  }
}
void rop_eval_fwd(double const *__restrict__ conc, double const *__restrict__ kf, double *__restrict__ rop_fwd)
{
  double rop_temp;
  int spec_ind;
  int spec_offset;
  int spec_offset_next;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    rop_temp = kf[i];
    spec_offset_next = net_reac_to_spec_offsets[1 + i];
    spec_offset = net_reac_to_spec_offsets[i];
    for (int ispec = spec_offset; ispec <= -1 + spec_offset_next; ++ispec)
    {
      spec_ind = rxn_to_spec[ispec];
      rop_temp = rop_temp * fast_powi(conc[spec_ind], reac_to_spec_nu[1 + 2 * ispec]);
    }
    rop_fwd[i] = rop_temp;
  }
}
void rop_eval_rev(double const *__restrict__ conc, double const *__restrict__ kr, double *__restrict__ rop_rev)
{
  double rop_temp;
  int spec_ind;
  int spec_offset;
  int spec_offset_next;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    rop_temp = kr[i];
    spec_offset_next = net_reac_to_spec_offsets[1 + i];
    spec_offset = net_reac_to_spec_offsets[i];
    for (int ispec = spec_offset; ispec <= -1 + spec_offset_next; ++ispec)
    {
      spec_ind = rxn_to_spec[ispec];
      rop_temp = rop_temp * fast_powi(conc[spec_ind], reac_to_spec_nu[2 * ispec]);
    }
    rop_rev[i] = rop_temp;
  }
}
void rop_net_fixed(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ pres_mod, double *__restrict__ rop_net)
{
  int i_0;
  double net_rate;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    net_rate = rop_fwd[i];
    net_rate = net_rate + -1.0 * rop_rev[i];
    i_0 = thd_mask[i];
    if (i_0 >= 0)
      net_rate = net_rate * pres_mod[i_0];
    rop_net[i] = net_rate;
  }
}
void spec_rates(double const *__restrict__ rop_net, double *__restrict__ wdot)
{
  double net_rate;
  int nu;
  int offset;
  int offset_next;
  int spec_ind;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    net_rate = rop_net[i];
    offset_next = net_reac_to_spec_offsets[1 + i];
    offset = net_reac_to_spec_offsets[i];
    for (int ispec = offset; ispec <= -1 + offset_next; ++ispec)
    {
      nu = reac_to_spec_nu[2 * ispec] + -1 * reac_to_spec_nu[1 + 2 * ispec];
      spec_ind = rxn_to_spec[ispec];
      wdot[spec_ind] = wdot[spec_ind] + nu * net_rate;
    }
  }
}
void get_molar_rates(double const *__restrict__ phi, double *__restrict__ dphi, double const *__restrict__ wdot)
{
  double V_val;

  {
    int const j = 0;

    V_val = fmax(1e-30, phi[1]);
    for (int i = 0; i <= 7; ++i)
      dphi[2 + i] = V_val * wdot[i];
  }
}
void temperature_rate(double const *__restrict__ h, double const *__restrict__ cp, double const *__restrict__ conc, double *__restrict__ dphi, double const *__restrict__ wdot)
{
  double lower;
  double upper;

  {
    int const j = 0;

    upper = 0.0;
    lower = 0.0;
    for (int i = 0; i <= 8; ++i)
    {
      lower = lower + conc[i] * cp[i];
      upper = upper + h[i] * wdot[i];
    }
    dphi[0] = dphi[0] + -1.0 * upper / lower;
  }
}
void get_extra_var_rates(double const *__restrict__ wdot, double *__restrict__ dphi, double const *__restrict__ phi, double const *__restrict__ P_arr)
{
  double T;
  double V_val;
  double dE;

  {
    int const j = 0;

    V_val = fmax(1e-30, phi[1]);
    T = fmin(10000.0, fmax(100.0, phi[0]));
    dphi[1] = V_val * dphi[0] / T;
    dE = 0.0;
    for (int i = 0; i <= 7; ++i)
      dE = dE + (1.0 + -1.0 * mw_factor[i]) * wdot[i];
    dphi[1] = dphi[1] + V_val * dE * T * 8314.4621 / P_arr[0];
  }
}

void species_rates(double const *__restrict__ t, double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ dphi, double *__restrict__ rwk)
{
    double* __restrict__ Atroe = rwk + 0 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ Btroe = rwk + 1 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ Fcent = rwk + 2 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ Fi = rwk + 3 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ Kc = rwk + 4 * work_size + 28 * omp_get_thread_num();
    double* __restrict__ Pr = rwk + 32 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ b = rwk + 33 * work_size + 9 * omp_get_thread_num();
    double* __restrict__ conc = rwk + 42 * work_size + 9 * omp_get_thread_num();
    double* __restrict__ cp = rwk + 51 * work_size + 9 * omp_get_thread_num();
    double* __restrict__ h = rwk + 89 * work_size + 9 * omp_get_thread_num();
    double* __restrict__ kf = rwk + 98 * work_size + 28 * omp_get_thread_num();
    double* __restrict__ kf_fall = rwk + 126 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ kr = rwk + 127 * work_size + 28 * omp_get_thread_num();
    double* __restrict__ pres_mod = rwk + 155 * work_size + 6 * omp_get_thread_num();
    double* __restrict__ rop_fwd = rwk + 161 * work_size + 28 * omp_get_thread_num();
    double* __restrict__ rop_net = rwk + 189 * work_size + 28 * omp_get_thread_num();
    double* __restrict__ rop_rev = rwk + 217 * work_size + 28 * omp_get_thread_num();
    double* __restrict__ thd_conc = rwk + 245 * work_size + 6 * omp_get_thread_num();
    double* __restrict__ wdot = rwk + 251 * work_size + 9 * omp_get_thread_num();
    ndot_reset(dphi);
    
    wdot_reset(wdot);
    
    get_concentrations(phi, P_arr, conc);
    
    a_only_simple(phi, kf);
    
    beta_int_simple(phi, kf);
    
    rateconst_fullsimple(phi, kf);
    
    eval_thd_body_concs(P_arr, phi, conc, thd_conc);
    
    rateconst_fullfall(phi, kf_fall);
    
    red_pres(phi, thd_conc, kf, kf_fall, Pr);
    
    fall_troe(Pr, phi, Fi, Fcent, Atroe, Btroe);
    
    eval_b(phi, b);
    
    rateconst_Kc(b, Kc, kf, kr);
    
    ci_thd(thd_conc, pres_mod);
    
    ci_fall(Fi, Pr, pres_mod);
    
    rop_eval_fwd(conc, kf, rop_fwd);
    
    rop_eval_rev(conc, kr, rop_rev);
    
    rop_net_fixed(rop_fwd, rop_rev, pres_mod, rop_net);
    
    spec_rates(rop_net, wdot);
    
    get_molar_rates(phi, dphi, wdot);
    
    eval_h(phi, h);
    
    eval_cp(phi, cp);
    
    temperature_rate(h, cp, conc, dphi, wdot);
    
    get_extra_var_rates(wdot, dphi, phi, P_arr);
}
