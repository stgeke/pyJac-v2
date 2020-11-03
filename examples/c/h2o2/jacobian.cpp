#include "jacobian.hpp"
#include <math.h>
#include <stdlib.h>




void reset_arrays(double *__restrict__ jac)
{
  int col;
  int row;

  for (int i = 0; i <= 99; ++i)
  {
    int const j = 0;

    col = (i % 10);
    row = (i / 10);
    jac[10 * row + col] = 0.0;
  }
}
void dRopidnj(double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ conc, double *__restrict__ jac)
{
  double Sj_fwd;
  double Sj_rev;
  double ci;
  int i_0;
  int inner_offset;
  int inner_offset_next;
  double kf_i;
  double kr_i;
  int net_offset;
  int net_offset_next;
  int nu_fwd;
  int nu_k;
  int nu_rev;
  int spec_inner;
  int spec_j;
  int spec_k;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    kr_i = 0.0;
    if (rev_mask[i] >= 0)
      kr_i = kr[i];
    kf_i = kf[i];
    inner_offset_next = net_reac_to_spec_offsets[1 + i];
    inner_offset = net_reac_to_spec_offsets[i];
    net_offset_next = net_reac_to_spec_offsets[1 + i];
    net_offset = net_reac_to_spec_offsets[i];
    i_0 = thd_mask[i];
    ci = 1.0;
    if (thd_mask[i] >= 0)
      ci = pres_mod[i_0];
    for (int net_ind_k = net_offset; net_ind_k <= -1 + net_offset_next; ++net_ind_k)
      if (-1 + -1 * inner_offset + inner_offset_next >= 0)
      {
        spec_k = rxn_to_spec[net_ind_k];
        if (spec_k != 8)
        {
          nu_k = reac_to_spec_nu[2 * net_ind_k] + -1 * reac_to_spec_nu[1 + 2 * net_ind_k];
          for (int net_ind_j = inner_offset; net_ind_j <= -1 + inner_offset_next; ++net_ind_j)
          {
            spec_j = rxn_to_spec[net_ind_j];
            if (spec_j != 8)
            {
              Sj_rev = reac_to_spec_nu[2 * net_ind_j];
              Sj_fwd = reac_to_spec_nu[1 + 2 * net_ind_j];
              for (int net_ind_inner = inner_offset; net_ind_inner <= -1 + inner_offset_next; ++net_ind_inner)
              {
                nu_rev = reac_to_spec_nu[2 * net_ind_inner];
                nu_fwd = reac_to_spec_nu[1 + 2 * net_ind_inner];
                spec_inner = rxn_to_spec[net_ind_inner];
                if (spec_inner == spec_j)
                {
                  nu_rev = nu_rev + -1;
                  nu_fwd = nu_fwd + -1;
                }
                Sj_rev = Sj_rev * fast_powi(fmax(1e-300, conc[spec_inner]), nu_rev);
                Sj_fwd = Sj_fwd * fast_powi(fmax(1e-300, conc[spec_inner]), nu_fwd);
              }
              jac[10 * (spec_k + 2) + spec_j + 2] = jac[10 * (spec_k + 2) + spec_j + 2] + (kf_i * Sj_fwd + -1.0 * kr_i * Sj_rev) * ci * nu_k;
            }
          }
        }
      }
  }
}
void dRopidnj_ns(double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ conc, double *__restrict__ jac)
{
  double Sns_fwd;
  double Sns_rev;
  double ci;
  int inner_offset;
  int inner_offset_next;
  double jac_updater;
  double kf_i;
  double kr_i;
  int net_offset;
  int net_offset_next;
  int nu_fwd;
  int nu_k;
  int nu_rev;
  int spec_inner;
  int spec_k;

  {
    int const i = 8;

    {
      int const j = 0;

      kr_i = 0.0;
      if (rev_mask[8] >= 0)
        kr_i = kr[rev_mask[8]];
      kf_i = kf[8];
      inner_offset_next = net_reac_to_spec_offsets[9];
      inner_offset = net_reac_to_spec_offsets[8];
      net_offset_next = net_reac_to_spec_offsets[9];
      net_offset = net_reac_to_spec_offsets[8];
      ci = 1.0;
      if (thd_mask[8] >= 0)
        ci = pres_mod[thd_mask[8]];
      for (int net_ind_k = net_offset; net_ind_k <= -1 + net_offset_next; ++net_ind_k)
        if (-1 + -1 * inner_offset + inner_offset_next >= 0)
        {
          spec_k = rxn_to_spec[net_ind_k];
          if (spec_k != 8)
          {
            nu_k = reac_to_spec_nu[2 * net_ind_k] + -1 * reac_to_spec_nu[1 + 2 * net_ind_k];
            Sns_rev = 1.0;
            Sns_fwd = 1.0;
            for (int net_ind_inner = inner_offset; net_ind_inner <= -1 + inner_offset_next; ++net_ind_inner)
            {
              nu_rev = reac_to_spec_nu[2 * net_ind_inner];
              nu_fwd = reac_to_spec_nu[1 + 2 * net_ind_inner];
              spec_inner = rxn_to_spec[net_ind_inner];
              if (spec_inner == 8)
              {
                Sns_rev = Sns_rev * nu_rev;
                nu_rev = nu_rev + -1;
              }
              Sns_rev = Sns_rev * fast_powi(fmax(1e-300, conc[spec_inner]), nu_rev);
              if (spec_inner == 8)
              {
                Sns_fwd = Sns_fwd * nu_fwd;
                nu_fwd = nu_fwd + -1;
              }
              Sns_fwd = Sns_fwd * fast_powi(fmax(1e-300, conc[spec_inner]), nu_fwd);
            }
            jac_updater = (kr_i * Sns_rev + -1.0 * kf_i * Sns_fwd) * ci * nu_k;
            for (int spec_j = 0; spec_j <= 7; ++spec_j)
              jac[10 * (spec_k + 2) + 2 + spec_j] = jac[10 * (spec_k + 2) + 2 + spec_j] + jac_updater;
          }
        }
    }
  }
}
void dci_thd_dnj(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double *__restrict__ jac)
{
  double dci;
  int i_0;
  int i_1;
  int nu_k;
  double ropi;
  int rxn_off;
  int rxn_off_next;
  int spec_j;
  int spec_k;
  int thd_off;
  int thd_off_next;

  for (int i = 0; i <= 4; ++i)
  {
    int const j = 0;

    thd_off_next = thd_offset[1 + i];
    thd_off = thd_offset[i];
    i_0 = thd_map[i];
    ropi = rop_fwd[i_0];
    rxn_off_next = net_reac_to_spec_offsets[i_0 + 1];
    rxn_off = net_reac_to_spec_offsets[i_0];
    i_1 = rev_mask[i_0];
    if (rev_mask[i_0] >= 0)
      ropi = ropi + -1.0 * rop_rev[i_1];
    for (int spec_j_ind = thd_off; spec_j_ind <= -1 + thd_off_next; ++spec_j_ind)
      if (-1 + rxn_off_next + -1 * rxn_off >= 0)
      {
        spec_j = thd_spec[spec_j_ind];
        if (spec_j != 8)
        {
          dci = 0.0;
          if (thd_type[i] == 2)
            dci = 1.0;
          if (thd_type[i] == 1)
            dci = thd_eff[spec_j_ind] + -1.0;
          for (int spec_k_ind = rxn_off; spec_k_ind <= -1 + rxn_off_next; ++spec_k_ind)
            if (rxn_to_spec[spec_k_ind] != 8)
            {
              spec_k = rxn_to_spec[spec_k_ind];
              nu_k = reac_to_spec_nu[2 * spec_k_ind] + -1 * reac_to_spec_nu[1 + 2 * spec_k_ind];
              jac[10 * (spec_k + 2) + spec_j + 2] = jac[10 * (spec_k + 2) + spec_j + 2] + nu_k * dci * ropi;
            }
        }
      }
  }
}
void dci_thd_dnj_ns(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double *__restrict__ jac)
{
  double dci;
  int i_0;
  int i_1;
  double ns_thd_eff;
  int nu_k;
  double ropi;
  int rxn_off;
  int rxn_off_next;
  int spec_k;

  for (int i = 0; i <= 4; ++i)
  {
    int const j = 0;

    ns_thd_eff = thd_eff[thd_offset[1 + i] + -1];
    i_0 = thd_map[i];
    ropi = rop_fwd[i_0];
    rxn_off_next = net_reac_to_spec_offsets[i_0 + 1];
    rxn_off = net_reac_to_spec_offsets[i_0];
    i_1 = rev_mask[i_0];
    if (rev_mask[i_0] >= 0)
      ropi = ropi + -1.0 * rop_rev[i_1];
    dci = 0.0;
    if (thd_type[i] == 2)
      dci = -1.0;
    if (thd_type[i] == 1)
      dci = 1.0 + -1.0 * ns_thd_eff;
    for (int spec_k_ind = rxn_off; spec_k_ind <= -1 + rxn_off_next; ++spec_k_ind)
    {
      spec_k = rxn_to_spec[spec_k_ind];
      nu_k = reac_to_spec_nu[2 * spec_k_ind] + -1 * reac_to_spec_nu[1 + 2 * spec_k_ind];
      if (spec_k != 8)
        for (int spec_j = 0; spec_j <= 7; ++spec_j)
          jac[10 * (spec_k + 2) + 2 + spec_j] = jac[10 * (spec_k + 2) + 2 + spec_j] + nu_k * dci * ropi;
    }
  }
}
void dci_troe_dnj(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ Pr, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Fi, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac)
{
  double Fi_fac;
  double dFi;
  double dci;
  int i_2;
  double k0;
  double kinf;
  int nu_k;
  double ropi;
  int rxn_off;
  int rxn_off_next;
  int spec_j;
  int spec_k;
  int thd_off;
  int thd_off_next;

  {
    int const i = 0;

    {
      int const j = 0;

      ropi = rop_fwd[20];
      if (!fall_type[0])
        kinf = kf[20];
      if (fall_type[0])
        kinf = kf_fall[0];
      if (!fall_type[0])
        k0 = kf_fall[0];
      if (fall_type[0])
        k0 = kf[20];
      rxn_off_next = net_reac_to_spec_offsets[21];
      rxn_off = net_reac_to_spec_offsets[20];
      thd_off_next = thd_offset[6];
      thd_off = thd_offset[5];
      i_2 = rev_mask[20];
      if (rev_mask[20] >= 0)
        ropi = ropi + -1.0 * rop_rev[i_2];
      dFi = Atroe[0] * Atroe[0] + Btroe[0] * Btroe[0];
      dFi = -2.0 * Atroe[0] * Btroe[0] * (0.14 * Atroe[0] + Btroe[0]) * log(fmax(1e-300, Fcent[0])) / (fmax(1e-300, Pr[0]) * dFi * dFi * 2.302585092994046);
      Fi_fac = dFi;
      if (!fall_type[0])
        Fi_fac = Pr[0] * Fi_fac + 1.0;
      dFi = k0 * (Fi[0] * Fi_fac + -1.0 * pres_mod[5]) / (kinf * (Pr[0] + 1.0));
      for (int spec_j_ind = thd_off; spec_j_ind <= -1 + thd_off_next; ++spec_j_ind)
        if (-1 + rxn_off_next + -1 * rxn_off >= 0)
        {
          spec_j = thd_spec[spec_j_ind];
          if (spec_j != 8)
          {
            dci = 0.0;
            if (thd_type[5] == 2)
              dci = 1.0;
            if (thd_type[5] == 1)
              dci = thd_eff[spec_j_ind] + -1.0;
            for (int spec_k_ind = rxn_off; spec_k_ind <= -1 + rxn_off_next; ++spec_k_ind)
              if (rxn_to_spec[spec_k_ind] != 8)
              {
                spec_k = rxn_to_spec[spec_k_ind];
                nu_k = reac_to_spec_nu[2 * spec_k_ind] + -1 * reac_to_spec_nu[1 + 2 * spec_k_ind];
                jac[10 * (spec_k + 2) + spec_j + 2] = jac[10 * (spec_k + 2) + spec_j + 2] + nu_k * dci * ropi * dFi;
              }
          }
        }
    }
  }
}
void dci_troe_dnj_ns(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ Pr, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Fi, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac)
{
  double Fi_fac;
  double dFi;
  double dci;
  int i_2;
  double k0;
  double kinf;
  double ns_thd_eff;
  int nu_k;
  double ropi;
  int rxn_off;
  int rxn_off_next;
  int spec_k;

  {
    int const i = 0;

    {
      int const j = 0;

      ropi = rop_fwd[20];
      if (!fall_type[0])
        kinf = kf[20];
      if (fall_type[0])
        kinf = kf_fall[0];
      if (!fall_type[0])
        k0 = kf_fall[0];
      if (fall_type[0])
        k0 = kf[20];
      ns_thd_eff = thd_eff[thd_offset[6] + -1];
      rxn_off_next = net_reac_to_spec_offsets[21];
      rxn_off = net_reac_to_spec_offsets[20];
      i_2 = rev_mask[20];
      if (rev_mask[20] >= 0)
        ropi = ropi + -1.0 * rop_rev[i_2];
      dFi = Atroe[0] * Atroe[0] + Btroe[0] * Btroe[0];
      dFi = -2.0 * Atroe[0] * Btroe[0] * (0.14 * Atroe[0] + Btroe[0]) * log(fmax(1e-300, Fcent[0])) / (fmax(1e-300, Pr[0]) * dFi * dFi * 2.302585092994046);
      Fi_fac = dFi;
      if (!fall_type[0])
        Fi_fac = Pr[0] * Fi_fac + 1.0;
      dFi = k0 * (Fi[0] * Fi_fac + -1.0 * pres_mod[5]) / (kinf * (Pr[0] + 1.0));
      dci = 0.0;
      if (thd_type[5] == 2)
        dci = -1.0;
      if (thd_type[5] == 1)
        dci = 1.0 + -1.0 * ns_thd_eff;
      for (int spec_k_ind = rxn_off; spec_k_ind <= -1 + rxn_off_next; ++spec_k_ind)
      {
        spec_k = rxn_to_spec[spec_k_ind];
        nu_k = reac_to_spec_nu[2 * spec_k_ind] + -1 * reac_to_spec_nu[1 + 2 * spec_k_ind];
        if (spec_k != 8)
          for (int spec_j = 0; spec_j <= 7; ++spec_j)
            jac[10 * (spec_k + 2) + 2 + spec_j] = jac[10 * (spec_k + 2) + 2 + spec_j] + nu_k * dci * ropi * dFi;
      }
    }
  }
}
void cp_total(double const *__restrict__ cp, double const *__restrict__ conc, double *__restrict__ cp_tot)
{
  double spec_tot;

  {
    int const j = 0;

    spec_tot = 0.0;
    cp_tot[0] = 0.0;
    for (int i = 0; i <= 8; ++i)
      spec_tot = spec_tot + cp[i] * conc[i];
    cp_tot[0] = cp_tot[0] + spec_tot;
  }
}
void eval_dcp(double const *__restrict__ phi, double *__restrict__ dcp)
{
  double T;

  {
    int const j = 0;

    T = fmin(10000.0, fmax(100.0, phi[0]));
    for (int k = 0; k <= 8; ++k)
    {
      if (T < T_mid[k])
        dcp[k] = 8314.4621 * (T * (T * (4.0 * T * a_lo[4 + 7 * k] + 3.0 * a_lo[3 + 7 * k]) + 2.0 * a_lo[2 + 7 * k]) + a_lo[1 + 7 * k]);
      if (!(T < T_mid[k]))
        dcp[k] = 8314.4621 * (T * (T * (4.0 * T * a_hi[4 + 7 * k] + 3.0 * a_hi[3 + 7 * k]) + 2.0 * a_hi[2 + 7 * k]) + a_hi[1 + 7 * k]);
    }
  }
}
void eval_db(double const *__restrict__ phi, double *__restrict__ db)
{
  double T;
  double Tinv;

  {
    int const j = 0;

    Tinv = 1.0 / fmin(10000.0, fmax(100.0, phi[0]));
    T = fmin(10000.0, fmax(100.0, phi[0]));
    for (int k = 0; k <= 8; ++k)
    {
      if (T < T_mid[k])
        db[k] = T * (T * (T * a_lo[4 + 7 * k] / 5.0 + a_lo[3 + 7 * k] / 4.0) + a_lo[2 + 7 * k] / 3.0) + a_lo[1 + 7 * k] / 2.0 + Tinv * (a_lo[7 * k] + -1.0 + a_lo[5 + 7 * k] * Tinv);
      if (!(T < T_mid[k]))
        db[k] = T * (T * (T * a_hi[4 + 7 * k] / 5.0 + a_hi[3 + 7 * k] / 4.0) + a_hi[2 + 7 * k] / 3.0) + a_hi[1 + 7 * k] / 2.0 + Tinv * (a_hi[7 * k] + -1.0 + a_hi[5 + 7 * k] * Tinv);
    }
  }
}
void dTdot_dnj(double const *__restrict__ cp, double const *__restrict__ h, double const *__restrict__ cp_tot, double const *__restrict__ phi, double const *__restrict__ dphi, double *__restrict__ jac)
{
  double sum;

  for (int i = 0; i <= 7; ++i)
  {
    int const j = 0;

    sum = 0.0;
    for (int i_spec_k = 0; i_spec_k <= 7; ++i_spec_k)
      sum = sum + (h[net_nonzero_spec_no_ns[i_spec_k]] + -1.0 * h[8] * mw_factor[net_nonzero_spec_no_ns[i_spec_k]]) * jac[10 * (net_nonzero_spec_no_ns[i_spec_k] + 2) + 2 + i];
    jac[2 + i] = -1.0 * (sum + dphi[0] * (cp[i] + -1.0 * cp[8])) / (phi[1] * cp_tot[0]);
  }
}
void dVdot_dnj(double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ jac)
{
  double T_inv;
  double T_val;
  double sum;

  {
    int const j = 0;

    T_val = phi[0];
    T_inv = 1.0 / phi[0];
    for (int i = 0; i <= 7; ++i)
    {
      sum = 0.0;
      for (int i_spec_k = 0; i_spec_k <= 7; ++i_spec_k)
        sum = sum + (1.0 + -1.0 * mw_factor[net_nonzero_spec_no_ns[i_spec_k]]) * jac[10 * (net_nonzero_spec_no_ns[i_spec_k] + 2) + 2 + i];
      jac[12 + i] = jac[12 + i] + T_val * 8314.4621 * sum / P_arr[0] + phi[1] * jac[2 + i] * T_inv;
    }
  }
}
void dRopi_dT(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ db, double *__restrict__ jac)
{
  double Tinv;
  double ci;
  double dBk_sum;
  double dRopidT;
  double dkf;
  int i_0;
  int offset;
  int offset_next;

  {
    int const j = 0;

    Tinv = 1.0 / phi[0];
    for (int i = 0; i <= 27; ++i)
    {
      offset_next = net_reac_to_spec_offsets[1 + i];
      offset = net_reac_to_spec_offsets[i];
      i_0 = thd_mask[i];
      dkf = (simple_beta[i] + simple_Ta[i] * Tinv) * Tinv;
      dRopidT = rop_fwd[i] * dkf;
      if (rev_mask[i] >= 0)
        dBk_sum = 0.0;
      ci = 1.0;
      if (thd_mask[i] >= 0)
        ci = pres_mod[i_0];
      if (rev_mask[i] >= 0)
      {
        for (int net_ind = offset; net_ind <= -1 + offset_next; ++net_ind)
          dBk_sum = dBk_sum + (reac_to_spec_nu[2 * net_ind] + -1.0 * reac_to_spec_nu[1 + 2 * net_ind]) * db[rxn_to_spec[net_ind]];
        dRopidT = dRopidT + -1.0 * rop_rev[i] * (dkf + -1.0 * dBk_sum);
      }
      dRopidT = dRopidT * ci * phi[1];
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2)] = jac[10 * (rxn_to_spec[k_ind] + 2)] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dRopidT;
    }
  }
}
void dRopi_dT_ns(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ P_arr, double const *__restrict__ conc, double *__restrict__ jac)
{
  double Sns_fwd;
  double Sns_rev;
  double ci;
  double dRopidT;
  int i_1;
  int i_2;
  double kr_i;
  int net_spec;
  int nu_fwd;
  int nu_rev;
  int offset;
  int offset_next;

  {
    int const i = 0;

    {
      int const j = 0;

      kr_i = 0.0;
      offset_next = net_reac_to_spec_offsets[9];
      offset = net_reac_to_spec_offsets[8];
      i_2 = rev_mask[8];
      if (rev_mask[8] >= 0)
        kr_i = kr[i_2];
      i_1 = thd_mask[8];
      ci = 1.0;
      if (thd_mask[8] >= 0)
        ci = pres_mod[i_1];
      Sns_rev = reac_to_spec_nu[offset_next + offset_next + -2];
      Sns_fwd = reac_to_spec_nu[offset_next + offset_next + -2 + 1];
      for (int net_ind = offset; net_ind <= -1 + offset_next; ++net_ind)
      {
        nu_rev = reac_to_spec_nu[2 * net_ind];
        nu_fwd = reac_to_spec_nu[1 + 2 * net_ind];
        net_spec = rxn_to_spec[net_ind];
        if (net_spec == 8)
        {
          nu_rev = nu_rev + -1;
          nu_fwd = nu_fwd + -1;
        }
        Sns_rev = Sns_rev * fast_powi(fmax(1e-300, conc[net_spec]), nu_rev);
        Sns_fwd = Sns_fwd * fast_powi(fmax(1e-300, conc[net_spec]), nu_fwd);
      }
      dRopidT = (Sns_rev * kr_i + -1.0 * Sns_fwd * kf[8]) * phi[1] * ci * P_arr[0] / (8314.4621 * phi[0] * phi[0]);
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2)] = jac[10 * (rxn_to_spec[k_ind] + 2)] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dRopidT;
    }
  }
}
void dci_thd_dT(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ jac)
{
  double Tinv;
  double dci_thd_dT_fac;
  int i_0;
  int i_1;
  double mod;
  int offset;
  int offset_next;
  double rop_net;

  {
    int const j = 0;

    Tinv = 1.0 / phi[0];
    for (int i = 0; i <= 4; ++i)
    {
      mod = 1.0;
      if (thd_type[i] == 2)
        mod = thd_spec[thd_offset[1 + i] + -1] == 8;
      if (thd_type[i] == 1 && thd_spec[thd_offset[1 + i] + -1] == 8)
        mod = thd_eff[thd_offset[1 + i] + -1];
      i_0 = thd_map[i];
      rop_net = rop_fwd[i_0];
      offset_next = net_reac_to_spec_offsets[i_0 + 1];
      offset = net_reac_to_spec_offsets[i_0];
      i_1 = rev_mask[i_0];
      if (i_1 >= 0)
        rop_net = rop_net + -1.0 * rop_rev[i_1];
      dci_thd_dT_fac = -1.0 * P_arr[0] * mod * 0.00012027236253804078 * Tinv * Tinv * phi[1] * rop_net;
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2)] = jac[10 * (rxn_to_spec[k_ind] + 2)] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dci_thd_dT_fac;
    }
  }
}
void dci_troe_dT(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ pres_mod, double const *__restrict__ Fi, double const *__restrict__ Pr, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac)
{
  double Ta_0;
  double Ta_inf;
  double Tinv;
  double Tval;
  double absq;
  double absqsq;
  double beta_0;
  double beta_inf;
  double dFcent;
  double dFi;
  double dci_fall_dT;
  double dci_thd_dT_fac;
  int i_1;
  int i_2;
  int i_3;
  double kf_0;
  double kf_inf;
  double logFcent;
  double mod;
  int offset;
  int offset_next;
  double pmod;
  double rop_net;
  double theta_Pr;
  double theta_no_Pr;

  {
    int const j = 0;

    Tval = phi[0];
    Tinv = 1.0 / phi[0];
    {
      int const i = 0;

      mod = 1.0;
      if (thd_type[5] == 2)
        mod = thd_spec[thd_offset[6] + -1] == 8;
      if (thd_type[5] == 1 && thd_spec[thd_offset[6] + -1] == 8)
        mod = thd_eff[thd_offset[6] + -1];
      if (fall_type[0])
        kf_inf = kf_fall[0];
      if (!fall_type[0])
        kf_0 = kf_fall[0];
      logFcent = log(fmax(1e-300, Fcent[0]));
      dFcent = -1.0 * troe_a[0] * troe_T1[0] * exp(fmin(690.775527898, -1.0 * Tval * troe_T1[0])) + (troe_a[0] + -1.0) * troe_T3[0] * exp(fmin(690.775527898, -1.0 * Tval * troe_T3[0])) + troe_T2[0] * Tinv * Tinv * exp(fmin(690.775527898, -1.0 * troe_T2[0] * Tinv));
      pmod = pres_mod[5];
      i_1 = thd_map[5];
      rop_net = rop_fwd[i_1];
      if (!fall_type[0])
        kf_inf = kf[i_1];
      if (fall_type[0])
        kf_0 = kf[i_1];
      offset_next = net_reac_to_spec_offsets[i_1 + 1];
      offset = net_reac_to_spec_offsets[i_1];
      i_3 = rev_mask[i_1];
      if (i_3 >= 0)
        rop_net = rop_net + -1.0 * rop_rev[i_3];
      dci_thd_dT_fac = -1.0 * P_arr[0] * mod * 0.00012027236253804078 * Tinv * Tinv;
      theta_no_Pr = dci_thd_dT_fac * kf_0 / kf_inf;
      i_2 = simple_mask[i_1];
      if (!fall_type[0])
        beta_inf = simple_beta[i_2];
      if (fall_type[0])
        beta_inf = fall_beta[0];
      if (!fall_type[0])
        beta_0 = fall_beta[0];
      if (fall_type[0])
        beta_0 = simple_beta[i_2];
      absq = Atroe[0] * Atroe[0] + Btroe[0] * Btroe[0];
      absqsq = absq * absq;
      if (!fall_type[0])
        Ta_inf = simple_Ta[i_2];
      if (fall_type[0])
        Ta_inf = fall_Ta[0];
      if (!fall_type[0])
        Ta_0 = fall_Ta[0];
      if (fall_type[0])
        Ta_0 = simple_Ta[i_2];
      theta_Pr = Tinv * (beta_0 + -1.0 * beta_inf + (Ta_0 + -1.0 * Ta_inf) * Tinv);
      dFi = -1.0 * Btroe[0] * (2.0 * Atroe[0] * Fcent[0] * (0.14 * Atroe[0] + Btroe[0]) * (Pr[0] * theta_Pr + theta_no_Pr) * logFcent + Pr[0] * dFcent * (2.0 * Atroe[0] * (1.1762 * Atroe[0] + -1.0 * 0.67 * Btroe[0]) * logFcent + -1.0 * Btroe[0] * absq * 2.302585092994046)) / (Fcent[0] * fmax(1e-300, Pr[0]) * absqsq * 2.302585092994046);
      dci_fall_dT = pmod * (-1.0 * (Pr[0] * theta_Pr + theta_no_Pr) / (Pr[0] + 1.0) + dFi);
      if (!fall_type[0])
        dci_fall_dT = dci_fall_dT + theta_Pr * pmod + Fi[0] * theta_no_Pr / (Pr[0] + 1.0);
      dci_fall_dT = dci_fall_dT * phi[1] * rop_net;
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2)] = jac[10 * (rxn_to_spec[k_ind] + 2)] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dci_fall_dT;
    }
  }
}
void dTdot_dT(double const *__restrict__ cp_tot, double const *__restrict__ dphi, double const *__restrict__ cp, double const *__restrict__ dcp, double const *__restrict__ h, double const *__restrict__ conc, double const *__restrict__ phi, double const *__restrict__ wdot, double *__restrict__ jac)
{
  double Tinv;
  double Vinv;
  double dTsum;
  double rate_sum;

  {
    int const j = 0;

    rate_sum = 0.0;
    Tinv = 1.0 / phi[0];
    dTsum = (cp[8] * Tinv + -1.0 * dcp[8]) * conc[8];
    Vinv = 1.0 / phi[1];
    for (int i = 0; i <= 7; ++i)
    {
      rate_sum = rate_sum + Vinv * jac[20 + 10 * i] * (-1.0 * h[i] + h[8] * mw_factor[i]);
      rate_sum = rate_sum + wdot[i] * (-1.0 * cp[i] + mw_factor[i] * cp[8]);
      dTsum = dTsum + (cp[8] * Tinv + -1.0 * dcp[i]) * conc[i];
    }
    jac[0] = jac[0] + (dphi[0] * dTsum + rate_sum) / cp_tot[0];
  }
}
void dVdotdT(double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ wdot, double const *__restrict__ dphi, double *__restrict__ jac)
{
  double Tinv;
  double Vinv;
  double sum;

  {
    int const j = 0;

    Vinv = 1.0 / phi[1];
    Tinv = 1.0 / phi[0];
    sum = 0.0;
    for (int i = 0; i <= 7; ++i)
      sum = sum + (1.0 + -1.0 * mw_factor[i]) * (Vinv * jac[20 + 10 * i] + Tinv * wdot[i]);
    jac[10] = jac[10] + phi[1] * Tinv * (jac[0] + -1.0 * Tinv * dphi[0]);
    jac[10] = jac[10] + 8314.4621 * phi[0] * phi[1] * sum / P_arr[0];
  }
}
void dRopi_dV(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double *__restrict__ jac)
{
  double dRopi_dE;
  int i_0;
  int nu_fwd;
  int nu_rev;
  int offset;
  int offset_next;

  for (int i = 0; i <= 27; ++i)
  {
    int const j = 0;

    nu_rev = -1;
    nu_fwd = -1;
    offset_next = net_reac_to_spec_offsets[1 + i];
    offset = net_reac_to_spec_offsets[i];
    i_0 = thd_mask[i];
    for (int net_ind = offset; net_ind <= -1 + offset_next; ++net_ind)
    {
      nu_rev = nu_rev + reac_to_spec_nu[2 * net_ind];
      nu_fwd = nu_fwd + reac_to_spec_nu[1 + 2 * net_ind];
    }
    dRopi_dE = -1.0 * nu_fwd * rop_fwd[i];
    dRopi_dE = dRopi_dE + nu_rev * rop_rev[i];
    if (i_0 >= 0)
      dRopi_dE = dRopi_dE * pres_mod[i_0];
    for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
      if (rxn_to_spec[k_ind] != 8)
        jac[10 * (rxn_to_spec[k_ind] + 2) + 1] = jac[10 * (rxn_to_spec[k_ind] + 2) + 1] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dRopi_dE;
  }
}
void dRopi_dV_ns(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ P_arr, double const *__restrict__ conc, double *__restrict__ jac)
{
  double Sns_fwd;
  double Sns_rev;
  double ci;
  double dRopi_dE;
  double fac;
  int i_1;
  int i_2;
  double kr_i;
  int net_spec;
  int nu_fwd;
  int nu_rev;
  int offset;
  int offset_next;

  {
    int const j = 0;

    fac = P_arr[0] / (8314.4621 * phi[0]);
    {
      int const i = 0;

      kr_i = 0.0;
      offset_next = net_reac_to_spec_offsets[9];
      offset = net_reac_to_spec_offsets[8];
      i_2 = rev_mask[8];
      if (i_2 >= 0)
        kr_i = kr[i_2];
      i_1 = thd_mask[8];
      ci = 1.0;
      if (i_1 >= 0)
        ci = pres_mod[i_1];
      Sns_rev = reac_to_spec_nu[offset_next + offset_next + -2];
      Sns_fwd = reac_to_spec_nu[offset_next + offset_next + -2 + 1];
      for (int net_ind = offset; net_ind <= -1 + offset_next; ++net_ind)
      {
        nu_rev = reac_to_spec_nu[2 * net_ind];
        nu_fwd = reac_to_spec_nu[1 + 2 * net_ind];
        net_spec = rxn_to_spec[net_ind];
        if (net_spec == 8)
        {
          nu_rev = nu_rev + -1;
          nu_fwd = nu_fwd + -1;
        }
        Sns_rev = Sns_rev * fast_powi(fmax(1e-300, conc[net_spec]), nu_rev);
        Sns_fwd = Sns_fwd * fast_powi(fmax(1e-300, conc[net_spec]), nu_fwd);
      }
      dRopi_dE = (Sns_fwd * kf[8] + -1.0 * Sns_rev * kr_i) * ci * fac;
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2) + 1] = jac[10 * (rxn_to_spec[k_ind] + 2) + 1] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dRopi_dE;
    }
  }
}
void dci_thd_dE(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ pres_mod, double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ jac)
{
  double dci_thd_dE_fac;
  int i_0;
  int i_1;
  double mod;
  int offset;
  int offset_next;
  double rop_net;
  double rt_inv;

  {
    int const j = 0;

    rt_inv = 1.0 / (8314.4621 * phi[0]);
    for (int i = 0; i <= 4; ++i)
    {
      mod = thd_type[i] == 1;
      if (thd_type[i] == 2)
        mod = thd_spec[thd_offset[1 + i] + -1] == 8;
      if (thd_type[i] == 1 && thd_spec[thd_offset[1 + i] + -1] == 8)
        mod = thd_eff[thd_offset[1 + i] + -1];
      if (thd_type[i] != 3)
        mod = mod * P_arr[0] * rt_inv + -1.0 * pres_mod[i];
      i_0 = thd_map[i];
      rop_net = rop_fwd[i_0];
      offset_next = net_reac_to_spec_offsets[i_0 + 1];
      offset = net_reac_to_spec_offsets[i_0];
      i_1 = rev_mask[i_0];
      if (i_1 >= 0)
        rop_net = rop_net + -1.0 * rop_rev[i_1];
      dci_thd_dE_fac = mod * rop_net;
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2) + 1] = jac[10 * (rxn_to_spec[k_ind] + 2) + 1] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dci_thd_dE_fac;
    }
  }
}
void dci_troe_dE(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ pres_mod, double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ Fi, double const *__restrict__ Pr, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac)
{
  double absqsq;
  double dFi;
  double dci_fall_dE;
  int i_1;
  int i_2;
  double kf_0;
  double kf_inf;
  double mod;
  int not_unity;
  int offset;
  int offset_next;
  double rop_net;
  double rt_inv;

  {
    int const j = 0;

    rt_inv = 1.0 / (8314.4621 * phi[0]);
    {
      int const i = 0;

      mod = thd_type[5] == 1;
      if (thd_type[5] == 2)
        mod = thd_spec[thd_offset[6] + -1] == 8;
      if (thd_type[5] == 1 && thd_spec[thd_offset[6] + -1] == 8)
        mod = thd_eff[thd_offset[6] + -1];
      if (fall_type[0])
        kf_inf = kf_fall[0];
      if (!fall_type[0])
        kf_0 = kf_fall[0];
      not_unity = thd_type[5] != 3;
      i_1 = thd_map[5];
      rop_net = rop_fwd[i_1];
      if (!fall_type[0])
        kf_inf = kf[i_1];
      if (fall_type[0])
        kf_0 = kf[i_1];
      mod = mod * P_arr[0] * rt_inv * kf_0 / kf_inf;
      offset_next = net_reac_to_spec_offsets[i_1 + 1];
      offset = net_reac_to_spec_offsets[i_1];
      i_2 = rev_mask[i_1];
      if (i_2 >= 0)
        rop_net = rop_net + -1.0 * rop_rev[i_2];
      absqsq = Atroe[0] * Atroe[0] + Btroe[0] * Btroe[0];
      absqsq = absqsq * absqsq;
      dFi = -2.0 * Atroe[0] * Btroe[0] * log(fmax(1e-300, Fcent[0])) * (0.14 * Atroe[0] + Btroe[0]) * (mod + -1.0 * Pr[0] * not_unity) / (fmax(1e-300, Pr[0]) * absqsq * 2.302585092994046);
      dci_fall_dE = pres_mod[5] * ((-1.0 * mod + Pr[0] * not_unity) / (Pr[0] + 1.0) + dFi);
      if (!fall_type[0])
        dci_fall_dE = dci_fall_dE + Fi[0] * mod / (Pr[0] + 1.0) + -1.0 * not_unity * pres_mod[5];
      dci_fall_dE = dci_fall_dE * rop_net;
      for (int k_ind = offset; k_ind <= -1 + offset_next; ++k_ind)
        if (rxn_to_spec[k_ind] != 8)
          jac[10 * (rxn_to_spec[k_ind] + 2) + 1] = jac[10 * (rxn_to_spec[k_ind] + 2) + 1] + (reac_to_spec_nu[2 * k_ind] + -1.0 * reac_to_spec_nu[1 + 2 * k_ind]) * dci_fall_dE;
    }
  }
}
void dTdotdV(double const *__restrict__ cp_tot, double const *__restrict__ cp, double const *__restrict__ h, double const *__restrict__ wdot, double const *__restrict__ dphi, double *__restrict__ jac, double const *__restrict__ conc, double const *__restrict__ phi)
{
  double dTsum;
  double spec_inv;
  double specsum;

  {
    int const j = 0;

    specsum = 0.0;
    spec_inv = 1.0 / (cp_tot[0] * phi[1]);
    dTsum = 0.0;
    for (int i = 0; i <= 7; ++i)
    {
      dTsum = dTsum + (cp[i] + -1.0 * cp[8]) * conc[i];
      specsum = specsum + (h[i] + -1.0 * h[8] * mw_factor[i]) * (jac[21 + 10 * i] + -1.0 * wdot[i]);
    }
    jac[1] = jac[1] + (dphi[0] * dTsum + -1.0 * specsum) * spec_inv;
  }
}
void dVdotdV(double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ dphi, double *__restrict__ jac)
{
  double sum;

  {
    int const j = 0;

    sum = 0.0;
    for (int i = 0; i <= 7; ++i)
      sum = sum + (1.0 + -1.0 * mw_factor[i]) * jac[21 + 10 * i];
    jac[11] = jac[11] + (phi[1] * jac[1] + dphi[0]) / phi[0];
    jac[11] = jac[11] + 8314.4621 * phi[0] * sum / P_arr[0];
  }
}

void jacobian(double const *__restrict__ t, double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ jac, double *__restrict__ rwk)
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
    double* __restrict__ cp_tot = rwk + 60 * work_size + 1 * omp_get_thread_num();
    double* __restrict__ db = rwk + 61 * work_size + 9 * omp_get_thread_num();
    double* __restrict__ dcp = rwk + 70 * work_size + 9 * omp_get_thread_num();
    double* __restrict__ dphi = rwk + 79 * work_size + 10 * omp_get_thread_num();
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
    
    reset_arrays(jac);
    
    dRopidnj(pres_mod, kf, kr, conc, jac);
    
    dRopidnj_ns(pres_mod, kf, kr, conc, jac);
    
    dci_thd_dnj(rop_fwd, rop_rev, jac);
    
    dci_thd_dnj_ns(rop_fwd, rop_rev, jac);
    
    dci_troe_dnj(rop_fwd, rop_rev, Pr, pres_mod, kf, kf_fall, Fi, Atroe, Btroe, Fcent, jac);
    
    dci_troe_dnj_ns(rop_fwd, rop_rev, Pr, pres_mod, kf, kf_fall, Fi, Atroe, Btroe, Fcent, jac);
    
    cp_total(cp, conc, cp_tot);
    
    eval_dcp(phi, dcp);
    
    eval_db(phi, db);
    
    dTdot_dnj(cp, h, cp_tot, phi, dphi, jac);
    
    dVdot_dnj(phi, P_arr, jac);
    
    dRopi_dT(phi, pres_mod, rop_fwd, rop_rev, db, jac);
    
    dRopi_dT_ns(phi, pres_mod, kf, kr, P_arr, conc, jac);
    
    dci_thd_dT(rop_fwd, rop_rev, phi, P_arr, jac);
    
    dci_troe_dT(rop_fwd, rop_rev, phi, P_arr, pres_mod, Fi, Pr, kf, kf_fall, Atroe, Btroe, Fcent, jac);
    
    dTdot_dT(cp_tot, dphi, cp, dcp, h, conc, phi, wdot, jac);
    
    dVdotdT(phi, P_arr, wdot, dphi, jac);
    
    dRopi_dV(phi, pres_mod, rop_fwd, rop_rev, jac);
    
    dRopi_dV_ns(phi, pres_mod, kf, kr, P_arr, conc, jac);
    
    dci_thd_dE(rop_fwd, rop_rev, pres_mod, phi, P_arr, jac);
    
    dci_troe_dE(rop_fwd, rop_rev, pres_mod, phi, P_arr, Fi, Pr, kf, kf_fall, Atroe, Btroe, Fcent, jac);
    
    dTdotdV(cp_tot, cp, h, wdot, dphi, jac, conc, phi);
    
    dVdotdV(phi, P_arr, dphi, jac);
}
