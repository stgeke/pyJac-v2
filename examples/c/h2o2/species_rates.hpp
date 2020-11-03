#ifndef SPECIES_RATES_HPP
#define SPECIES_RATES_HPP
#ifdef _OPENMP
 #include <omp.h>
#else
 #warning 'OpenMP not found! Unexpected results may occur if using more than one thread.'
 #define omp_get_num_threads() (1)
#endif
#include "mechanism.hpp"
#include "species_rates.hpp"
#include "chem_utils.hpp"
void ndot_reset(double *__restrict__ dphi);
void wdot_reset(double *__restrict__ wdot);
void get_concentrations(double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ conc);
void a_only_simple(double const *__restrict__ phi, double *__restrict__ kf);

static inline double fast_powi(double val, int pow)
{
     // account for negatives
     if (pow < 0)
     {
         val = 1.0 / val;
         pow = -pow;
     }
     // switch for speed
     switch(pow)
     {
         case 0:
             return 1;
         case 1:
             return val;
         case 2:
             return val * val;
         case 3:
             return val * val * val;
         case 4:
             return val * val * val * val;
         case 5:
             return val * val * val * val * val;
     }
     double retval = val * val * val * val * val * val;
     for (int i = 6; i < pow; ++i)
     {
         retval *= val;
     }
     return retval;
}

void beta_int_simple(double const *__restrict__ phi, double *__restrict__ kf);
void rateconst_fullsimple(double const *__restrict__ phi, double *__restrict__ kf);
void eval_thd_body_concs(double const *__restrict__ P_arr, double const *__restrict__ phi, double const *__restrict__ conc, double *__restrict__ thd_conc);
void rateconst_fullfall(double const *__restrict__ phi, double *__restrict__ kf_fall);
void red_pres(double const *__restrict__ phi, double const *__restrict__ thd_conc, double const *__restrict__ kf, double const *__restrict__ kf_fall, double *__restrict__ Pr);
void fall_troe(double const *__restrict__ Pr, double const *__restrict__ phi, double *__restrict__ Fi, double *__restrict__ Fcent, double *__restrict__ Atroe, double *__restrict__ Btroe);
void rateconst_Kc(double const *__restrict__ b, double *__restrict__ Kc, double const *__restrict__ kf, double *__restrict__ kr);
void ci_thd(double const *__restrict__ thd_conc, double *__restrict__ pres_mod);
void ci_fall(double const *__restrict__ Fi, double const *__restrict__ Pr, double *__restrict__ pres_mod);
void rop_eval_fwd(double const *__restrict__ conc, double const *__restrict__ kf, double *__restrict__ rop_fwd);
void rop_eval_rev(double const *__restrict__ conc, double const *__restrict__ kr, double *__restrict__ rop_rev);
void rop_net_fixed(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ pres_mod, double *__restrict__ rop_net);
void spec_rates(double const *__restrict__ rop_net, double *__restrict__ wdot);
void get_molar_rates(double const *__restrict__ phi, double *__restrict__ dphi, double const *__restrict__ wdot);
void temperature_rate(double const *__restrict__ h, double const *__restrict__ cp, double const *__restrict__ conc, double *__restrict__ dphi, double const *__restrict__ wdot);
void get_extra_var_rates(double const *__restrict__ wdot, double *__restrict__ dphi, double const *__restrict__ phi, double const *__restrict__ P_arr);
static double const fall_A[1] = { 28.463930238863654 };
static double const fall_Ta[1] = { -855.4732602605767 };
static double const fall_beta[1] = { -0.9 };
static double const mw_factor[8] = { 0.05046260138179634, 0.02523130069089817, 0.4005056573545609, 0.8010113147091218, 0.4257369580454591, 0.45096825873635726, 0.8262426154000201, 0.8514739160909182 };
static double const simple_A[28] = { 120000000000.0, 500000000000.0, 3.655839600035736, 20000000000.0, 9.172638504792172, 28.660640533109706, 30.665974102635822, 30.052277738640093, 27.274346171989816, 30.908165848920724, 1000000000000.0, 25.223075507276675, 31.72536567815065, 2.2e+16, 22.10203193164551, 24.52547397636735, 25.154082635789724, 9.400960731584833, 23.025850929940457, 12.283033686666302, 25.02733093015058, 3.5751506887855933, 23.39741448637294, 21.416413017506358, 35.06940464597285, 18.683045008419857, 26.763520548223823, 29.240459028362647 };
static double const simple_Ta[28] = { 0.0, 0.0, 3150.1544760183588, 0.0, 2012.878259436651, 0.0, 0.0, 0.0, 0.0, 8575.364604764993, 0.0, 0.0, 0.0, 0.0, 337.6603280204982, 537.4384952695858, 319.54442368556835, 2616.7417372676464, 1811.590433492986, 1726.0431074669282, 0.0, -1061.7932818528334, -251.60978242958137, 214.8747541948625, 14799.687402507976, -820.2478907204353, 6038.634778309953, 8720.79505900929 };
static double const simple_beta[28] = { -1.0, -1.0, 2.7, 0.0, 2.0, -0.86, -1.24, -0.76, -0.8, -0.6707, -1.0, -0.6, -1.25, -2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.51, -0.37, 2.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static double const thd_eff[18] = { 2.4, 15.4, 0.83, 2.0, 6.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.63, 0.73, 3.65, 0.38, 2.0, 6.0, 0.7 };
static double const troe_T1[1] = { 0.0005694760820045558 };
static double const troe_T2[1] = { 5182.0 };
static double const troe_T3[1] = { 0.010638297872340425 };
static double const troe_a[1] = { 0.7346 };
static int const fall_type[1] = { 0 };
static int const net_reac_to_spec_offsets[29] = { 0, 2, 5, 9, 13, 17, 20, 23, 27, 31, 35, 37, 39, 42, 45, 49, 53, 56, 60, 64, 68, 70, 73, 77, 81, 85, 88, 91, 95 };
static int const nu_sum[28] = { -1, -1, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0 };
static int const reac_to_spec_nu[190] = { 0, 2, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 2, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 2, 2, 1, 0, 2, 1, 0, 0, 2, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 1, 0, 1, 0, 0, 2, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 2, 1, 0, 1, 0, 0, 2, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1 };
static int const rxn_to_spec[95] = { 2, 3, 1, 2, 4, 0, 1, 2, 4, 2, 3, 4, 6, 2, 4, 6, 7, 1, 3, 6, 1, 3, 6, 1, 3, 5, 6, 1, 3, 6, 8, 1, 2, 3, 4, 0, 1, 0, 1, 0, 1, 5, 1, 4, 5, 1, 2, 5, 6, 0, 1, 3, 6, 1, 4, 6, 0, 1, 6, 7, 1, 4, 5, 7, 0, 1, 4, 5, 4, 7, 2, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 4, 5, 6, 7, 3, 6, 7, 3, 6, 7, 3, 4, 5, 6 };
static int const simple_rtype_1_inds[4] = { 0, 1, 10, 13 };
static int const simple_rtype_1_map[4] = { 0, 1, 10, 13 };
static int const simple_rtype_2_inds[23] = { 2, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
static int const simple_rtype_2_map[23] = { 2, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
static int const thd_mask[28] = { 0, 1, -1, -1, -1, 2, -1, -1, -1, -1, 3, -1, -1, 4, -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1 };
static int const thd_offset[7] = { 0, 3, 6, 9, 12, 15, 18 };
static int const thd_spec[18] = { 0, 5, 8, 0, 5, 8, 3, 5, 8, 0, 5, 8, 0, 5, 8, 0, 5, 8 };
static int const thd_type[6] = { 1, 1, 1, 1, 1, 1 };
void species_rates(double const *__restrict__ t, double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ dphi, double *__restrict__ rwk);
#endif