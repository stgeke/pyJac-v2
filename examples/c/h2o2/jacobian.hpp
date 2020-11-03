#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP
#ifdef _OPENMP
 #include <omp.h>
#else
 #warning 'OpenMP not found! Unexpected results may occur if using more than one thread.'
 #define omp_get_num_threads() (1)
#endif
#include "mechanism.hpp"
#include "jacobian.hpp"
#include "species_rates.hpp"
#include "chem_utils.hpp"
void reset_arrays(double *__restrict__ jac);
void dRopidnj(double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ conc, double *__restrict__ jac);
void dRopidnj_ns(double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ conc, double *__restrict__ jac);
void dci_thd_dnj(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double *__restrict__ jac);
void dci_thd_dnj_ns(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double *__restrict__ jac);
void dci_troe_dnj(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ Pr, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Fi, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac);
void dci_troe_dnj_ns(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ Pr, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Fi, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac);
void cp_total(double const *__restrict__ cp, double const *__restrict__ conc, double *__restrict__ cp_tot);
void eval_dcp(double const *__restrict__ phi, double *__restrict__ dcp);
void eval_db(double const *__restrict__ phi, double *__restrict__ db);
void dTdot_dnj(double const *__restrict__ cp, double const *__restrict__ h, double const *__restrict__ cp_tot, double const *__restrict__ phi, double const *__restrict__ dphi, double *__restrict__ jac);
void dVdot_dnj(double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ jac);
void dRopi_dT(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ db, double *__restrict__ jac);
void dRopi_dT_ns(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ P_arr, double const *__restrict__ conc, double *__restrict__ jac);
void dci_thd_dT(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ jac);
void dci_troe_dT(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ pres_mod, double const *__restrict__ Fi, double const *__restrict__ Pr, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac);
void dTdot_dT(double const *__restrict__ cp_tot, double const *__restrict__ dphi, double const *__restrict__ cp, double const *__restrict__ dcp, double const *__restrict__ h, double const *__restrict__ conc, double const *__restrict__ phi, double const *__restrict__ wdot, double *__restrict__ jac);
void dVdotdT(double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ wdot, double const *__restrict__ dphi, double *__restrict__ jac);
void dRopi_dV(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double *__restrict__ jac);
void dRopi_dV_ns(double const *__restrict__ phi, double const *__restrict__ pres_mod, double const *__restrict__ kf, double const *__restrict__ kr, double const *__restrict__ P_arr, double const *__restrict__ conc, double *__restrict__ jac);
void dci_thd_dE(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ pres_mod, double const *__restrict__ phi, double const *__restrict__ P_arr, double *__restrict__ jac);
void dci_troe_dE(double const *__restrict__ rop_fwd, double const *__restrict__ rop_rev, double const *__restrict__ pres_mod, double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ Fi, double const *__restrict__ Pr, double const *__restrict__ kf, double const *__restrict__ kf_fall, double const *__restrict__ Atroe, double const *__restrict__ Btroe, double const *__restrict__ Fcent, double *__restrict__ jac);
void dTdotdV(double const *__restrict__ cp_tot, double const *__restrict__ cp, double const *__restrict__ h, double const *__restrict__ wdot, double const *__restrict__ dphi, double *__restrict__ jac, double const *__restrict__ conc, double const *__restrict__ phi);
void dVdotdV(double const *__restrict__ phi, double const *__restrict__ P_arr, double const *__restrict__ dphi, double *__restrict__ jac);
static int const net_nonzero_spec_no_ns[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static int const rev_mask[28] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
static int const simple_mask[28] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
static int const thd_map[6] = { 0, 1, 5, 10, 13, 20 };
void jacobian(double const *__restrict__ t, double const *__restrict__ P_arr, double const *__restrict__ phi, double *__restrict__ jac, double *__restrict__ rwk);
#endif