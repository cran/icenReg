//
//  intCoxFast.h
//  
//
//  Created by Cliff Anderson Bergman on 4/10/15.
//
//

#ifndef _intCoxFast_h
#define _intCoxFast_h

#include "Eigen_local/Dense"
#include <stdio.h>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
using namespace std;



class NN_info{
public:
    int nn_id;
    vector<int> g1, g2; //NOTE TO SELF: MAKE SURE THESE ARE UPDATED WHENEVER NN'S ARE SWITCHED!!!!
};

class VEM_info{
public:
    vector<int> g_r;
    vector<bool> g_l;       //this will be of same length as g_r.
                            //It corresponds to whether the p_mass is also less l[g_r], given that it must be less than r to be in g_r
};

class ICPH_OptimInfo{
public:
    // Probability mass objects
    vector<double> p_mass;  //vector of probability masses
    vector<double> S;   // survival function from probability masses. This is the baseline survival

    //observation objects
    vector<int> lind;   //indices of left side of each observation
    vector<int> rind;   //indices of right side of each observation
    Eigen::VectorXd eta; //linear function of covariates
    vector<double> expEta;  //exponentiated linear function of covariates
    vector<double> P_obs;   //probability of each observation: i.e. pow(S[lind[i]], expEta[i]) - pow(S[rind[i], expEta[i])
    
    vector<double> d_ob_cont;    //contribution each observation makes to derivative of likelihood function
    
    //nearest neighbor objects
    vector<NN_info> nn_info;
    vector<VEM_info> vem_info;
    vector<double> nn_ders;
    vector<double> vdm_ders;
    double h;
    int hSteps;
    
    Eigen::VectorXd propVec;
    Eigen::VectorXd betas;
    Eigen::MatrixXd covars;
    
    Eigen::VectorXd d_cov;
    Eigen::MatrixXd d2_mat;
    
    ICPH_OptimInfo(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars);
};


//void update_p_ob(int i, ICPH_OptimInfo &optinfo);
//calculates p_ob for i'th observation

//double icph_llk(vector<double> &S, vector<double> &expEta, vector<int> &lind, vector<int> &rind);
//calculates llk by summing up all P_obs
// Does not update any of the vectors

//double icph_llk(ICPH_OptimInfo &optinfo){return(icph_llk(optinfo.S, optinfo.expEta, optinfo.lind, optinfo.rind));}
//wrapper for above function, but plugs in the ICPH_optimInfo fields

double icph_llk_sumOnly(ICPH_OptimInfo &optinfo);
// Faster than the above functions: just sums up the already
// computed P_obs

double min_icph_llk(vector<double> &P_obs, vector<int> &g1, vector<int> &g2);
//same as above, but only uses the observations in g1 and g2

//void getNecInd(int p_ind1, int p_ind2, vector<int> &l_inds, vector<int> &r_inds, vector<int> &g1, vector<int> &g2);
//given indices p_ind1 and p_ind2,
//finds necessary g1 and g2

//void NNE_exchange(int p1_ind, int p2_ind, double alpha, ICPH_OptimInfo &optInfo);
//exchanges probability mass between p1_ind and p2_ind
//updates S, P_obs, p_mass accordingly

void check_NN_id(int id1, int propID, ICPH_OptimInfo &optInfo);
//checks that propID is listed as nearest neighbor to id1.
//if false, updates NN_info so that it is

void update_p_ob(int i, ICPH_OptimInfo &optinfo);
// updates the P_obs vector: requires S and expEta being up to date

//void numericDerv_NNE(int p1_ind, int p2_ind, ICPH_OptimInfo &optinfo);
// numeric dervs for nearest neighbor exchange

//void analyticDerv_NNE(int p1_ind, int p2_ind, ICPH_OptimInfo &optinfo);
// analytic dervs for nne, not used

void NNE_optim(int p1, int p2, ICPH_OptimInfo &optinfo);
//optimizes a NN exchange. Currently works with numeric dervs

void NNE_all(ICPH_OptimInfo &optinfo);
//Does all the NN optimizing exchanges for points with positive mass

void analyticDerv_VEM(vector<int> &minmax, vector<double> &output_dervs, ICPH_OptimInfo  &optinfo);
// Analtyic derivatives for VEM step

//void vem_update(ICPH_OptimInfo &optinfo);
void vem_update(ICPH_OptimInfo &optinfo, bool recPosd);
// VEM update

double max(double a, double b);
double min(double a, double b);

void vec_delta(vector<double> &delta, ICPH_OptimInfo &optinfo);
// This alters p_mass by delta and updates S and p_mass

void mult_vec(double a, vector<double> &vec);
// multiplies vec by a

void analyticDerv_ICM(vector<double> &d1, vector<double> &d2, ICPH_OptimInfo  &optinfo);
// computes the derivatives for the ICM step

void ICM_step(ICPH_OptimInfo &optinfo);
// computes the ICM step

void update_etas(ICPH_OptimInfo &optinfo);
// updates the etas and expEta

void covar_d_analytic(Eigen::VectorXd &d1, Eigen::MatrixXd &dmat, ICPH_OptimInfo &optinfo);
// computes the analytic derivatives for covariates

void updateCovars(ICPH_OptimInfo &optinfo);
// updates regression parameters



double signVal(double x){
    if(x > 0) return 1.0;
    return -1.0;
}

extern "C" {
    
//    SEXP test_icph_llk(SEXP S, SEXP expEta, SEXP lind, SEXP rind);
    SEXP test_nne(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars);
    SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals);
//    SEXP findLR_ind(SEXP lVal, SEXP rVal, SEXP MI_l, SEXP MI_r);
}
#endif
