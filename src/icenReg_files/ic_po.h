//
//  ic_po.h
//  
//
//  Created by Cliff Anderson Bergman on 5/25/15.
//
//

#ifndef ____ic_po__
#define ____ic_po__

class ICPO_OptimInfo{
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
    
    ICPO_OptimInfo(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars);
};

double icph_llk_sumOnly(ICPO_OptimInfo &optinfo);       //done
double icph_llk_keepEta(ICPO_OptimInfo &optinfo);       //done
void update_p_ob(int i, ICPO_OptimInfo &optinfo);       //done
void easyExchange(int p1, int p2, double delta, ICPO_OptimInfo &optinfo); //done
void NNE_numeric_der(int p1, int p2, vector<double> &dvec, ICPO_OptimInfo &optinfo); //done
void NNE_optim(int p1, int p2, ICPO_OptimInfo &optinfo);    //done
void NNE_all(ICPO_OptimInfo &optinfo);                      //done
void vem_update(ICPO_OptimInfo &optinfo, bool recPosd);     //done
void vec_delta(vector<double> &delta, ICPO_OptimInfo &optinfo); //done
void analyticDerv_ICM(vector<double> &d1, vector<double> &d2, ICPO_OptimInfo  &optinfo); //done
void ICM_step(ICPO_OptimInfo &optinfo);    //done
void update_etas(ICPO_OptimInfo &optinfo); //done
void covar_d_analytic(Eigen::VectorXd &d1, Eigen::MatrixXd &dmat, ICPO_OptimInfo &optinfo);
void updateCovars(ICPO_OptimInfo &optinfo);



#endif /* defined(____ic_po__) */
