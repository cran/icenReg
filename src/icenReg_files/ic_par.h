//
//  ic_par.h
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#ifndef ____ic_par__
#define ____ic_par__
class linkFun{
//Abstract class with necessary information about link function
public:
    virtual double con_s(double b_s, double nu) =0;
    //Computes conditional survival probability from baseline S and nu
    virtual double con_d(double b_d, double b_s, double nu) =0;
    //Computes condtional density from baseline density, baseline S and nu
};

class propOdd : public linkFun{
// Proportional Odds link function
public:
    double con_s(double b_s, double nu) { return(nu * b_s / (b_s * (nu -1) + 1));}
    double con_d(double b_d, double b_s, double nu){
        double sqrt_denom = b_s * nu - b_s + 1;
        return( b_d * nu / (sqrt_denom * sqrt_denom));
    }
};

class propHaz : public linkFun{
// Proportional Hazards link function
public:
    double con_s(double b_s, double nu) { return(pow(b_s, nu));}
    double con_d(double b_d, double b_s, double nu){return( b_d * nu * pow(b_s, nu -1));}
};


class baseLineInfo{
//Abstract class with necessary information about baseline distirbution
public:
    virtual double base_d(double x, Eigen::VectorXd &par)=0;
    //computes baseline density
    virtual double base_s(double x, Eigen::VectorXd &par)=0;
    //computes baseline survival
    virtual void update_baseline_vals(Eigen::VectorXd &s_t, Eigen::VectorXd &d_t,
                                      Eigen::VectorXd &s_vals, Eigen::VectorXd &d_vals, Eigen::VectorXd &par)=0;
    //updates the baseline survival probabilities and densities required
};

class parBLInfo : public baseLineInfo{
//Abstract class with necessary info about parametric baseline distribution
    //Really just a placeholder: at some point I would like to add a spline-based baseline distribution
public:
    void update_baseline_vals(Eigen::VectorXd &s_t, Eigen::VectorXd &d_t,
                              Eigen::VectorXd &s_vals, Eigen::VectorXd &d_vals,
                              Eigen::VectorXd &par);
};

class gammaInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(dgamma(x, exp(par[0]), exp(par[1]), 0));}
    double base_s(double x, Eigen::VectorXd &par){return(pgamma(x, exp(par[0]), exp(par[1]), 0, 0));}
};

class weibullInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(dweibull(x, exp(par[0]), exp(par[1]), 0));}
    double base_s(double x, Eigen::VectorXd &par){return(pweibull(x, exp(par[0]), exp(par[1]), 0, 0));}
};

class lnormInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(ic_dlnorm(x, par[0], exp(par[1])));}
    double base_s(double x, Eigen::VectorXd &par){return(ic_plnorm(x, par[0], exp(par[1])));}
};

class expInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(dexp(x, exp(par[0]), 0));}
    double base_s(double x, Eigen::VectorXd &par){return(pexp(x, exp(par[0]), 0, 0));}
};

class loglogisticInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(ic_dloglogistic(x, exp(par[0]), exp(par[1])));}
    double base_s(double x, Eigen::VectorXd &par){return(1 - ic_ploglogistic(x, exp(par[0]), exp(par[1])));}
};



struct dinf{
//structure for containing all the necessary info for uncensored data
    //values d, s and nu are the necessary indicators for looking up the density, survival and nu
public:
    int d;
    int s;
    int nu;
};

struct intInf{
// for general interval censoring
public:
    int l;
    int r;
    int nu;
};

struct rinf{
    //for right censored data
public:
    int l;
    int nu;
};

struct linf{
    //for left censored data
public:
    int r;
    int nu;
};

class IC_parOpt{
public:
    
    baseLineInfo* blInf;
    // Information about the baseline distribution
    
    linkFun* lnkFn;
    // Information about link function
    
    Eigen::VectorXd b_pars;
    //baseline distribution parameters
    Eigen::VectorXd d_b_pars;
    //derivatives of baseline parameters
    Eigen::MatrixXd d2_b_pars;
    //Hessian of baseline parameters
    Eigen::VectorXd betas;
    //regression parameters
    Eigen::VectorXd d_betas;
    // derivatives of betas
    Eigen::MatrixXd d2_betas;
    //Hessian of regression parameters
    Eigen::MatrixXd covars;
    //covariates
    Eigen::VectorXd eta;
    //vector of linear predictions
    Eigen::VectorXd expEta;
    //vector of exponential of linear predictors
    
    Eigen::VectorXd s_t;
    //vector of times associated with survival probabilities to calculate
    Eigen::VectorXd d_t;
    // vector of times associated with densities to caculate
    Eigen::VectorXd s_v;
    //vector of survival probabilities calculated at s_t
    Eigen::VectorXd d_v;
    // vector of densities calculated at d_t
    
    vector<dinf> uc;
    // vector of indicator information about uncensored observations
    // i.e. i's element is a pair of indicators for s_t and d_t values for subject i
    vector<intInf> gic;
    // same as above, but for general interval censored, i.e. l ind and r ind
    vector<linf> lc;
    // same as above, but for left censoring (i.e. indicator for right side of interval)
    vector<rinf> rc;
    // same as above, but for right censoring (i.e. indicator for left side of interval)

    double calcLike_baseReady();
    //calculates the likelihood. Assumes baseline probs are ready
    void calculate_baseline_probs(){
        blInf->update_baseline_vals(s_t,d_t,s_v, d_v,b_pars);};
    double calcLike_all(){
        calculate_baseline_probs();
        return(calcLike_baseReady());
    };
    
    void NR_baseline_pars();
    void NR_reg_pars();

    void update_etas();
    
    void calc_baseline_dervs();
    void numericCovar_dervs();

    
    void fillFullHessianAndScore(SEXP r_mat, SEXP score);
    //for filling out the full hessian and score at the MLE

    double h;


    IC_parOpt(SEXP R_s_t, SEXP R_d_t, SEXP covars,
              SEXP uncenInd, SEXP gicInd, SEXP lInd, SEXP rInd,
              SEXP parType, SEXP linkType);
};

extern "C" {
    SEXP ic_par(SEXP R_s_t,     // this is a vector of times for which the baseline survival will be calculated
                SEXP R_d_t,     // this is a vector of times for which the baseline density will be calculated
                SEXP covars,    // matrix of covariates. Will need to be ordered according to group
                SEXP uncenInd,  // a n_1 x 2 matrix of indicators. First will be index to look up corresponding d
                                // (i.e. for d_t), second is index to look up s_t
                SEXP gicInd,    // a n_2 x 2 matrix of indicators. First will be index to look up left side of interval,
                                // second is to look up right side of interval
                SEXP lInd,      // n_3 vector of indicators of right side of interval
                SEXP rInd,      // n_4 vector of indicators of left side of interval
                SEXP parType,   // integer of paramteric type. 1 = gamma, 2 = weib, 3 = lnorm, 4 = exp, 5 = loglogistic
                SEXP linkType,  // integer of link type. 1 = proportional odds, 2 = proportional hazards
                SEXP outHessian // hessian matrix at MLE. Easier to pass this than to create it in C++
                );
}

#endif /* defined(____ic_par__) */