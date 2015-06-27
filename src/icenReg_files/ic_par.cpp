//
//  ic_par.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#include "ic_par.h"

double IC_parOpt::calcLike_baseReady(){
    double ans = 0;
    int w_ind = -1;
    int thisSize = uc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(lnkFn->con_d(d_v[uc[i].d], s_v[uc[i].s], expEta[uc[i].nu])) * w[w_ind] ;
    }
    thisSize = gic.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(lnkFn->con_s(s_v[gic[i].l], expEta[gic[i].nu])
                   -lnkFn->con_s(s_v[gic[i].r], expEta[gic[i].nu]) ) * w[w_ind];
    }
    thisSize = lc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(1.0 - lnkFn->con_s(s_v[lc[i].r], expEta[lc[i].nu])) * w[w_ind];
    }
    thisSize = rc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(lnkFn->con_s(s_v[rc[i].l], expEta[rc[i].nu])) * w[w_ind];
    }
    
    if(isnan(ans)) ans = R_NegInf;
    return(ans);
}

void parBLInfo::update_baseline_vals(Eigen::VectorXd &s_t, Eigen::VectorXd &d_t,
                                     Eigen::VectorXd &s_vals, Eigen::VectorXd &d_vals,
                                     Eigen::VectorXd &par){
    for(int i = 0; i < s_t.size(); i++){s_vals[i] = base_s(s_t[i], par);}
    for(int i = 0; i < d_t.size(); i++){d_vals[i] = base_d(d_t[i], par);}
}

void IC_parOpt::calc_baseline_dervs(){
    int k = b_pars.size();
    vector<double> lk_l(k);
    vector<double> lk_h(k);
    d_b_pars.resize(k);
    d2_b_pars.resize(k,k);

    double lk_0 = calcLike_all();
    double org_h = h;
    bool bad_derv = true;
    int tries = 0;
    while(tries < 4 && bad_derv){
        bad_derv = false;
        tries++;
        for(int i = 0; i < k; i++){
            b_pars[i] += h;
            lk_h[i] = calcLike_all();
            b_pars[i] -= 2 * h;
            lk_l[i] = calcLike_all();
            b_pars[i] += h;
            d_b_pars[i] = (lk_h[i] - lk_l[i])/(2 * h);
            d2_b_pars(i,i) = (lk_h[i] + lk_l[i] - 2*lk_0) / (h*h);
        
            if(lk_h[i] == R_NegInf || lk_l[i] == R_NegInf){
                bad_derv = true;
                h = h/4;
            }
        }
    }

    if(bad_derv){Rprintf("error: was not able to calculate derivative of baseline parameters!\n");}
        
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            if(i != j){
                b_pars[i] += h;
                b_pars[j] += h;
                lk_hh = calcLike_all();
                b_pars[i] -= 2 * h;
                b_pars[j] -= 2 * h;
                lk_ll = calcLike_all();
                b_pars[i] += h;
                b_pars[j] += h;
                rho = (lk_hh + lk_ll + 2 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2 * h * h);
                d2_b_pars(i,j) = rho;
                d2_b_pars(j,i) = rho;
            }
        }
    }
    
    calculate_baseline_probs();
    h = org_h;
}


void IC_parOpt::NR_baseline_pars(){
    calc_baseline_dervs();

    double lk_0 = calcLike_baseReady();
    int k = b_pars.size();
    Eigen::VectorXd propVec(k);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolve(d2_b_pars);
    Eigen::VectorXd evals(1);
    evals[0] = 1;

    if(esolve.info() == Eigen::Success)
        evals = esolve.eigenvalues();
    
    if(max(evals) < -0.001) { propVec = -d2_b_pars.ldlt().solve(d_b_pars); }
    else{
        for(int i = 0; i < k; i++){
                if(d2_b_pars(i,i) < -0.001) propVec[i] = -d_b_pars[i]/d2_b_pars(i,i);
                else propVec[i] = signVal(d_b_pars[i]) * 0.1;
            
                if(isnan(propVec[i])) propVec[i] = 0;
        }
    }
    
    int tries = 0;
    b_pars += propVec;
    propVec *= -1;
    double lk_new = calcLike_all();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        b_pars += propVec;
        lk_new = calcLike_all();
    }
    if(lk_new < lk_0){
        b_pars += propVec;
        lk_new = calcLike_all();
    }
}


void IC_parOpt::update_etas(){
    eta = covars * betas;
    for(int i = 0; i < eta.size(); i++)
        expEta[i] = exp(eta[i]);
}

void IC_parOpt::numericCovar_dervs(){
    int k = betas.size();
    vector<double> lk_l(k);
    vector<double> lk_h(k);
    d_betas.resize(k);
    d2_betas.resize(k,k);
    
    double lk_0 = calcLike_baseReady();
    
    
    for(int i = 0; i < k; i++){
        betas[i] += h;
        update_etas();
        lk_h[i] = calcLike_baseReady();
        betas[i] -= 2 * h;
        update_etas();
        lk_l[i] = calcLike_baseReady();
        betas[i] += h;
        d_betas[i] = (lk_h[i] - lk_l[i])/(2 * h);
        d2_betas(i,i) = (lk_h[i] + lk_l[i] - 2*lk_0) / (h*h);
    }
    
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            if(i != j){
                betas[i] += h;
                betas[j] += h;
                update_etas();
                lk_hh = calcLike_baseReady();
                betas[i] -= 2 * h;
                betas[j] -= 2 * h;
                update_etas();
                lk_ll = calcLike_baseReady();
                betas[i] += h;
                betas[j] += h;
                rho = (lk_hh + lk_ll + 2 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2 * h * h);
                d2_betas(i,j) = rho;
                d2_betas(j,i) = rho;
            }
        }
    }
    update_etas();
}

void IC_parOpt::fillFullHessianAndScore(SEXP r_mat, SEXP score){
    int k_base = b_pars.size();
    int k_reg = betas.size();
    int k_tot = k_base + k_reg;
    double lk_0 = calcLike_all();
    vector<double> lk_l(k_tot);
    vector<double> lk_h(k_tot);
    for(int i = 0; i < k_base; i++){
        b_pars[i] += h;
        lk_h[i] = calcLike_all();
        b_pars[i] -= 2 * h;
        lk_l[i] = calcLike_all();
        b_pars[i] += h;
        REAL(r_mat)[i + i * k_tot] = (lk_h[i] + lk_l[i] - 2 * lk_0)/(h*h);
        REAL(score)[i] = (lk_h[i] - lk_l[i])/(2*h);
    }
    calculate_baseline_probs();
    int i_tot;
    for(int i = 0; i < k_reg; i++){
        i_tot = i + k_base;
        betas[i] += h;
        update_etas();
        lk_h[i_tot] = calcLike_baseReady();
        betas[i] -= 2*h;
        update_etas();
        lk_l[i_tot] = calcLike_baseReady();
        betas[i] += h;
        REAL(r_mat)[i_tot + i_tot * k_tot] = (lk_l[i_tot] + lk_h[i_tot] - 2*lk_0)/(h * h);
        REAL(score)[i_tot] = (lk_h[i_tot] - lk_l[i_tot])/(2*h);

    }
    update_etas();
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k_tot; i++){
        for(int j = 0; j < k_tot; j++){
            if(i != j){
                if(i <k_base)   {   b_pars[i] += h;}
                else            {   betas[i - k_base] += h;}
                
                if(j < k_base)  {   b_pars[j] += h;}
                else            {   betas[j - k_base] += h;};
                update_etas();
                lk_hh = calcLike_all();

                if(i <k_base)   {   b_pars[i] -= 2*h;}
                else            {   betas[i - k_base] -= 2*h;}
                
                if(j < k_base)  {   b_pars[j] -= 2*h;}
                else            {   betas[j - k_base] -= 2*h;};
                update_etas();
                lk_ll = calcLike_all();

                if(i <k_base)   {   b_pars[i] += h;}
                else            {   betas[i - k_base] += h;}
                
                if(j < k_base)  {   b_pars[j] += h;}
                else            {   betas[j - k_base] += h;};
                
                rho = (lk_hh + lk_ll + 2 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2 * h * h);
                REAL(r_mat)[i + j * k_tot] = rho;
                REAL(r_mat)[j + i * k_tot] = rho;
            }
        }
    }
    update_etas();
}

void IC_parOpt::NR_reg_pars(){
    numericCovar_dervs();
    
    double lk_0 = calcLike_baseReady();
    int k = betas.size();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolve(d2_betas);
    Eigen::VectorXd evals(1);
    evals[0] = 1;
    
    if(esolve.info() == Eigen::Success)
        evals = esolve.eigenvalues();
    int tries = 0;
    double delta = 1;
    while(max(evals) > -0.000001 && tries < 10){
        tries++;
        for(int i = 0; i < k; i++)
            d2_betas(i,i) -= delta;
        delta *= 2;
        esolve.compute(d2_betas);
        if(esolve.info() == Eigen::Success)
            evals = esolve.eigenvalues();
    }
   
    Eigen::VectorXd propVec(k);
    if(max(evals) > 0)  {propVec = -d2_betas.ldlt().solve(d_betas);}
        else{
        for(int i = 0; i < k; i++){
            propVec[i] = 0;
            if(d2_betas(i,i) < 0)   propVec[i] = -d_betas[i] / d2_betas(i,i);
            else propVec[i] = signVal(d_betas[i]) * 0.01;
        
        if(isnan(propVec[i])) propVec[i] = 0;
        }
    }
    tries = 0;
    betas += propVec;
    propVec *= -1;
    update_etas();
    double lk_new = calcLike_all();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        betas += propVec;
        update_etas();
        lk_new = calcLike_all();
    }
    if(lk_new < lk_0){
        betas += propVec;
        update_etas();
        lk_new = calcLike_all();
    }
}




IC_parOpt::IC_parOpt(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
                     SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
                     SEXP R_parType, SEXP R_linkType, SEXP R_w){
    blInf = NULL;
    if(INTEGER(R_parType)[0] == 1) {
        blInf = new gammaInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 2){
        blInf = new weibullInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 3){
        blInf = new lnormInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 4){
        blInf = new expInfo();
        b_pars.resize(1);
        b_pars[0] = 0;
    }
    else if(INTEGER(R_parType)[0] == 5){
        blInf = new loglogisticInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else{Rprintf("warning: parameter type not supported!\n");}
    
    lnkFn = NULL;
    if(INTEGER(R_linkType)[0] == 1) {lnkFn = new propOdd;}
    else if(INTEGER(R_linkType)[0] == 2) {lnkFn = new propHaz;}
    else{Rprintf("warning: link type not supported!\n");}
    
    Rvec2eigen(R_s_t, s_t);
    Rvec2eigen(R_d_t, d_t);
    s_v.resize(s_t.size());
    d_v.resize(d_t.size());
    copyRmatrix_intoEigen(R_covars, covars);
    int k = covars.cols();
    betas.resize(k);
    for(int i = 0; i < k; i++)  betas[i] = 0;
    d_betas.resize(k);
    d2_betas.resize(k, k);
    
    SEXP RuncenDim = getAttrib(R_uncenInd, R_DimSymbol);
    PROTECT(RuncenDim);
    SEXP RgicDim = getAttrib(R_gicInd, R_DimSymbol);
    PROTECT(RgicDim);
    
    int n_1 = INTEGER(RuncenDim)[0];
    int n_2 = INTEGER(RgicDim)[0];
    int n_3 = LENGTH(R_lInd);
    int n_4 = LENGTH(R_rInd);
    
    int tot_n = n_1 + n_2 + n_3 + n_4;
    eta.resize(tot_n);
    expEta.resize(tot_n);
    w.resize(tot_n);
    
    for(int i = 0; i < tot_n; i++){
        eta[i] = 0;
        expEta[i] = 1;
        w[i] = REAL(R_w)[i];
    }
    
    uc.resize(n_1);
    for(int i = 0; i < n_1; i++){
        uc[i].d = INTEGER(R_uncenInd)[i] - 1;
        uc[i].s = INTEGER(R_uncenInd)[i + n_1] - 1;
        uc[i].nu = i;
    }
    
    gic.resize(n_2);
    for(int i = 0; i < n_2; i++){
        gic[i].l = INTEGER(R_gicInd)[i] - 1;
        gic[i].r = INTEGER(R_gicInd)[i + n_2] - 1;
        gic[i].nu = i + n_1;
    }
    
    lc.resize(n_3);
    for(int i = 0; i < n_3; i++){
        lc[i].r = INTEGER(R_lInd)[i] - 1;
        lc[i].nu = i + n_1 + n_2;
    }
    
    rc.resize(n_4);
    for(int i = 0; i < n_4; i++){
        rc[i].l = INTEGER(R_rInd)[i] - 1;
        rc[i].nu = i + n_1 + n_2 + n_3;
    }
    
    h = 0.0001;
    UNPROTECT(2);
}


SEXP ic_par(SEXP R_s_t, SEXP R_d_t, SEXP covars,
            SEXP uncenInd, SEXP gicInd, SEXP lInd, SEXP rInd,
            SEXP parType, SEXP linkType,
            SEXP outHessian, SEXP R_w){
    IC_parOpt optObj = IC_parOpt(R_s_t, R_d_t, covars, uncenInd, gicInd, lInd, rInd, parType, linkType, R_w);
    if(optObj.blInf == NULL) return(R_NilValue);
    if(optObj.lnkFn == NULL) return(R_NilValue);
    double lk_old = R_NegInf;
    int iter = 0;
    int maxIter = 1000;
    double tol = pow(10, -10);
    double lk_new = optObj.calcLike_all();

    if(lk_new == R_NegInf){
        int bk = optObj.b_pars.size();
        int tries = 0;
        double delta = 1;
        while(tries < 10 && lk_new == R_NegInf){
            tries++;
            for(int i = 0; i < bk; i++){
                if(lk_new == R_NegInf){
                    optObj.b_pars[i] = delta;
                    lk_new = optObj.calcLike_all();
                    if(lk_new == R_NegInf)   optObj.b_pars[i] = 0;
                }
            }
            delta *= 5;
        }
    }
    
    if(lk_new == R_NegInf){
        int bk = optObj.b_pars.size();
        int tries = 0;
        double delta = -1;
        while(tries < 10 && lk_new == R_NegInf){
            tries++;
            for(int i = 0; i < bk; i++){
                if(lk_new == R_NegInf){
                    optObj.b_pars[i] = delta;
                    lk_new = optObj.calcLike_all();
                    if(lk_new == R_NegInf)   optObj.b_pars[i] = 0;
                }
            }
            delta *= 5;
        }
    }
    if(lk_new == R_NegInf){
        Rprintf("failed to find adequate starting point! Please contact maintainer of package\n");
        return(R_NilValue);
    }
    
    while(iter < maxIter && lk_new - lk_old > tol){
        lk_old = lk_new;
        iter++;
        optObj.NR_baseline_pars();
        
        optObj.NR_reg_pars();
        lk_new = optObj.calcLike_baseReady();
        
    }
    SEXP score = PROTECT(allocVector(REALSXP, optObj.betas.size() + optObj.b_pars.size() ) );
    optObj.fillFullHessianAndScore(outHessian, score);
    SEXP reg_est = PROTECT(allocVector(REALSXP, optObj.betas.size()));
    SEXP base_est = PROTECT(allocVector(REALSXP, optObj.b_pars.size()));
    SEXP final_llk = PROTECT(allocVector(REALSXP, 1));
    SEXP iters = PROTECT(allocVector(REALSXP, 1));
    for(int i = 0; i < LENGTH(reg_est); i++)    REAL(reg_est)[i] = optObj.betas[i];
    for(int i = 0; i < LENGTH(base_est); i++)   REAL(base_est)[i] = optObj.b_pars[i];
    REAL(final_llk)[0] = optObj.calcLike_baseReady();
    REAL(iters)[0] = iter;
    
    SEXP ans = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, reg_est);
    SET_VECTOR_ELT(ans, 1, base_est);
    SET_VECTOR_ELT(ans, 2, final_llk);
    SET_VECTOR_ELT(ans, 3, iters);
    SET_VECTOR_ELT(ans, 4, outHessian);
    SET_VECTOR_ELT(ans, 5, score);
    UNPROTECT(6);
    
    if(INTEGER(parType)[0] == 1){
        gammaInfo* deleteObj = static_cast<gammaInfo*>(optObj.blInf);
        delete deleteObj;
    }
    if(INTEGER(parType)[0] == 2){
        weibullInfo* deleteObj = static_cast<weibullInfo*>(optObj.blInf);
        delete deleteObj;
    }
    if(INTEGER(parType)[0] == 3){
        lnormInfo* deleteObj = static_cast<lnormInfo*>(optObj.blInf);
        delete deleteObj;
    }
    if(INTEGER(parType)[0] == 4){
        expInfo* deleteObj = static_cast<expInfo*>(optObj.blInf);
        delete deleteObj;
    }
    if(INTEGER(parType)[0] == 5){
        loglogisticInfo* deleteObj = static_cast<loglogisticInfo*>(optObj.blInf);
        delete deleteObj;
    }
    if(INTEGER(linkType)[0] == 1){
        propOdd* deleteObj = static_cast<propOdd*>(optObj.lnkFn);
        delete deleteObj;
    }
    if(INTEGER(linkType)[0] == 2){
        propHaz* deleteObj = static_cast<propHaz*>(optObj.lnkFn);
        delete deleteObj;
    }
    return(ans);
}
