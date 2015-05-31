//
//  ic_sp_cm.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/25/15.
//
//

#include "ic_sp_ch.h"
//#include "../icenReg_files/basicUtilities.cpp"


/*      LIKELIHOOD TOOLS        */
void actSet_Abst::update_p_ob(int i){
    double chl = baseCH[ obs_inf[i].l ];
    double chr = baseCH[ obs_inf[i].r +1 ];
    double eta = etas[i];
    obs_inf[i].pob = basHaz2CondS(chl, eta) - basHaz2CondS(chr, eta);
}

double actSet_Abst::sum_llk(){
    int n = obs_inf.size();
    double ans = 0;
    for(int i = 0; i < n; i++){
        update_p_ob(i);
        ans += log(obs_inf[i].pob);
    }
    if(isnan(ans)) ans = R_NegInf;
    return(ans);
}

double actSet_Abst::par_llk(int a_ind){
    int par_n = actIndex[a_ind].dep_obs.size();
    double ans = 0;
    int thisInd;
    for(int i = 0; i < par_n; i++){
        thisInd = actIndex[a_ind].dep_obs[i];
        update_p_ob(thisInd);
        ans += log(obs_inf[thisInd].pob);
    }
    if(isnan(ans)) ans = R_NegInf;
    return(ans);
}

void actSet_Abst::update_etas(){
    etas = covars * reg_par;
    for(int i = 0; i < etas.size(); i++)
        expEtas[i] = exp(etas[i]);
}


/*      ACTIVE SET MANIPULATION TOOLS       */
int actSet_Abst::getNextRawActInd(int act_i){
    int ak = actIndex.size();
    if(act_i >= ak || act_i < 0){Rprintf("error in getNextRawActInd: invalid act_i!\n"); return -1;}
    if(act_i == (ak-1))    return(baseCH.size() - 1);//return(actIndex[act_i].ind + 1);
    return(actIndex[act_i + 1].ind);
}

void actSet_Abst::act_setPar(int act_i, double val){
    int l = actIndex[act_i].ind;
    int r = getNextRawActInd(act_i) - 1;
    actIndex[act_i].par = val;
    for(int i = l; i <= r; i++){
        baseCH[i] = val;
    }
}

void actSet_Abst::act_addPar(vector<double> &delta){
    int p_k = delta.size();
    int a_k = actIndex.size();
    if(p_k != a_k){Rprintf("in act_addPar, delta is not the same length as actIndex!\n");return;}
    for(int i = 0; i < p_k; i++)
        act_addPar(i, delta[i]);
}

void actSet_Abst::addDepNodes(vector<int> &intoVec, int l, int r){
    int tot_size = 0;
    for(int i = l; i <= r; i++){
        tot_size += node_inf[i].l.size();
    }
    tot_size *= 2;
    
    intoVec.clear();
    intoVec.reserve(tot_size);
    for(int i = l; i <= r; i++){
        intoVec.insert(intoVec.end(), node_inf[i].l.begin(), node_inf[i].l.end() );
        intoVec.insert(intoVec.end(), node_inf[i].r.begin(), node_inf[i].r.end() );
    }
    sort(intoVec.begin(), intoVec.end());
    intoVec.erase(unique(intoVec.begin(), intoVec.end()), intoVec.end());
}

void actSet_Abst::addActive(int raw_ind, double par){
    if(raw_ind == 0 || raw_ind == (baseCH.size() - 1)) {return;}
        if(actIndex.size() == 0) {
        Rprintf("warning: trying to add actIndex, but haven't initialized actIndex with a first active point!\n");
        return;
    }
    
    int k_act = actIndex.size();
    actPointInf* nextActPt;
    
    nextActPt = &actIndex[0];
    if(raw_ind < nextActPt->ind){
        
        actPointInf newActPt;
        newActPt.par = nextActPt->par;
        newActPt.ind = raw_ind;
        int lower_ind = nextActPt->ind;
        int upcoming_ind = actIndex[1].ind;
        
        addDepNodes(newActPt.dep_obs, raw_ind, lower_ind -1);
        addDepNodes(nextActPt->dep_obs, lower_ind, upcoming_ind -1);
        actIndex.insert(actIndex.begin(), newActPt);
        act_setPar(0, par);
        return;
    }
    for(int i = 0; i < (k_act); i++){
        nextActPt = &actIndex[i];
        if(raw_ind == nextActPt->ind) {return;}
        if(raw_ind >= nextActPt->ind && raw_ind < getNextRawActInd(i)){
            
            actPointInf newActPt;
            newActPt.ind = raw_ind;
            
            int lower_ind = nextActPt->ind;
            int upper_ind = getNextRawActInd(i);
            addDepNodes(nextActPt->dep_obs, lower_ind, raw_ind - 1);
            addDepNodes(newActPt.dep_obs, raw_ind, upper_ind -1);
            actIndex.insert(actIndex.begin() + i + 1, newActPt);
            act_setPar(i + 1, par);
            act_addPar(i, 0);
            return;
        }
    }
    
    nextActPt = &actIndex[k_act-1];
    
    actPointInf newActPt;
    newActPt.ind = raw_ind;
    
    int lower_ind = nextActPt->ind;
    int upper_ind = getNextRawActInd(actIndex.size()- 1);
    if(lower_ind == upper_ind) return;
    addDepNodes(nextActPt->dep_obs, lower_ind, raw_ind - 1);
    addDepNodes(newActPt.dep_obs, raw_ind, upper_ind);

    actIndex.insert(actIndex.end(), newActPt);
    act_setPar(k_act, par);
    act_addPar(k_act - 1, 0);
}

void actSet_Abst::removeActive(int act_ind){
    vector<actPointInf>::iterator iter;
    int ak = actIndex.size();
    if(act_ind <= 0 || act_ind >= (ak) ){
        Rprintf("invalid act_ind for removeActive\n");
        return;
    }
    
    iter = actIndex.begin() + act_ind;
    actIndex.erase(iter);
    
    int upper_ind = getNextRawActInd(act_ind-1);
    addDepNodes(actIndex[act_ind-1].dep_obs, actIndex[act_ind-1].ind, upper_ind);
    act_addPar(act_ind-1, 0.0);
}

void printIndexInfo(actSet_Abst &actSet, int i){
    actPointInf* a = &actSet.actIndex[i];
    Rprintf("ind = %d, par = %f, dep_obs = ", a->ind, a->par);
    for(int i = 0; i < a->dep_obs.size(); i++)  { Rprintf(" %d ", a->dep_obs[i]);};
    Rprintf("\n");
}


int actSet_Abst::getActInd(int raw_ind){
    for(int i = 0; i < actIndex.size(); i++){ if(raw_ind == actIndex[i].ind) return(i);}
    return(-1);
}

void actSet_Abst::checkIfActShouldDrop(int act_ind){
    if(act_ind == 0) return;
    double thispar = actIndex[act_ind].par;
    double lowerpar = actIndex[act_ind-1].par;
    if( abs(thispar -lowerpar) < pow(10, -12))   removeActive(act_ind);
}

int actSet_Abst::getNextActRawInd(int raw_ind){
    for(int i = 0; i < actIndex.size(); i++){
        if(raw_ind < actIndex[i].ind)   return(actIndex[i].ind);
    }
    return(baseCH.size()-1);
};



/*      INITIALIZATION TOOLS    */
void setupActSet(SEXP Rlind, SEXP Rrind, SEXP RCovars, actSet_Abst* actSet){
    actSet->h = 0.0001;
    int n = LENGTH(Rlind);
    if(n != LENGTH(Rrind)){Rprintf("length of Rlind and Rrind not equal\n"); return;}
    actSet->base_p_obs.resize(n);
    actSet->etas.resize(n);
    actSet->expEtas.resize(n);
    
    for(int i = 0; i < n; i++){
        actSet->etas[i] = 0;
        actSet->expEtas[i] = 1;
    }
    
    copyRmatrix_intoEigen(RCovars, actSet->covars);
    int reg_k = actSet->covars.cols();
    if(reg_k == 0) actSet->hasCovars = false; else actSet->hasCovars = true;
    if(reg_k > 0){
        if(n != actSet->covars.rows()) {Rprintf("covar rows not equal to n!\n"); return;}
    }
    actSet->reg_d1.resize(reg_k);
    actSet->reg_d2.resize(reg_k, reg_k);
    actSet->reg_par.resize(reg_k);
    
    int maxInd = 0;
    for(int i = 0; i < n; i++){
        maxInd = max(maxInd, INTEGER(Rrind)[i]);
    }

    

    actSet->baseCH.resize(maxInd + 2);
    for(int i = 0; i < maxInd; i++)
        actSet->baseCH[i] = R_NegInf;
    actSet->baseCH[maxInd+1] = R_PosInf;
    
    actSet->H_d1.resize(maxInd - 1);
    actSet->H_d2.resize(maxInd - 1, maxInd - 1);
    
    vector<int> minActPoints;
    int this_l, this_r;
    actSet->obs_inf.resize(n);
    actSet->node_inf.resize(maxInd + 2);
    
    
    for(int i = 0; i < n; i++){
        this_l = INTEGER(Rlind)[i];
        this_r = INTEGER(Rrind)[i];
        addIfNeeded(minActPoints, this_l, this_r, maxInd);
        actSet->obs_inf[i].l = this_l;
        actSet->obs_inf[i].r = this_r;
        actSet->node_inf[this_l].l.push_back(i);
        actSet->node_inf[this_r + 1].r.push_back(i);
    }
    
    sort(minActPoints.begin(), minActPoints.end());
    
    
    double stepSize = 4.0/minActPoints.size();
    double curVal = -2.0;
    
    actPointInf firstActPt;
    
    firstActPt.ind = minActPoints[0];
    firstActPt.par = curVal;
    actSet->addDepNodes(firstActPt.dep_obs, firstActPt.ind, (actSet->baseCH.size() - 2));
    actSet->actIndex.insert(actSet->actIndex.begin(), firstActPt);
    actSet->act_setPar(0, curVal);
    
    for(int i = 1; i < minActPoints.size(); i++){
        curVal += stepSize;
        actSet->addActive( minActPoints[i], curVal);
    }
    
}



/*      OPTIMIZATION TOOLS      */
void actSet_Abst::numericBaseDervsOne(int raw_ind, vector<double> &dvec){
    dvec.resize(2);
    if(raw_ind <= 0 || raw_ind >= (baseCH.size()- 1)){Rprintf("warning: inappropriate choice of ind for numericBaseDervs ind = %d\n", raw_ind); return;}

    bool needToDropWhenDone = getActInd(raw_ind) == -1;
    addActive(raw_ind);
    
    int act_ind = getActInd(raw_ind);
    double llk_st = par_llk(act_ind);
    act_addPar(act_ind, h);
    double llk_h = par_llk(act_ind);
    act_addPar(act_ind, -2*h);
    double llk_l = par_llk(act_ind);
    act_addPar(act_ind, h);
    
    dvec[0] = (llk_h - llk_l)/(2*h);
    dvec[1] = (llk_h + llk_l - 2 * llk_st) / (h * h);
    if(needToDropWhenDone){
        removeActive(act_ind);
    }
}

void actSet_Abst::numericBaseDervsAllAct(vector<double> &d1, vector<double> &d2){
    int ak = actIndex.size();
    d1.resize(ak);
    d2.resize(ak);
    vector<double> ind_dervs(2);
    for(int i = 0; i < ak; i++){
        numericBaseDervsOne(actIndex[i].ind, ind_dervs);
        d1[i] = ind_dervs[0];
        d2[i] = ind_dervs[1];
    }
}

void actSet_Abst::numericBaseDervsAllRaw(vector<double> &d1, vector<double> &d2){
    int k = baseCH.size() - 1;
    d1.resize(k);
    d2.resize(k);
    vector<double> ind_dervs(2);
    for(int i = 1; i < (k-1); i++){
        numericBaseDervsOne(i, ind_dervs);
        d1[i] = ind_dervs[0];
        d2[i] = ind_dervs[1];
    }
}

void actSet_Abst::analyticBase1stDerv(vector<double> &d1){
    int k = baseCH.size() - 1;
    vector<double> raw_dervs(k);
    raw_dervs[0] = R_NegInf;
    for(int i = 1; i < k; i++)
        raw_dervs[i] = 0;
    int n = obs_inf.size();
    
    for(int i = 0; i < n; i++)
        update_p_ob(i);
    
    double this_pob, this_ch, this_cont, this_eta;
    int l, r;
    for(int i = 0; i < n; i++){
        this_pob = obs_inf[i].pob;
        this_eta = etas[i];
        l = obs_inf[i].l;
        r = obs_inf[i].r;
        if(l > 0){
            this_ch = baseCH[l];
            this_cont = base_d1_contr(this_ch, this_pob, this_eta);
            raw_dervs[l] += this_cont;
        }
        if(r < (k-1)){
            this_ch = baseCH[r+1];
            this_cont = -base_d1_contr(this_ch, this_pob, this_eta);
            raw_dervs[r+1] += this_cont;
        }
    }
    d1.resize(k);
    d1[0] = R_NegInf;
    vector<int> surInds(2);
    
  
    for(int i = 1; i < k; i++)
        d1[i] = 0;
    
    int nextInd;
    for(int i = 1; i < (k - 1); i++){
        d1[i] += raw_dervs[i];
        nextInd = getNextActRawInd(i);
        for(int j = i+1; j < nextInd; j++) {
            d1[i] += raw_dervs[j];
        }
    }
    d1[k-1] = raw_dervs[k-1];
}

void actSet_Abst::uniActiveOptim(int raw_ind){
    addActive(raw_ind);
    int act_ind = getActInd(raw_ind);
    vector<double> dv(2);
    numericBaseDervsOne(raw_ind, dv);
    double upperLim = R_PosInf;
    double lowerLim = R_NegInf;
    double curVal = actIndex[act_ind].par;
    if(act_ind > 0)                     lowerLim = actIndex[act_ind-1].par - curVal;
    if(act_ind < actIndex.size() - 1)   upperLim = actIndex[act_ind+1].par - curVal;
    double prop = 0;
    if(isnan(dv[0] || isnan(dv[1] || dv[0] == R_NegInf || dv[0] == R_PosInf || dv[1] == R_NegInf || dv[1] == R_PosInf))) {Rprintf("error: degenerate derivative estimated! quiting uniActiveOptim\n");}
    if(dv[1] < -.00001){
        prop = -dv[0]/dv[1];
    }
    
    else    prop = dv[0];
    
    prop = max(prop, lowerLim);
    prop = min(prop, upperLim);
    double llk_st = par_llk(act_ind);
    act_addPar(act_ind, prop);
    double llk_nw = par_llk(act_ind);
    
    int tries = 0;
    prop *= -1;
    while(llk_nw < llk_st && tries < 10){
        tries++;
        prop *= 0.5;
        act_addPar(act_ind, prop);
        llk_nw = par_llk(act_ind);
    }
    if(llk_nw < llk_st){
        act_addPar(act_ind, prop);
        llk_nw = par_llk(act_ind);
        if( (llk_nw - llk_st) < -pow(10, -12)){
            Rprintf("warning: likelihood decreased in uniActiveOptim. Difference = %f\n", llk_nw - llk_st);
        }
        prop = 0;
    }
    
    prop = -prop;
    checkIfActShouldDrop(act_ind);
    checkIfActShouldDrop(act_ind+1);
}

void actSet_Abst::vem_step(){
    vector<double> d1;
    analyticBase1stDerv(d1);
    int max_d_ind = getMaxIndex<vector<double> >(d1);
    uniActiveOptim(max_d_ind);
}

void actSet_Abst::icm_step(){
    vector<double> d1;
    vector<double> d2;
    numericBaseDervsAllAct(d1, d2);
    for(int i = 0; i < d1.size(); i ++){
        if(isnan(d2[i]))    {Rprintf("warning: d2 isnan!\n"); return;}
        if(d2[i] >= 0)      {Rprintf("warning: d2 >= 0 in icm step\n"); return;}
    }
    vector<double> x(d1.size());
    if(x.size() != actIndex.size()){Rprintf("warning: x.size()! = actIndex.size()\n"); return;}
    for(int i = 0; i < actIndex.size(); i++) x[i] = actIndex[i].par;
    vector<double> prop(d1.size());
    
    double llk_st = sum_llk();
    pavaForOptim(d1, d2, x, prop);
    
    act_addPar(prop);
    double llk_new = sum_llk();
    mult_vec(-1.0, prop);
    int tries = 0;
    while(llk_st > llk_new && tries < 5){
        tries++;
        mult_vec(0.5, prop);
        act_addPar(prop);
        llk_new = sum_llk();
    }
    if(llk_new < llk_st){
        act_addPar(prop);
        llk_new = sum_llk();
    }
    
    for(int i = 0; i < actIndex.size(); i++)
        checkIfActShouldDrop(i);
}

void actSet_Abst::numericRegDervs(){
    int k = reg_par.size();
    vector<double> lk_l(k);
    vector<double> lk_h(k);
    reg_d1.resize(k);
    reg_d2.resize(k,k);
    
    double lk_0 = sum_llk();
    
    
    for(int i = 0; i < k; i++){
        reg_par[i] += h;
        update_etas();
        lk_h[i] = sum_llk();
        reg_par[i] -= 2 * h;
        update_etas();
        lk_l[i] = sum_llk();
        reg_par[i] += h;
        reg_d1[i] = (lk_h[i] - lk_l[i])/(2 * h);
        reg_d2(i,i) = (lk_h[i] + lk_l[i] - 2*lk_0) / (h*h);
    }
    
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            if(i != j){
                reg_par[i] += h;
                reg_par[j] += h;
                update_etas();
                lk_hh = sum_llk();
                reg_par[i] -= 2 * h;
                reg_par[j] -= 2 * h;
                update_etas();
                lk_ll = sum_llk();
                reg_par[i] += h;
                reg_par[j] += h;
                rho = (lk_hh + lk_ll + 2 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2 * h * h);
                reg_d2(i,j) = rho;
                reg_d2(j,i) = rho;
            }
        }
    }
    update_etas();
}

void actSet_Abst::covar_nr_step(){
    numericRegDervs();
    
    double lk_0 = sum_llk();
    int k = reg_par.size();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolve(reg_d2);
    Eigen::VectorXd evals(1);
    evals[0] = 1;
    
    if(esolve.info() == Eigen::Success)
        evals = esolve.eigenvalues();
    int tries = 0;
    double delta = 1;
    while(max(evals) > -0.000001 && tries < 10){
        tries++;
        for(int i = 0; i < k; i++)
            reg_d2(i,i) -= delta;
        delta *= 2;
        esolve.compute(reg_d2);
        if(esolve.info() == Eigen::Success)
            evals = esolve.eigenvalues();
    }
    
    if(max(evals) > 0){
        return;
    }
    Eigen::VectorXd propVec = -reg_d2.ldlt().solve(reg_d1);
    tries = 0;
    reg_par += propVec;
    propVec *= -1;
    update_etas();
    double lk_new = sum_llk();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        reg_par += propVec;
        update_etas();
        lk_new = sum_llk();
    }
}


/*      CALLING ALGORITHM FROM R     */
SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType){
    actSet_Abst* optObj;
    if(INTEGER(fitType)[0] == 1){
        optObj = new actSet_ph;
    }
    else if(INTEGER(fitType)[0] == 2){
        optObj = new actSet_po;
    }
    else { Rprintf("fit type not supported\n");return(R_NilValue);}

    setupActSet(Rlind, Rrind, Rcovars, optObj);

    double llk_old = R_NegInf;
    double llk_new = optObj->sum_llk();
    int tries = 0;
    
    while(tries < 500 && (llk_new - llk_old) > pow(10, -10)){
        tries++;
        llk_old = llk_new;
        
        optObj->vem_step();
        if(optObj->hasCovars)       optObj->covar_nr_step();
        for(int i = 0; i < 5; i++)  optObj->icm_step();
        llk_new = optObj->sum_llk();
        }

 
    if(llk_new < llk_old){
        Rprintf("warning: likelihood decreased! difference = %f\n", llk_new - llk_old);
    }
    
    vector<double> p_hat;
    cumhaz2p_hat(optObj->baseCH, p_hat);
    
    
    SEXP ans = PROTECT(allocVector(VECSXP, 5));
    SEXP R_pans = PROTECT(allocVector(REALSXP,p_hat.size()));
    SEXP R_coef = PROTECT(allocVector(REALSXP, optObj->reg_par.size()));
    SEXP R_fnl_llk = PROTECT(allocVector(REALSXP, 1));
    SEXP R_its = PROTECT(allocVector(REALSXP, 1));
    SEXP R_score = PROTECT(allocVector(REALSXP, optObj->reg_par.size()));
    for(int i = 0; i < p_hat.size(); i++)
        REAL(R_pans)[i] = p_hat[i];
    for(int i = 0; i < optObj->reg_par.size(); i++){
        REAL(R_coef)[i] = optObj->reg_par[i];
        REAL(R_score)[i] = optObj->reg_d1[i];
    }
    REAL(R_fnl_llk)[0] = llk_new;
    REAL(R_its)[0] = tries;
    
    SET_VECTOR_ELT(ans, 0, R_pans);
    SET_VECTOR_ELT(ans, 1, R_coef);
    SET_VECTOR_ELT(ans, 2, R_fnl_llk);
    SET_VECTOR_ELT(ans, 3, R_its);
    SET_VECTOR_ELT(ans, 4, R_score);
    
    UNPROTECT(6);
    return(ans);

}


void cumhaz2p_hat(Eigen::VectorXd &ch, vector<double> &p){
    int k = ch.size();
    vector<double> S(k);
    p.resize(k-1);
    for(int i = 0; i < k; i++)
        S[i] = exp(-exp(ch[i]));
    
    for(int i = 0; i < (k-1); i++)
        p[i] = S[i+1] - S[i];
}