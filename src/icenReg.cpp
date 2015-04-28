//
//  intCoxFast.cpp
//  
//
//  Created by Cliff Anderson Bergman on 4/10/15.
//
//

#include "intCoxFast.h"

double icph_llk(vector<double> &S, vector<double> &expEta, vector<int> &lind, vector<int> &rind){
    int n = lind.size();
    if(n != rind.size()){
        Rprintf("Warning: length r != length l\n");
        return(R_NegInf);
    }
    int k = S.size();
    int l, r;
    double llk = 0;
    for(int i = 0; i < n; i++){
        l = lind[i];
        r = rind[i];
        
        if(l >= k){
            Rprintf("Error: l >= k!\n");
            return(R_NegInf);
        }
        if(r >= k){
            Rprintf("Error: r >= k!\n");
            return(R_NegInf);
        }
        if(l > r){
            Rprintf("Error: l > r!\n");
            return(R_NegInf);
        }
        
        llk += log( pow(S[l], expEta[i]) - pow(S[r+1], expEta[i]) );
        if(isnan(llk)) llk = R_NegInf;
    }
    return(llk);
}

double icph_llk_sumOnly(ICPH_OptimInfo &optinfo){
    int n = optinfo.P_obs.size();
    double llk = 0;
    for(int i = 0; i < n; i++)
        llk += log(optinfo.P_obs[i]);
    return(llk);
}

void getNecInd(int p_ind1, int p_ind2,
               vector<int> &l_inds, vector<int> &r_inds,
               vector<int> &g1, vector<int> &g2){
    g1.clear();
    g2.clear();
    bool in_g1;
    bool in_g2;
    for(int i = 0; i < l_inds.size(); i++){
        in_g1 = l_inds[i] <= p_ind1 && r_inds[i] >= p_ind1;
        in_g2 = l_inds[i] <= p_ind2 && r_inds[i] >= p_ind2;
        if(in_g1 && !in_g2) g1.push_back(i);
        else if(!in_g1 && in_g2) g2.push_back(i);
        else if(l_inds[i] > p_ind1 && r_inds[i] <= p_ind2)  g1.push_back(i);
    }
    if(g1.size() == 0 || g2.size() == 0) {
        Rprintf("Error: sizes of necessary inds = 0! This shouldn't happen if we only have Turnbull intervals\n");
        Rprintf("p1 = %d, p2 = %d\n",p_ind1, p_ind2);
    }
}

void check_NN_id(int id1, int propID, ICPH_OptimInfo &optInfo){
    if(optInfo.nn_info[id1].nn_id == propID) return;
    getNecInd(id1, propID, optInfo.lind, optInfo.rind, optInfo.nn_info[id1].g1, optInfo.nn_info[id1].g2);
}

double min_icph_llk(vector<double> &P_obs, vector<int> &g1, vector<int> &g2){
    int n1 = g1.size();
    int n2 = g2.size();
    int thisInd;
    double llk = 0;

    for(int i = 0; i < n1; i++){
        thisInd = g1[i];
        llk += log( P_obs[thisInd] );
    }

    for(int i = 0; i < n2; i++){
        thisInd = g2[i];
        llk += log( P_obs[thisInd] );
    }
    if(isnan(llk)) llk = R_NegInf;
    return(llk);
}

void update_p_ob(int i, ICPH_OptimInfo &optinfo){
    int l,r;
    l = optinfo.lind[i];
    r = optinfo.rind[i];
    double thisEE = optinfo.expEta[i];
    double this_p = pow(optinfo.S[l], thisEE);
    if(optinfo.S[r+1] > 0)
        this_p -= pow(optinfo.S[r+1], thisEE);
    optinfo.P_obs[i] = this_p;
    
//    if(this_p < 0) Rprintf("l = %d, r = %d, S[l] = %f, S[r+1] = %f\n", l, r, optinfo.S[l], optinfo.S[r+1]);
}


void orderedExchange(int p1_ind, int p2_ind, double alpha, vector<int> &g1, vector<int> &g2, ICPH_OptimInfo &optInfo){
    optInfo.p_mass[p1_ind] += alpha;
    optInfo.p_mass[p2_ind] -= alpha;
    for(int i = p1_ind; i < p2_ind; i++)    optInfo.S[i+1] -= alpha;
    int thisInd;
    for(int i = 0; i < g1.size(); i++){
        thisInd = g1[i];
        update_p_ob(thisInd, optInfo);
    }
    for(int i = 0; i < g2.size(); i++){
        thisInd = g2[i];
        update_p_ob(thisInd, optInfo);
    }
}

void NNE_exchange(int p1_ind, int p2_ind, double alpha,
                  ICPH_OptimInfo &optinfo){
    orderedExchange(p1_ind, p2_ind, alpha, optinfo.nn_info[p1_ind].g1, optinfo.nn_info[p1_ind].g2, optinfo);
}


void numericDerv_NNE(int p1_ind, int p2_ind, ICPH_OptimInfo &optinfo){
    double h = optinfo.h;
    double lk0, lkl, lkh;
    lk0 = min_icph_llk(optinfo.P_obs, optinfo.nn_info[p1_ind].g1, optinfo.nn_info[p1_ind].g2);
    NNE_exchange(p1_ind, p2_ind, h, optinfo);
    lkh = min_icph_llk(optinfo.P_obs, optinfo.nn_info[p1_ind].g1, optinfo.nn_info[p1_ind].g2);
    NNE_exchange(p1_ind, p2_ind, -2*h, optinfo);
    lkl = min_icph_llk(optinfo.P_obs, optinfo.nn_info[p1_ind].g1, optinfo.nn_info[p1_ind].g2);
    NNE_exchange(p1_ind, p2_ind, h, optinfo);
    
//    Rprintf("lk0 = %f, lkh = %f, lkl = %f, \n", lk0, lkh, lkl);
    
    optinfo.nn_ders[0] = (lkh - lkl)/(2*h);
    optinfo.nn_ders[1] = (lkh + lkl - 2 * lk0)/(h*h);
}

void analyticDerv_VEM(vector<int> &minmax, vector<double> &output_dervs, ICPH_OptimInfo  &optinfo){
    vector<double>* beta = &(optinfo.expEta);
    vector<double>* S = &(optinfo.S);
    
    int n = optinfo.P_obs.size();
    int k = optinfo.p_mass.size();
    vector<double> l_d_cont(n);
    vector<double> r_d_cont(n);
    vector<double> l_d2_cont(n);
    vector<double> r_d2_cont(n);
    vector<double> p_ders(k);
    int l, r;
    double  pob, b, leftSide, rightSide;
    bool r0;
    for(int i = 0; i < n; i++){
        l = optinfo.lind[i];
        r = optinfo.rind[i];
        pob = optinfo.P_obs[i];
        b = (*beta)[i];
        leftSide = -b * pow( (*S)[l], b - 1 )/pob;
        if((*S)[r+1] <= 0)  r0 = true;
        else                r0 = false;
        if(!r0)             rightSide = b * pow( (*S)[r+1], b - 1 ) / pob;
        else                rightSide = 0;
        l_d_cont[i] = leftSide;
        r_d_cont[i] = rightSide;
        
        l_d2_cont[i] = -leftSide *(b-1) / (*S)[l];
        if(!r0)             r_d2_cont[i] = -rightSide * (b-1) / (*S)[r+1];
        else                r_d2_cont[i] = 0;
        }
    
    int ind;
    vector<bool>* gl;
    vector<int>* gr;
    for(int i = 0; i < k; i++){
        p_ders[i] = 0;
        gl = &(optinfo.vem_info[i].g_l);
        gr = &(optinfo.vem_info[i].g_r);
        for(int j = 0; j < (*gr).size(); j++){
            ind = (*gr)[j];
            rightSide = r_d_cont[ind];
            if((*gl)[j]) leftSide = l_d_cont[ind];
            else leftSide = 0;
            p_ders[i] += leftSide + rightSide;
        }
    }
    int minInd = 0;
    int maxInd = 1;
    double minD = R_PosInf;
    double maxD = R_NegInf;
    for(int i = 0; i < k; i++){
        if(minD > p_ders[i] && optinfo.p_mass[i] > 0.000005){
            minInd = i;
            minD = p_ders[i];
        }
        if(maxD < p_ders[i]){
            maxInd = i;
            maxD = p_ders[i];
        }
    }
    
    double mind2 = 0;
    double maxd2 = 0;
    double leftSide2, rightSide2;
    
    gl = &(optinfo.vem_info[minInd].g_l);
    gr = &(optinfo.vem_info[minInd].g_r);
    for(int j = 0; j < (*gr).size(); j++){
        ind = (*gr)[j];
        rightSide = r_d_cont[ind];
        rightSide2 = r_d2_cont[ind];
        if((*gl)[j] == true) {
            leftSide = l_d_cont[ind];
            leftSide2 = l_d2_cont[ind];
        }
        else{
            leftSide = 0;
            leftSide2 = 0;
        }
        mind2 += -(leftSide + rightSide) * (leftSide + rightSide) + leftSide2 + rightSide2;
    }
    
    gl = &(optinfo.vem_info[maxInd].g_l);
    gr = &(optinfo.vem_info[maxInd].g_r);
    for(int j = 0; j < (*gr).size(); j++){
        ind = (*gr)[j];
        rightSide = r_d_cont[ind];
        rightSide2 = r_d2_cont[ind];
        if((*gl)[j] == true) {
            
            leftSide = l_d_cont[ind];
            leftSide2 = l_d2_cont[ind];
        }
        else{
            leftSide = 0;
            leftSide2 = 0;
        }
        maxd2 += -(leftSide + rightSide) * (leftSide + rightSide) + leftSide2 + rightSide2;
    }
    output_dervs.resize(4);
    output_dervs[0] = maxD;
    output_dervs[1] = maxd2;
    output_dervs[2] = minD;
    output_dervs[3] = mind2;
    
    minmax.resize(2);
    minmax[0] = maxInd;
    minmax[1] = minInd;
}

void vem_update(ICPH_OptimInfo &optinfo){
    vector<double> dervInfo(4);
    vector<int> minmaxInd(2);
    
    analyticDerv_VEM(minmaxInd, dervInfo, optinfo);
    
    int p1, p2, thisCase;
    p1 = 0;
    p2 = 0;
    double d1, d2, maxProp, minProp, analtyicD;
    maxProp = R_PosInf;
    minProp = R_NegInf;
    if(minmaxInd[0] > minmaxInd[1]){
        p2 = minmaxInd[0];
        p1 = minmaxInd[1];
        analtyicD = dervInfo[2] - dervInfo[0];
        d2 = dervInfo[3] - dervInfo[1];
        
        minProp = -optinfo.p_mass[p1]/2;
        maxProp = optinfo.p_mass[p2]/2;
//        Rprintf("case 1 p1 mass = %f, p2 mass = %f, maxProp = %f, minProp = %f \n",optinfo.p_mass[p1], optinfo.p_mass[p2], maxProp, minProp);
        thisCase = 1;
    }
    else if(minmaxInd[0] < minmaxInd[1]){
        p1 = minmaxInd[0];
        p2 = minmaxInd[1];
        analtyicD = dervInfo[0] - dervInfo[2];
        d2 = dervInfo[1] - dervInfo[3];
        minProp = -optinfo.p_mass[p1]/2;
        maxProp = optinfo.p_mass[p2]/2;
//        Rprintf("case 2 p1 mass = %f, p2 mass = %f\n",optinfo.p_mass[p1], optinfo.p_mass[p2]);
        thisCase = 2;
    }
    else{
       // Rprintf("Error in VEM step: p1 = p2! p1 = %d, p2 = %d, d_p1 = %f, d_p2 = %f Skipping step\n", minmaxInd[0], minmaxInd[1], dervInfo[0], dervInfo[2]);
        return;
    }
    
    vector<int> g1;
    vector<int> g2;
    
    getNecInd(p1, p2, optinfo.lind, optinfo.rind, g1, g2);
 //   double prop = -d1/d2;
 //   prop = min(prop, maxProp);
//    Rprintf("prop = %f,\n d1 = %f d2 = %f\n", prop, d1, d2);

//    Rprintf("p1 = %d, p2 = %d\n", p1, p2);
   
    
//    orderedExchange(p1, p2, 0, g1, g2, optinfo);

    double lk0 = min_icph_llk(optinfo.P_obs, g1, g2);
    
    
    double h = optinfo.h;
    
 //   double maxPos_h1 = optinfo.p_mass[p2];
 //   double maxPos_h2 = optinfo.p_mass[p2];
    
    if(h > maxProp){
        h = max(minProp, -h);
    }
    
    orderedExchange(p1, p2, h, g1, g2, optinfo);
    double lkh = min_icph_llk(optinfo.P_obs, g1, g2);
    orderedExchange(p1, p2, h, g1, g2, optinfo);
    double lkhh = min_icph_llk(optinfo.P_obs, g1, g2);

    orderedExchange(p1, p2, -2*h, g1, g2, optinfo);

   
//    Rprintf("Analytic d1 = %f, d2 = %f\n", d1, d2);
    
    d1 = (lkhh - lk0)/(h+h);
    
    
    d2 = (lkhh + lk0 - 2 * lkh) / (h*h);
    
//    Rprintf("Numeric  d1 = %f, d2 = %f\n", d1, d2);
    
//    Rprintf("d1 = %f, d2 = %f, llk0 = %f, llkh = %f, llkhh = %f, h = %f\n", d1, d2, lk0, lkh, lkhh, h);
//    Rprintf("denom = %f\n",(lkhh + lk0 - 2 * lkh) );
    
    double prop = -d1/d2;
    
    if(isnan(prop)){
        if(analtyicD < 0)
            prop = minProp;
        else
            prop = maxProp;
    }
    
//    Rprintf("unadjusted prop = %f\n", prop);
    
    prop = min(prop, maxProp);
    prop = max(prop, minProp);
    
    orderedExchange(p1, p2, prop, g1, g2, optinfo);

    if(optinfo.p_mass[p1] < -optinfo.h) Rprintf("Warning: p_mass[p1] < 0! Case = %d p_mass[p1] = %f\n\n", thisCase, optinfo.p_mass[p1]);
    if(optinfo.p_mass[p2] < -optinfo.h) Rprintf("Warning: p_mass[p2] < 0! Case = %d p_mass[p2] = %f\n\n", thisCase, optinfo.p_mass[p2]);
    
    
    double newllk = min_icph_llk(optinfo.P_obs, g1, g2);
 //   Rprintf("oldllk = %f, newllk = %f\n",oldllk, newllk);
    int tries = 0;
    prop *= -1;
    
    bool keepDividing = false;
    if(isnan(newllk)) keepDividing = true;
    else if(lk0 > newllk) keepDividing = true;
    
    while(tries < 10 && keepDividing){
        tries++;
        prop = prop/2;
        orderedExchange(p1, p2, prop, g1, g2, optinfo);
        newllk = min_icph_llk(optinfo.P_obs, g1, g2);
        
        keepDividing = false;
        if(isnan(newllk)) keepDividing = true;
        else if(lk0 > newllk) keepDividing = true;
    }
    if( newllk < lk0){
        orderedExchange(p1, p2, prop, g1, g2, optinfo);
        newllk = min_icph_llk(optinfo.P_obs, g1, g2);

    }
//    Rprintf("In VEM step, change in llk = %f, number of tries = %d, p1 = %d, p2 = %d, prop = %f\n\n", newllk - lk0, tries, p1, p2, prop);
    if(isnan(prop)){
        Rprintf("isnan(prop)! lk0 = %f, newllk = %f, \n", lk0, newllk);
    }
}



void analyticDerv_NNE(int p1_ind, int p2_ind, ICPH_OptimInfo &optinfo){
    vector<double>* beta = &(optinfo.expEta);
    vector<double>* S = &(optinfo.S);
    vector<int>* g1 = &(optinfo.nn_info[p1_ind].g1);
    vector<int>* g2 = &(optinfo.nn_info[p1_ind].g2);
    optinfo.nn_ders[0] = 0;
    optinfo.nn_ders[1] = 0;
    
    int ind, l ,r;
    double dc, dc2, pob, b;
    for(int i = 0; i < (*g1).size(); i++){
        ind = (*g1)[i];
        l = optinfo.lind[ind];
        r = optinfo.rind[ind] + 1;
        pob = optinfo.P_obs[ind];
        b = (*beta)[ind];
        dc = b * ( pow( (*S)[l], b - 1 ) + pow( (*S)[r], b - 1 ) );
        dc2 = b * (b-1) * ( pow( (*S)[l], b - 2 ) - pow( (*S)[r], b - 2 ) );
        optinfo.nn_ders[0] += dc/pob;
        optinfo.nn_ders[1] += -dc * dc /(pob*pob) - dc2 / pob;
        
    }

    for(int i = 0; i < (*g2).size(); i++){
        ind = (*g2)[i];
        l = optinfo.lind[ind];
        r = optinfo.rind[ind] + 1;
        pob = optinfo.P_obs[ind];
        b = (*beta)[ind];
        dc = b * ( pow( (*S)[l], b - 1 ) + pow( (*S)[r], b - 1 ) );
        dc2 = b * (b-1) * ( pow( (*S)[l], b - 2 ) - pow( (*S)[r], b - 2 ) );
        
        optinfo.nn_ders[0] -= dc/pob;
        optinfo.nn_ders[1] -= dc * dc /(pob*pob) - dc2/pob;
    }
}

void NNE_optim(int p1, int p2, ICPH_OptimInfo &optinfo){
    check_NN_id(p1, p2, optinfo);

    numericDerv_NNE(p1, p2, optinfo);
//    analyticDerv_NNE(p1, p2, optinfo);
    double d1 = optinfo.nn_ders[0];
    double d2 = optinfo.nn_ders[1];
    if(d2 >= 0){
//        Rprintf("Warning: d2 = %f in NNE_optim! p1 = %f, p2 = %f, Quiting NNE_optim\n",
  //              d2, optinfo.p_mass[p1], optinfo.p_mass[p2]);
        return;
    }
    
    double delta = -d1/d2;
    if(isnan(delta)){
        if(optinfo.p_mass[p1] < optinfo.p_mass[p2])
            delta = optinfo.p_mass[p2];
        else
            delta = -optinfo.p_mass[p1];
    }
    double maxVal = optinfo.p_mass[p2];
    double minVal = -optinfo.p_mass[p1];
    
    if(delta > maxVal) delta = maxVal;
    if(delta < minVal) delta = minVal;


    
    double lk_old = min_icph_llk(optinfo.P_obs, optinfo.nn_info[p1].g1, optinfo.nn_info[p1].g2);
    
    NNE_exchange(p1, p2, delta, optinfo);
    
    double lk_new = min_icph_llk(optinfo.P_obs, optinfo.nn_info[p1].g1, optinfo.nn_info[p1].g2);
    int tries = 0;
    delta = -delta;
    while(tries < optinfo.hSteps && lk_new < lk_old){
        tries++;
        delta = delta/2;
        NNE_exchange(p1, p2, delta, optinfo);
        lk_new = min_icph_llk(optinfo.P_obs, optinfo.nn_info[p1].g1, optinfo.nn_info[p1].g2);
    }
//    Rprintf("lk_new - lk_old = %f, tries = %d\n", lk_new - lk_old, tries);
    if(lk_new < lk_old) NNE_exchange(p1, p2, delta, optinfo);
}

void NNE_all(ICPH_OptimInfo &optinfo){
    int p1, p2, k;
    p1 = 0;
    p2 = 1;
    k = optinfo.p_mass.size();
    while(p2 < k){
        if(optinfo.p_mass[p2] == 0) p2++;
        else{
            NNE_optim(p1, p2, optinfo);
            p1 = p2;
            p2++;
            while(optinfo.p_mass[p1] <= 0 && p2 <k){
                p1 = p2;
                p2++;
            }
        }
    }
}


double max(double a, double b){
    if(a > b) return(a);
    return(b);
}

double min(double a, double b){
    if(a < b) return (a);
    return(b);
}

void vec_delta(vector<double> &delta, ICPH_OptimInfo &optinfo){
    int n = optinfo.P_obs.size();
    int k = delta.size();
    vector<double>* S;
    vector<double>* p_mass;
    S = &(optinfo.S);
    p_mass = &(optinfo.p_mass);
    for(int i = 0; i < k; i++)  {(*p_mass)[i] += delta[i];}
    for(int i = 0; i < k; i++)  {(*S)[i+1] = (*S)[i] - (*p_mass)[i];}
    for(int i = 0; i < n; i++)  {update_p_ob(i, optinfo);}
}

void mult_vec(double a, vector<double> &vec){
    for(int i = 0; i < vec.size(); i++)
        vec[i] *= a;
}


void analyticDerv_ICM(vector<double> &d1, vector<double> &d2, ICPH_OptimInfo  &optinfo){
    vector<double>* beta = &(optinfo.expEta);
    vector<double>* S = &(optinfo.S);
    int n = optinfo.P_obs.size();
    int k = optinfo.p_mass.size();
    d1.resize(k);
    d2.resize(k);
    vector<double> l_d_cont(n);
    vector<double> r_d_cont(n);
    vector<double> l_d2_cont(n);
    vector<double> r_d2_cont(n);
    int l, r;
    double  pob, b, leftSide, rightSide, leftSide2, rightSide2;
    bool r0;
    for(int i = 0; i < n; i++){
        l = optinfo.lind[i];
        r = optinfo.rind[i];
        pob = optinfo.P_obs[i];
        b = (*beta)[i];
        leftSide = -b * pow( (*S)[l], b - 1 )/pob;
        if((*S)[r+1] <= 0)  r0 = true;
        else                r0 = false;
        if(!r0)             rightSide = b * pow( (*S)[r+1], b - 1 ) / pob;
        else                rightSide = 0;
        l_d_cont[i] = leftSide;
        r_d_cont[i] = rightSide;
        l_d2_cont[i] = -leftSide *(b-1) / (*S)[l];
        if(!r0)             r_d2_cont[i] = -rightSide * (b-1) / (*S)[r+1];
        else                r_d2_cont[i] = 0;
    }
    int ind;
    vector<bool>* gl;
    vector<int>* gr;
    for(int i = 0; i < k; i++){
        d1[i] = 0;
        gl = &(optinfo.vem_info[i].g_l);
        gr = &(optinfo.vem_info[i].g_r);
        for(int j = 0; j < (*gr).size(); j++){
            ind = (*gr)[j];
            rightSide = r_d_cont[ind];
            rightSide2 = r_d2_cont[ind];
            if((*gl)[j]) {
                leftSide = l_d_cont[ind];
                leftSide2 = l_d2_cont[ind];
            }
            else leftSide = 0;
            d1[i] += leftSide + rightSide;
            d2[i] += -(leftSide + rightSide) * (leftSide + rightSide) + leftSide2 + rightSide2;
        }
    }
}


void ICM_step(ICPH_OptimInfo &optinfo){
    int k = optinfo.p_mass.size();
    vector<double> d1(k);
    vector<double> d2(k);
    vector<double> prop(k);
    analyticDerv_ICM(d1, d2, optinfo);
    
    double sumDelta = 0;
    for(int i = 0; i < k; i++){
//        prop[i] = -d1[i]/d2[i];
        prop[i] = d1[i];
        if(isnan(prop[i])){
            prop[i] = 0;
        }
        sumDelta += prop[i];
    }
    sumDelta = sumDelta/k;
    for(int i = 0; i < k; i++){
        prop[i] -= sumDelta;
    }
    double outside = 1;
    for(int i = 0; i < k; i++){
        if(prop[i] < - optinfo.p_mass[i]){
            outside += -prop[i] - optinfo.p_mass[i];
            prop[i] = -optinfo.p_mass[i];
        }
    }
    for(int i = 0; i < k; i++){
        prop[i] = (optinfo.p_mass[i] + prop[i]) / outside - optinfo.p_mass[i];
    }
    
    double lk_old = icph_llk_sumOnly(optinfo);
    vec_delta(prop, optinfo);

    double lk_new = icph_llk_sumOnly(optinfo);
    
    mult_vec(-1.0, prop);
    int tries = 0;
    bool lk_bad = false;
    if(isnan(lk_new)) lk_bad = true;
    else if(lk_new < lk_old) lk_bad = true;
    
    while(tries < 20 && lk_bad){
        tries++;
        mult_vec(0.5, prop);
        vec_delta(prop, optinfo);
        lk_new = icph_llk(optinfo);
        lk_bad = false;
        if(isnan(lk_new)) lk_bad = true;
        else if(lk_new < lk_old) lk_bad = true;
    }
    
    if(lk_bad){
//        Rprintf("should be back to 0...\n");
        vec_delta(prop, optinfo);
        lk_new = icph_llk(optinfo);
    }
//    Rprintf("After ICM step, change in llk = %f, tries = %d\n", lk_new - lk_old, tries);
    
}


ICPH_OptimInfo::ICPH_OptimInfo(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars){
    int n = LENGTH(Rlind);
    if(n != LENGTH(Rrind)){
        Rprintf("Warning: length of l ind and r ind not equal!\n");
        n = 0;
    }
    int k = LENGTH(Rp_mass);
    
    lind.resize(n);
    rind.resize(n);
    eta.resize(n);
    expEta.resize(n);
    P_obs.resize(n);
    for(int i = 0; i < n; i++){
        lind[i] = INTEGER(Rlind)[i];
        rind[i] = INTEGER(Rrind)[i];
        eta[i] = 0;
        expEta[i] = 1;
    }
    p_mass.resize(k);
    S.resize(k+1);
    nn_info.resize(k);
    vem_info.resize(k);
    S[0] = 1;
    for(int j = 0; j < k; j++){
        p_mass[j] = REAL(Rp_mass)[j];
        S[j+1] = S[j] - p_mass[j];
        nn_info[j].nn_id = -1;
        for(int i = 0; i < n; i++){
            if(rind[i] >= j){
                vem_info[j].g_r.push_back(i);
                if(lind[i] > j) vem_info[j].g_l.push_back(true);
                else vem_info[j].g_l.push_back(false);
            }
        }
    }
    if(S[k] != 0){
        double thisSum = S[0] - S[k];
        for(int i = 0; i < k; i++){
            p_mass[i] = p_mass[i] / thisSum;
            S[i+1] = S[i] - p_mass[i];
        }
    }
    
    vdm_ders.resize(2);
    nn_ders.resize(2);
    hSteps = 10;
    h = 0.00001;
    SEXP Rdims;
    Rdims = getAttrib(Rcovars, R_DimSymbol);
    PROTECT(Rdims);
    UNPROTECT(1);
    int numRows = INTEGER(Rdims)[0];
    int numCols = INTEGER(Rdims)[1];
    betas.resize(numCols);
    d_cov.resize(numCols);
    d2_mat.resize(numCols, numCols);
    propVec.resize(numCols);
    
    for(int i = 0; i < betas.size(); i++) betas[i] = 0;
    covars.resize(numRows, numCols);
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++)
            covars(i, j) = REAL(Rcovars)[i + j*numRows];
    }
}

SEXP test_icph_llk(SEXP cdf, SEXP expEta, SEXP lind, SEXP rind){
    int k = LENGTH(cdf);
    vector<double> cCdf(k);
    for (int i = 0; i < k; i++) cCdf[i] = REAL(cdf)[i];
    k = LENGTH(expEta);
    vector<double> cExpEta(k);
    for(int i = 0; i < k; i++) cExpEta[i] = REAL(expEta)[i];
    k = LENGTH(lind);
    vector<int> cLind(k);
    vector<int> cRind(k);
    for(int i = 0; i < k; i++){cLind[i] = INTEGER(lind)[i]; cRind[i] = INTEGER(rind)[i];}
    double output = icph_llk(cCdf, cExpEta, cLind, cRind);
    SEXP ans = PROTECT(allocVector(REALSXP, 1) );
    REAL(ans)[0] = output;
    UNPROTECT(1);
    return(ans);
}


void update_etas(ICPH_OptimInfo &optinfo){
    optinfo.eta = optinfo.covars * optinfo.betas;
    int k = optinfo.eta.size();
    for(int i = 0; i < k; i++)
        optinfo.expEta[i] = exp(optinfo.eta[i]);
    int n = optinfo.P_obs.size();
    for(int i = 0; i < n; i++)
        update_p_ob(i, optinfo);
}

void covar_d_analytic(Eigen::VectorXd &d1, Eigen::MatrixXd &dmat, ICPH_OptimInfo &optinfo){
    int n = optinfo.P_obs.size();
    int k = optinfo.betas.size();
    d1.resize(k);
    dmat.resize(k,k);
    vector<double> dl_dmu(n);
    vector<double> d2l_dmu2(n);
    vector<double>* S = &(optinfo.S);
    vector<int>* rind = &(optinfo.rind);
    vector<int>* lind = &(optinfo.lind);
    double left, right, left2, right2, s, expEta, ls, ps, pob;
    int l, r;
    for(int i = 0; i < n; i++){
        l = (*lind)[i];
        r = (*rind)[i];
        expEta = optinfo.expEta[i];
        pob = optinfo.P_obs[i];
        s = (*S)[l];
        ls = log(s);
        ps = pow(s, expEta);
        left = ls * ps;
        left2 = ls * ls * ps;
        s = (*S)[r+1];
        ls = log(s);
        ps = pow(s, expEta);
        if(s <= 0) { right = 0; right2 = 0;}
        else       { right = ls * ps; right2 = ls * ls * ps; }
        dl_dmu[i] = (left - right)/pob;
        d2l_dmu2[i] =  (left2 - right2)/pob - (dl_dmu[i] * dl_dmu[i]);
    }
    
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++)
            dmat(i,j) = 0;
    }
    
    for(int j = 0; j < k; j++){
        d1[j] = 0;
        for(int i = 0; i < n; i++){
            d1[j] += optinfo.covars(i, j) * optinfo.expEta[i] * dl_dmu[i];
            for(int m = 0; m < k; m++){
                dmat(j,m) += optinfo.covars(i, m) * optinfo.covars(i,j) * optinfo.expEta[i] * (d2l_dmu2[i] * optinfo.expEta[i]  + dl_dmu[i] );
            }
        }
    }
}



double covar_d(int i, ICPH_OptimInfo &optinfo){
    double lk0 = icph_llk_sumOnly(optinfo);
    double h = optinfo.h;
    optinfo.betas[i] += h;
    update_etas(optinfo);
    double lkNew = icph_llk_sumOnly(optinfo);
    optinfo.betas[i] -= h;
    update_etas(optinfo);
    return( (lkNew - lk0)/h);
}

double covar_dd(int i, int j, ICPH_OptimInfo &optinfo){
    double h = optinfo.h;
    
    optinfo.betas[i] += h;
    optinfo.betas[j] += h;
    update_etas(optinfo);
    double llk_hh = icph_llk_sumOnly(optinfo);

    optinfo.betas[i] -= 2*h;
    update_etas(optinfo);
    double llk_lh = icph_llk_sumOnly(optinfo);

    optinfo.betas[j] -= 2*h;
    update_etas(optinfo);
    double llk_ll = icph_llk_sumOnly(optinfo);

    optinfo.betas[i] += 2*h;
    update_etas(optinfo);
    double llk_hl = icph_llk_sumOnly(optinfo);
    
    optinfo.betas[i] -= h;
    optinfo.betas[j] += h;
    update_etas(optinfo);

    double output = (llk_hh + llk_ll - llk_hl - llk_lh)/ (4 * h * h);
    return(output);
}

void updateCovars(ICPH_OptimInfo &optinfo){
    double llk_old = icph_llk_sumOnly(optinfo);
    int k = optinfo.betas.size();
/*    for(int i = 0; i < k; i++)
        optinfo.d_cov[i] = covar_d(i, optinfo); */
    
    Eigen::VectorXd d_cov_test(k);
    Eigen::MatrixXd d_mat_test(k, k);
    covar_d_analytic(optinfo.d_cov, optinfo.d2_mat, optinfo);
    
/*    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++)
            optinfo.d2_mat(i,j) = covar_dd(i, j, optinfo);
    }

    for(int i = 0; i < k; i++){
        Rprintf(" %f ,   ", d_cov_test[i] - optinfo.d_cov[i]);
        for(int j = 0; j < k; j++)
            Rprintf("%f ", optinfo.d2_mat(i,j) - d_mat_test(i,j));
        Rprintf("\n");
    }   */
    
    
    
    optinfo.propVec = -optinfo.d2_mat.inverse() * optinfo.d_cov;
    optinfo.betas += optinfo.propVec;
    
    update_etas(optinfo);
    double llk_new = icph_llk_sumOnly(optinfo);
    
    int tries = 0;
    optinfo.propVec *= -1;
    bool lk_bad = false;
    if(isnan(llk_new)) lk_bad = true;
    else if(llk_new < llk_old) lk_bad = true;
    
    while(lk_bad && tries < 10){
        tries++;
        optinfo.propVec *= 0.5;
        optinfo.betas += optinfo.propVec;
        update_etas(optinfo);
        llk_new = icph_llk_sumOnly(optinfo);
        lk_bad = false;
        if(isnan(llk_new)) lk_bad = true;
        else if(llk_new < llk_old) lk_bad = true;
    }

    if(lk_bad){
        optinfo.betas += optinfo.propVec;
        update_etas(optinfo);
        llk_new = icph_llk_sumOnly(optinfo);
    }
//    Rprintf("in update Covars, change in llk = %f, tries %d\n", llk_new - llk_old, tries);
}


SEXP test_nne(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars){
    ICPH_OptimInfo optObj (Rp_mass, Rlind, Rrind, Rcovars);
    update_etas(optObj);
    int k = LENGTH(Rp_mass);
    int n = LENGTH(Rlind);
    optObj.S[0] = 1;
    for(int i = 0; i < k; i++) optObj.S[i+1] = optObj.S[i] - optObj.p_mass[i];
    for(int i = 0; i < n; i++){ update_p_ob(i, optObj); }
    double llk_old = R_NegInf;
    double llk_new = icph_llk_sumOnly(optObj);
    int tries = 0;
    while(tries < 500 && llk_new - llk_old > 0.0000001) {
        tries++;
        llk_old = llk_new;
        for(int i = 0; i < 2; i++)
            updateCovars(optObj);
        ICM_step(optObj);
        NNE_all(optObj);
        for(int j = 0; j < 2; j++)vem_update(optObj);
        llk_new = icph_llk_sumOnly(optObj);
    }
    
    SEXP ans = PROTECT(allocVector(VECSXP, 5));
    SEXP R_pans = PROTECT(allocVector(REALSXP,k));
    SEXP R_coef = PROTECT(allocVector(REALSXP, optObj.betas.size()));
    SEXP R_fnl_llk = PROTECT(allocVector(REALSXP, 1));
    SEXP R_its = PROTECT(allocVector(REALSXP, 1));
    SEXP R_score = PROTECT(allocVector(REALSXP, optObj.betas.size()));
    for(int i = 0; i < k; i++)
        REAL(R_pans)[i] = optObj.p_mass[i];
    for(int i = 0; i < optObj.betas.size(); i++){
        REAL(R_coef)[i] = optObj.betas[i];
        REAL(R_score)[i] = optObj.d_cov[i];
    }
    REAL(R_fnl_llk)[0] = llk_new;
    REAL(R_its)[0] = tries;
    
    SET_VECTOR_ELT(ans, 0, R_pans);
    SET_VECTOR_ELT(ans, 1, R_coef);
    SET_VECTOR_ELT(ans, 2, R_fnl_llk);
    SET_VECTOR_ELT(ans, 3, R_its);
    SET_VECTOR_ELT(ans, 4, R_score);
//    Rprintf("final llk = %f\n", llk_new);
    
    UNPROTECT(6);
    return(ans);
}


SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals){
    //NOTE: R_AllVals MUST be sorted!!
    int k = LENGTH(R_AllVals);
    vector<double> mi_l;
    vector<double> mi_r;
    
    bool foundLeft = false;
    double last_left = R_NegInf;
    for(int i = 0; i < k; i++){
        if(!foundLeft)                      foundLeft = LOGICAL(isL)[i] == TRUE;
        if(LOGICAL(isL)[i] == TRUE)         last_left = REAL(R_AllVals)[i];
        if(foundLeft){
            if(LOGICAL(isR)[i] == TRUE){
                mi_l.push_back(last_left);
                mi_r.push_back(REAL(R_AllVals)[i]);
                foundLeft = false;
            }
        }
    }
    int tbulls = mi_l.size();
    
    int n = LENGTH(lVals);
    SEXP l_ind = PROTECT(allocVector(INTSXP, n));
    SEXP r_ind = PROTECT(allocVector(INTSXP, n));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < tbulls; j++){
            if(mi_l[j] >= REAL(lVals)[i]){
                INTEGER(l_ind)[i] = j;
                break;
            }
        }
        
        for(int j = tbulls-1; j >= 0; j--){
            if(mi_r[j] <= REAL(rVals)[i]){
                INTEGER(r_ind)[i] = j;
                break;
            }
        }
    }
    
    
    SEXP ans = PROTECT(allocVector(VECSXP, 4));
    SEXP Rl_mi = PROTECT(allocVector(REALSXP, tbulls));
    SEXP Rr_mi = PROTECT(allocVector(REALSXP, tbulls));
    for(int i = 0; i < tbulls; i++){
        REAL(Rl_mi)[i] = mi_l[i];
        REAL(Rr_mi)[i] = mi_r[i];
    }
    SET_VECTOR_ELT(ans, 0, l_ind);
    SET_VECTOR_ELT(ans, 1, r_ind);
    SET_VECTOR_ELT(ans, 2, Rl_mi);
    SET_VECTOR_ELT(ans, 3, Rr_mi);
    UNPROTECT(5);
    return(ans);
}

SEXP findLR_ind(SEXP lVal, SEXP rVal, SEXP MI_l, SEXP MI_r){
    int n = LENGTH(lVal);
    SEXP l = PROTECT(allocVector(INTSXP, n));
    SEXP r = PROTECT(allocVector(INTSXP, n));

    int k = LENGTH(MI_l);
    int last;
    for(int i = 0; i < n; i++){
        for(int j = 0; j <k; j++){
            last = j;
            if(REAL(lVal)[i] <= REAL(MI_l)[j]){
                break;
            }
        }
        INTEGER(l)[i] = last;
        
        for(int j = k-1; j >= 0; j--){
            last = j;
            if(REAL(rVal)[i] >= REAL(MI_r)[j]){
                break;
            }
        }
        INTEGER(r)[i] = last;
    }
    
    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, l);
    SET_VECTOR_ELT(ans, 1, r);
    UNPROTECT(3);
    return(ans);
}
