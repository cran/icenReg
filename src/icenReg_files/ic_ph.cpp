//
//  intCoxFast.cpp
//  
//
//  Created by Cliff Anderson Bergman on 4/10/15.
//
//

#include "ic_ph.h"

/*double icph_llk(vector<double> &S, vector<double> &expEta, vector<int> &lind, vector<int> &rind){
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
}   */

double icph_llk_sumOnly(ICPH_OptimInfo &optinfo){
    int n = optinfo.P_obs.size();
    double llk = 0;
    for(int i = 0; i < n; i++)
        llk += log(optinfo.P_obs[i]);
    if(isnan(llk))  llk = R_NegInf;
    return(llk);
}

double icph_llk_keepEta(ICPH_OptimInfo &optinfo){
    int k = optinfo.p_mass.size();
    int n = optinfo.P_obs.size();
    optinfo.S[0] = 1.0;
    for(int i = 0; i < k; i++){ optinfo.S[i+1] = optinfo.S[i] - optinfo.p_mass[i];}
    for(int i = 0; i < n; i++){ update_p_ob(i, optinfo);}
    return(icph_llk_sumOnly(optinfo));
}



/*double min_icph_llk(vector<double> &P_obs, vector<int> &g1, vector<int> &g2){
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
}   */

void update_p_ob(int i, ICPH_OptimInfo &optinfo){
    int l,r, k;
    l = optinfo.lind[i];
    r = optinfo.rind[i];
    
    if(l > r)   Rprintf("warning: l < r in update_p_ob!!!\n");
    
    k = optinfo.p_mass.size();
    double thisEE = optinfo.expEta[i];
    double this_p = R_pow(optinfo.S[l], thisEE);
    if( (r+1) < k )
        this_p -= R_pow(optinfo.S[r+1], thisEE);
    optinfo.P_obs[i] = this_p;
    
//    if(this_p < 0) Rprintf("l = %d, r = %d, S[l] = %f, S[r+1] = %f\n", l, r, optinfo.S[l], optinfo.S[r+1]);
}



void easyExchange(int p1, int p2, double delta, ICPH_OptimInfo &optinfo){
    optinfo.p_mass[p1] += delta;
    optinfo.p_mass[p2] -= delta;
}


void NNE_numeric_der(int p1, int p2, vector<double> &dvec, ICPH_OptimInfo &optinfo){
    double h = optinfo.h;
    double p1_s, p2_s;
    p1_s = optinfo.p_mass[p1];
    p2_s = optinfo.p_mass[p2];
    dvec.resize(2);
    
    double lk0 = icph_llk_keepEta(optinfo);
    easyExchange(p1, p2, h, optinfo);
    double lkh = icph_llk_keepEta(optinfo);
    easyExchange(p1, p2, -2*h, optinfo);
    double lkl = icph_llk_keepEta(optinfo);
    easyExchange(p1, p2, h, optinfo);
    
    dvec[0] = (lkh - lkl)/(2*h);
    dvec[1] = (lkh + lkl - 2*lk0) / (h*h);
    
    optinfo.p_mass[p1] = p1_s;
    optinfo.p_mass[p2] = p2_s;
}

void NNE_optim(int p1, int p2, ICPH_OptimInfo &optinfo){
    vector<double> dvec(2);
    NNE_numeric_der(p1, p2, dvec, optinfo);
    double ps1, ps2;
    ps1 = optinfo.p_mass[p1];
    ps2 = optinfo.p_mass[p2];
    
    
    double prop;
    if(dvec[1] < 0) prop = -dvec[0]/dvec[1];
    else prop = signVal(dvec[0]) * 0.01;
    
    double maxProp = ps2;
    double minProp = -ps1;
    
    prop = min(prop, maxProp);
    prop = max(prop, minProp);
    
    double lk0 = icph_llk_keepEta(optinfo);
    easyExchange(p1, p2, prop, optinfo);
    double lknew = icph_llk_keepEta(optinfo);
    
    int tries = 0;
    int maxtries = 10;
    prop *= -1.0;
    
    while(tries < maxtries && lknew < lk0){
        tries++;
        prop *= 0.5;
        easyExchange(p1, p2, prop, optinfo);
        lknew = icph_llk_keepEta(optinfo);
    }
    
    if(lknew < lk0){
        optinfo.p_mass[p1] = ps1;
        optinfo.p_mass[p2] = ps2;
        lknew = icph_llk_keepEta(optinfo);
        if(lknew < lk0) Rprintf("warning: likelihood decreased in NNE optim!\n");
    }
}

void NNE_all(ICPH_OptimInfo &optinfo){
    
    int k = optinfo.p_mass.size();
    vector<double> orgp(k);
    for(int i = 0; i < k; i++)
        orgp[i] = optinfo.p_mass[i];
    double startllk = icph_llk_keepEta(optinfo);

    
    int p1, p2;
    p1 = 0;
    p2 = 1;
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

    double final_llk = icph_llk_keepEta(optinfo);
    if(final_llk < startllk){
        for(int i = 0; i < k; i++) optinfo.p_mass[i] = orgp[i];
        final_llk = icph_llk_keepEta(optinfo);
    }
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
        leftSide = -b * R_pow( (*S)[l], b - 1 )/pob;
        if((*S)[r+1] <= 0)  r0 = true;
        else                r0 = false;
        if(r == false)             rightSide = b * R_pow( (*S)[r+1], b - 1 ) / pob;
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
            if((*gl)[j] == true) leftSide = l_d_cont[ind];
            else leftSide = 0;
            p_ders[i] += leftSide + rightSide;
        }
    }
    int minInd = 0;
    int maxInd = 1;
    double minD = R_PosInf;
    double maxD = R_NegInf;
    for(int i = 0; i < k; i++){
        if(minD > p_ders[i]){
            if(optinfo.p_mass[i] > 0.000005){
                minInd = i;
                minD = p_ders[i];
            }
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

void vem_update(ICPH_OptimInfo &optinfo, bool recPosd){
    double veryBeginllk = icph_llk_keepEta(optinfo);
    for(int i = 0; i < optinfo.P_obs.size(); i++)
        update_p_ob(i, optinfo);
    double secllk = icph_llk_keepEta(optinfo);
    if(abs(veryBeginllk - secllk) > 0.00000001) Rprintf("warning: p_obs were not updated properly entering vem\n");
    
    
    int k = optinfo.p_mass.size();
    vector<double> d1v(k);
    vector<double> d2v(k);
    analyticDerv_ICM(d1v, d2v, optinfo);
    vector<int> minmaxInd(2);
    double mnd = R_PosInf;
    double mxd = R_NegInf;
//    double mxd_bck = R_NegInf;
    for(int i = 0; i < k; i++){
        if(mnd > d1v[i]){
            if(optinfo.p_mass[i] > optinfo.h){
                mnd = d1v[i];
                minmaxInd[1] = i;
            }
        }
        if(optinfo.p_mass[i] > optinfo.h || !recPosd){
            if(mxd < d1v[i]){
                mxd = d1v[i];
                minmaxInd[0] = i;
            }
        }
  /*      if(mxd_bck < d1v[i]){
            if(optinfo.p_mass[i] > optinfo.h){
                mxd_bck = d1v[i];
                minmaxInd[2] = i;
            }
        }   */
    }
    
    int p1, p2, thisCase;
    p1 = 0;
    p2 = 0;
    double d1, d2, maxProp, minProp, analtyicD;
    maxProp = R_PosInf;
    minProp = R_NegInf;
    if(minmaxInd[0] > minmaxInd[1]){
        p2 = minmaxInd[0];
        p1 = minmaxInd[1];
        analtyicD = mnd - mxd;
        minProp = -optinfo.p_mass[p1];
        maxProp = optinfo.p_mass[p2];
        thisCase = 1;
    }
    else if(minmaxInd[0] < minmaxInd[1]){
        p1 = minmaxInd[0];
        p2 = minmaxInd[1];
        analtyicD = mxd - mnd;
        minProp = -optinfo.p_mass[p1];
        maxProp = optinfo.p_mass[p2];
        thisCase = 2;
    }
    else{
        return;
    }
    
 //   vector<int> g1;
 //   vector<int> g2;
    
//    getNecInd(p1, p2, optinfo.lind, optinfo.rind, g1, g2);
    double lk0 = icph_llk_keepEta(optinfo);             //icph_llk_sumOnly(optinfo);
    double h = optinfo.h;

    double p1m, p2m;
    p1m = optinfo.p_mass[p1];
    p2m = optinfo.p_mass[p2];
    
    double backh = -2 * h;
    easyExchange(p1, p2, h, optinfo);                   //orderedExchange(p1, p2, h, g1, g2, optinfo);
    double lkh =  icph_llk_keepEta(optinfo);            //icph_llk_sumOnly(optinfo);
    easyExchange(p1, p2, backh, optinfo);                   //orderedExchange(p1, p2, -2*h, g1, g2, optinfo);
    double lkl =  icph_llk_keepEta(optinfo);            //icph_llk_sumOnly(optinfo);
    easyExchange(p1, p2, h, optinfo);                   //orderedExchange(p1, p2, h, g1, g2, optinfo);

/*    double lk0_test = icph_llk_keepEta(optinfo);            //icph_llk_sumOnly(optinfo);
    if(abs(lk0_test - lk0) > 0.00001){
        Rprintf("error: lk0_test - lk0 *10^5 = %f. lk0 - veryBeginllk * 10^5 = %f\n, p1 = %d, p2 = %d, k = %d\n",
                (lk0_test - lk0) * pow(10, 5),
                (lk0 - veryBeginllk) * pow(10.0, 5.0),
                p1, p2, optinfo.p_mass.size());
        Rprintf("diff in pmasses = %f, %f, org pmasses = %f, %f\n",
               (p1m - optinfo.p_mass[p1]) * pow(10,13),
               (p2m - optinfo.p_mass[p2]) * pow(10,13),
               p1m, p2m);
    }   */
    
    d1 = (lkh - lkl)/(h+h);
    d2 = (lkh + lkl - 2 * lk0) / (h*h);
    
    double prop;
    if(d2 < 0)
        prop = -d1/d2;
    else{
        prop = d1;
    }
    if(isnan(prop)){
        if(analtyicD < 0)
            prop = minProp;
        else
            prop = maxProp;
    }

    if(d1 * analtyicD < 0){
        if(d1 > 0)
            prop = maxProp;
        else
            prop = minProp;
    }
    
    prop = min(prop, maxProp);
    prop = max(prop, minProp);
    
    easyExchange(p1, p2, prop, optinfo);                            //orderedExchange(p1, p2, prop, g1, g2, optinfo);

    if(optinfo.p_mass[p1] < -optinfo.h) Rprintf("Warning: p_mass[p1] < 0! Case = %d p_mass[p1] = %f\n\n", thisCase, optinfo.p_mass[p1]);
    if(optinfo.p_mass[p2] < -optinfo.h) Rprintf("Warning: p_mass[p2] < 0! Case = %d p_mass[p2] = %f\n\n", thisCase, optinfo.p_mass[p2]);
    
    
    double newllk =  icph_llk_keepEta(optinfo);     //icph_llk_sumOnly(optinfo);
    int tries = 0;
    prop *= -1.0;
    
    bool keepDividing = false;
    if(isnan(newllk)) keepDividing = true;
    else if(lk0 > newllk) keepDividing = true;
    int maxTries = 3;
    
    
    while(tries < maxTries && keepDividing){
        tries++;
        prop = prop/2;
         easyExchange(p1, p2, prop, optinfo);               //orderedExchange(p1, p2, prop, g1, g2, optinfo);
        newllk = icph_llk_keepEta(optinfo);         //icph_llk_sumOnly(optinfo);
        keepDividing = false;
        if(isnan(newllk)) keepDividing = true;
        else if(lk0 > newllk) keepDividing = true;
    }
    if( newllk < lk0){
        easyExchange(p1, p2, prop, optinfo);           //orderedExchange(p1, p2, prop, g1, g2, optinfo);
        newllk = icph_llk_keepEta(optinfo);                 //icph_llk_sumOnly(optinfo);     //min_icph_llk(optinfo.P_obs, g1, g2);
   /*     if(newllk < (lk0 - 0.00000001))
            Rprintf("warning: even after compltete half stepping, newllk still less than lk0 in vem step!llk diff * 10^5 = %f d1 = %f, d2 = %f \n",
                    (newllk-lk0) * pow(10.0, 5.0), d1, d2); */
    }
    if(isnan(prop)){
        Rprintf("isnan(prop)! lk0 = %f, newllk = %f, \n", lk0, newllk);
    }
    if(newllk  < lk0){
        optinfo.p_mass[p1] = p1m;
        optinfo.p_mass[p2] = p2m;
        newllk = icph_llk_keepEta(optinfo);
        if(newllk != lk0)   Rprintf("resting did NOT work\n");
    }
}





void vec_delta(vector<double> &delta, ICPH_OptimInfo &optinfo){
//    int n = optinfo.P_obs.size();
    int k = delta.size();
//    vector<double>* S;
    vector<double>* p_mass;
//    S = &(optinfo.S);
    p_mass = &(optinfo.p_mass);
//    double tot_p = 0;
    for(int i = 0; i < k; i++)  {(*p_mass)[i] += delta[i];}
//    for(int i = 0; i < k; i++)  {(*S)[i+1] = (*S)[i] - (*p_mass)[i];}
//    for(int i = 0; i < n; i++)  {update_p_ob(i, optinfo);}
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
        leftSide = -b * R_pow( (*S)[l], b - 1 )/pob;
        if((*S)[r+1] <= 0)  r0 = true;
        else                r0 = false;
        if(!r0)             rightSide = b * R_pow( (*S)[r+1], b - 1 ) / pob;
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
    vector<double> orgp(k);
    for(int i = 0; i < k; i++)
        orgp[i] = optinfo.p_mass[i];
    
    analyticDerv_ICM(d1, d2, optinfo);
    
    
    
/*    for(int i = 0; i < k; i++)
        prop[i] = 0;
    vec_delta(prop, optinfo);   */
    
    
    double sumDelta = 0;
    for(int i = 0; i < k; i++){
        if(d2[i] < 0)
            prop[i] = -d1[i]/d2[i];
        else
            prop[i] = signVal(d1[i]) * 0.01;
        if(isnan(prop[i])){
            prop[i] = 0;
        }
        sumDelta += prop[i];
    }
    
    sumDelta = sumDelta/k;
    for(int i = 0; i < k; i++){prop[i] -= sumDelta;}
    double outside = 1;
    for(int i = 0; i < k; i++){
        if(prop[i] < - optinfo.p_mass[i]){
            outside += -1 * (prop[i] + optinfo.p_mass[i]);
            prop[i] = -optinfo.p_mass[i];
        }
        if(prop[i] > 1 - optinfo.p_mass[i]){
            outside -= prop[i] - (1 - optinfo.p_mass[i]);
            prop[i] = 1 - optinfo.p_mass[i];
        }
    }
   
    sumDelta = 0;
    for(int i = 0; i < k; i++) { sumDelta += prop[i]; }
    
    for(int i = 0; i < k; i++){
//        prop[i] = ((optinfo.p_mass[i] + prop[i]) / outside) - optinfo.p_mass[i];
          prop[i] = ((optinfo.p_mass[i] + prop[i]) / (1 + sumDelta)) - optinfo.p_mass[i];
    }
    
    for(int i = 0; i < k; i++){if(isnan(prop[i])) return;}
    
    double lk_old = icph_llk_keepEta(optinfo);          //icph_llk_sumOnly(optinfo);
    vec_delta(prop, optinfo);

    double lk_new = icph_llk_keepEta(optinfo);          //icph_llk_sumOnly(optinfo);
    
    mult_vec(-1.0, prop);
    
    int tries = 0;
    bool lk_bad = false;
    if(isnan(lk_new)) lk_bad = true;
    else if(lk_new < lk_old) lk_bad = true;
    
    
    while(tries < 10 && lk_bad){
        tries++;
        mult_vec(0.5, prop);
        vec_delta(prop, optinfo);
        lk_new = icph_llk_keepEta(optinfo);         //icph_llk_sumOnly(optinfo);
        lk_bad = false;
        if(isnan(lk_new)) lk_bad = true;
        else if(lk_new < lk_old) lk_bad = true;
    }
    
    if(lk_bad){
        vec_delta(prop, optinfo);
        lk_new = icph_llk_keepEta(optinfo); //icph_llk_sumOnly(optinfo);
    }
    
    if(lk_new < lk_old){
        for(int i = 0; i < k; i++) optinfo.p_mass[i] = orgp[i];
        lk_new = icph_llk_keepEta(optinfo);
    }
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
    h = 0.001;
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

/*SEXP test_icph_llk(SEXP cdf, SEXP expEta, SEXP lind, SEXP rind){
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
}   */


void update_etas(ICPH_OptimInfo &optinfo){
    optinfo.eta = optinfo.covars * optinfo.betas;
    int n = optinfo.P_obs.size();
    for(int i = 0; i < n; i++)
        optinfo.expEta[i] = exp(optinfo.eta[i]);
//    for(int i = 0; i < n; i++)
//        update_p_ob(i, optinfo);
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
        ps = R_pow(s, expEta);
        left = ls * ps;
        left2 = ls * ls * ps;
        s = (*S)[r+1];
        ls = log(s);
        ps = R_pow(s, expEta);
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

double uni_cov_step(int i, double d1, ICPH_OptimInfo &optinfo){
    double delta = signVal(d1) * 0.05;
    
    int tries = 0;
    int maxTries = 10;
    
    double lk0 = icph_llk_keepEta(optinfo);
    optinfo.betas[i] += delta;
    update_etas(optinfo);
    
    double newlk = icph_llk_keepEta(optinfo);
    while(tries < maxTries && newlk < lk0){
        tries++;
        optinfo.betas[i] -= delta;
        delta *= 0.25;
        optinfo.betas[i] += delta;
        update_etas(optinfo);
        newlk = icph_llk_keepEta(optinfo);
    }
    
    if(newlk <lk0){
        optinfo.betas[i] -= delta;
        delta = -signVal(d1) * 0.1;
        tries = 0;
        maxTries = 10;
        
        optinfo.betas[i] += delta;
        newlk = icph_llk_keepEta(optinfo);
        while(tries < maxTries && newlk < lk0){
            tries++;
            optinfo.betas[i] -= delta;
            delta *= 0.25;
            optinfo.betas[i] += delta;
            update_etas(optinfo);
            newlk = icph_llk_keepEta(optinfo);
        }
    }
//    Rprintf("used uni covar, change (*10^5) = %f \n", (newlk - lk0) * pow(10, 5));
    return(newlk - lk0);
}


void updateCovars(ICPH_OptimInfo &optinfo){
    double llk_old = icph_llk_keepEta(optinfo);
    int k = optinfo.betas.size();
/*    for(int i = 0; i < k; i++)
        optinfo.d_cov[i] = covar_d(i, optinfo); */
    
    Eigen::VectorXd d_cov_test(k);
    Eigen::MatrixXd d_mat_test(k, k);
    covar_d_analytic(optinfo.d_cov, optinfo.d2_mat, optinfo);
    
    
    optinfo.propVec = -optinfo.d2_mat.inverse() * optinfo.d_cov;
    for(int i = 0; i < optinfo.propVec.size(); i++){
        if(isnan(optinfo.propVec[i])) { Rprintf("quiting reg step!\n"); return; }
    }
    optinfo.betas += optinfo.propVec;
    
    update_etas(optinfo);
    double llk_new = icph_llk_keepEta(optinfo);
    
    int tries = 0;
    optinfo.propVec *= -1.0;
    bool lk_bad = false;
    if(isnan(llk_new)) lk_bad = true;
    else if(llk_new < llk_old) lk_bad = true;
    
    while(lk_bad && tries < 10){
        tries++;
        optinfo.propVec *= 0.5;
        optinfo.betas += optinfo.propVec;
        update_etas(optinfo);
        llk_new = icph_llk_keepEta(optinfo);
        lk_bad = false;
        if(isnan(llk_new)) lk_bad = true;
        else if(llk_new < llk_old) lk_bad = true;
    }

    if(lk_bad){
        optinfo.betas += optinfo.propVec;
        update_etas(optinfo);
        llk_new = icph_llk_keepEta(optinfo);
        double uni_dlta;
        for(int i = 0; i < k; i++){
            uni_dlta = uni_cov_step(i, optinfo.d_cov[i], optinfo);
            if(uni_dlta > 0.00000001) {llk_new += uni_dlta; }//break;}
            else if (uni_dlta < -0.00000001) Rprintf("warning: likelihood decreased in uni_dlta step!\n");
        }
    }
 //   Rprintf("Change in llk in updateCovars = %f\n", llk_new - llk_old);
}


SEXP ic_ph(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars){
    ICPH_OptimInfo optObj (Rp_mass, Rlind, Rrind, Rcovars);
    update_etas(optObj);
    int k = LENGTH(Rp_mass);
    int n = LENGTH(Rlind);
    optObj.S[0] = 1;
    for(int i = 0; i < k; i++) optObj.S[i+1] = optObj.S[i] - optObj.p_mass[i];
    for(int i = 0; i < n; i++){ update_p_ob(i, optObj); }
    double llk_old = R_NegInf;
    double llk_new = icph_llk_keepEta(optObj);
    int tries = 0;
    double conv_crit = R_pow(10, -10);
    
    while(tries < 500 && llk_new - llk_old > conv_crit) {
        
        llk_old = llk_new;
        tries++;
        
        NNE_all(optObj);
      
/*        llk_new = icph_llk_keepEta(optObj);
        Rprintf("(NNE)change in llk = %f\n", llk_new - llk_old);
            */
        
        vem_update(optObj, true);

/*        llk_new = icph_llk_keepEta(optObj);
        Rprintf("(VEM)change in llk = %f\n", llk_new - llk_old);    */

        
        ICM_step(optObj);

/*        llk_new = icph_llk_keepEta(optObj);
        Rprintf("(ICM)change in llk = %f\n", llk_new - llk_old);    */

        updateCovars(optObj);
        
/*        llk_new = icph_llk_keepEta(optObj);
        Rprintf("(REG)change in llk = %f\n", llk_new - llk_old);    */

        vem_update(optObj, false);

/*        llk_new = icph_llk_keepEta(optObj);
        Rprintf("(VEM)change in llk = %f\n", llk_new - llk_old);    */

        
        NNE_all(optObj);

 /*       llk_new = icph_llk_keepEta(optObj);
        Rprintf("(NNE)change in llk = %f\n", llk_new - llk_old);    */

        
        updateCovars(optObj);
        llk_new = icph_llk_keepEta(optObj);
        
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

/*
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
}   */
