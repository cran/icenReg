//
//  icenReg.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#include "Eigen_local/Dense"
#include <stdio.h>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace std;

#include "icenReg_files/basicUtilities.cpp"
#include "icenReg_files/ic_ph.h"
#include "icenReg_files/ic_ph.cpp"
#include "icenReg_files/ic_po.h"
#include "icenReg_files/ic_po.cpp"
#include "icenReg_files/ic_par.cpp"
#include "icenReg_files/ic_sp_ch.cpp"

SEXP ic_sp(SEXP Rp_mass, SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fType){
    if(INTEGER(fType)[0] == 1){
        SEXP ans = ic_ph(Rp_mass, Rlind, Rrind, Rcovars);
        PROTECT(ans);
        UNPROTECT(1);
        return(ans);
    }
    if(INTEGER(fType)[0] == 2){
        SEXP ans = ic_po(Rp_mass, Rlind, Rrind, Rcovars);
        PROTECT(ans);
        UNPROTECT(1);
        return(ans);
    }


    Rprintf("model type not supported yet!\n");
    return(R_NilValue);
}
