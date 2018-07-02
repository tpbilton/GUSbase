
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include "functions.h"
//#include <omp.h>


double der_ep(double pep, int a, int b){
  if((a == 0) & (b==0))
    return 0;
  if(a == 0)
    return -b * powl(1-pep, b-1);
  else if(b == 0)
    return a * powl(pep, a-1);
  else
    return powl(pep, a-1) * powl(1-pep, b-1) * (a - (a+b)*pep);
}



SEXP pest_c(SEXP p, SEXP ep, SEXP v, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps){

  int ind, snp, x, nInd_c, nSnps_c, v_c, *pref, *palt, a, b;
  double sum, sum_der1, sum_der2, temp_y, ptemp, temp_ep, ep_c, *pp;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  v_c = INTEGER(v)[0];
  ep_c = REAL(ep)[0];
  // Define the pointers to the other input R variables
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  pp = REAL(p);
  // Define the output variables
  double score_c[nSnps_c+1], ll_c = 0;
  for(snp = 0; snp < nSnps_c+1; snp++){
    score_c[snp] = 0;
  }
  double *pll, *pscore, temp_x[v_c+1];
  SEXP out;
  PROTECT(out = allocVector(VECSXP, 2));
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  SEXP score;
  PROTECT(score = allocVector(REALSXP, nSnps_c+1));
  pll = REAL(ll);
  pscore = REAL(score);

  // compute p_ep and 1-p_ep as these are unchanged for specific value of ep
  double pep[v_c+1], pepNeg[v_c+1];
  for(x = 0; x < v_c+1; x++){
    pep[x] = x*(1-ep_c)/v_c + ((v_c-x)*ep_c)/v_c;
    pepNeg[x] = 1 - pep[x];
  }
  temp_ep = ep_c*(1 - 2*ep_c);

// #pragma omp parallel for
  for(snp = 0; snp < nSnps_c; snp++){
    ptemp = pp[snp];
    for(x = 0; x < v_c+1; x++){
      temp_x[x] = binomial(x, v_c-x) * powl(ptemp, x) * powl(1-ptemp, v_c-x);
    }
    for(ind = 0; ind < nInd_c; ind++){
      sum = 0;
      sum_der1 = 0;
      sum_der2 = 0;
      for(x = 0; x < v_c+1; x++){
        a = pref[ind + nInd_c*snp];
        b = palt[ind + nInd_c*snp];
        temp_y = powl(pep[x], a) * powl(pepNeg[x], b);
        sum += temp_y * temp_x[x];
        sum_der1 += temp_y * temp_x[x] * (x - v_c*ptemp);
        sum_der2 += ((v_c - 2*x)* temp_ep)/v_c * der_ep(pep[x], a, b) * temp_x[x];
      }
      //Rprintf("sum - %lf: at ind %i and snp %i\n", sum, ind, snp);
      // add contribution to likelihood
      ll_c += log(sum);
      // add contributions to score vector
      score_c[snp] += sum_der1/sum;
      score_c[nSnps_c] += sum_der2/sum;
    }
  }
  pll[0] = ll_c;
  for(snp = 0; snp < nSnps_c+1; snp++){
    pscore[snp] = score_c[snp];
  }
  SET_VECTOR_ELT(out, 0, ll);
  SET_VECTOR_ELT(out, 1, score);
  UNPROTECT(3);
  return out;
}
