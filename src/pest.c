/*
##########################################################################
# Genotyping Uncertainty with Sequencing data - Base package (GUSbase)
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
 */

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
      a = pref[ind + nInd_c*snp];
      b = palt[ind + nInd_c*snp];
      for(x = 0; x < v_c+1; x++){
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



// estimate the genotype frequencies
SEXP gest_c(SEXP geno, SEXP ep, SEXP v, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps){

  int ind, snp, x, nInd_c, nSnps_c, v_c, *pref, *palt, a, b;
  double sum, temp_y, temp_v, temp_ep, ep_c, *pgeno;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  v_c = INTEGER(v)[0];
  ep_c = REAL(ep)[0];
  // Define the pointers to the other input R variables
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  pgeno = REAL(geno);
  // Define the output variables
  double score_c[v_c*nSnps_c+1], ll_c = 0;
  for(snp = 0; snp < (v_c*nSnps_c+1); snp++){
    score_c[snp] = 0;
  }
  double *pll, *pscore;
  SEXP out;
  PROTECT(out = allocVector(VECSXP, 2));
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  SEXP score;
  PROTECT(score = allocVector(REALSXP, v_c*nSnps_c+1));
  pll = REAL(ll);
  pscore = REAL(score);

  // compute p_ep and 1-p_ep as these are unchanged for specific value of ep
  double pep[v_c+1], pepNeg[v_c+1];
  for(x = 0; x < v_c+1; x++){
    pep[x] = x*(1-ep_c)/v_c + ((v_c-x)*ep_c)/v_c;
    pepNeg[x] = 1 - pep[x];
  }
  temp_ep = ep_c*(1 - 2*ep_c);


  double sum_der1[v_c], sum_der2 = 0;
  // #pragma omp parallel for
  for(snp = 0; snp < nSnps_c; snp++){
    //gtemp = pgeno[snp];
    //for(x = 0; x < v_c+1; x++){
    //  temp_x[x] = binomial(x, v_c-x) * powl(ptemp, x) * powl(1-ptemp, v_c-x);
    //}
    for(ind = 0; ind < nInd_c; ind++){
      sum = 0;
      //compute ll
      a = pref[ind + nInd_c*snp];
      b = palt[ind + nInd_c*snp];
      temp_v = powl(pep[v_c], a) * powl(pepNeg[v_c], b);
      for(x = 0; x < v_c; x++){
        temp_y = powl(pep[x], a) * powl(pepNeg[x], b);
        sum += temp_y * pgeno[x + v_c*snp];
        sum_der1[x] += (temp_y - temp_v) * pgeno[x + v_c*snp]*(1-pgeno[x + v_c*snp]);  // geno
        sum_der2 += ((v_c - 2*x)* temp_ep)/v_c * der_ep(pep[x], a, b) * pgeno[x + v_c*snp]; // epsilon
      }
      sum += temp_v * pgeno[v_c + v_c*snp]; // contribution of ll at x=2_v
      sum_der2 += (-1.0 * temp_ep)/v_c * der_ep(pep[v_c], a, b) * pgeno[v_c + v_c*snp]; // derivative of epsilon at x=2_v
      // add contribution to likelihood
      ll_c += log(sum);
      // add contributions to score vector
      for(x = 0; x < v_c; x++){
        score_c[x + v_c*snp] += sum_der1[x]/sum;
      }
      //score_c[nSnps_c] += sum_der2;
    }
  }
  score_c[v_c*nSnps_c] = sum_der2/ll_c;
  pll[0] = ll_c;
  for(snp = 0; snp < (v_c*nSnps_c+1); snp++){
    pscore[snp] = score_c[snp];
  }
  SET_VECTOR_ELT(out, 0, ll);
  SET_VECTOR_ELT(out, 1, score);
  UNPROTECT(3);
  return out;
}

