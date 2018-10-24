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

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() { return 1; }
#endif

// estimate the genotype frequencies
SEXP pest_em_c(SEXP pinit, SEXP epinit, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps,
               SEXP nThreads, SEXP EMpara){

  int nInd_c, nSnps_c, ploid_c, *pref, *palt, nThreads_c, maxThreads;
  double *pepinit, *ppinit, nIter, delta;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  ploid_c = 2; //INTEGER(ploid)[0];
  // Define the pointers to the other input R variables
  pepinit = REAL(epinit);
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  ppinit = REAL(pinit);
  // Extract EM parameters
  nIter = REAL(EMpara)[0];
  delta = REAL(EMpara)[1];

    // set up number of threads for OPENMP
  nThreads_c = asInteger(nThreads);
  maxThreads = omp_get_max_threads();
  if (nThreads_c <= 0) {
    // if nThreads is set to zero then use everything
    nThreads_c = maxThreads;
  }
  else if (nThreads_c > maxThreads) {
    // don't allow more threads than the maximum available
    nThreads_c = maxThreads;
  }
  // Set up the output values
  double llvec[nSnps_c], pvec[nSnps_c], epvec[nSnps_c];

  // Perform the EM algorithm.
  int ind, snp;
  #pragma omp parallel for num_threads(nThreads_c) private(ind)
  for(snp = 0; snp < nSnps_c; snp++){
    // define variables
    double llval, prellval, sum, sum2, sum_temp;
    double pfreq, pfreq_new, ep, sumA, sumB;
    int iter = 0, x, a, b;
    double pep[ploid_c + 1], pepNeg[ploid_c + 1], zi[ploid_c + 1], px[ploid_c + 1];
    // initialize some variables
    llval = 0;
    prellval = 0;
    iter = 0;
    sum = 0;
    // Set some parameter values
    pfreq = ppinit[snp];
    ep    = pepinit[snp];
    // Run the EM algorithm
    while((iter < 2) || ((iter < nIter) & ((llval - prellval) > delta))){
      // update algorithm parameters
      iter = iter + 1;
      prellval = llval;
      llval = 0;
      // Set the geno probablities for the new iteration
      pfreq_new = 0;
      sumA = 0;
      sumB = 0;
      // compute p_ep and 1-p_ep
      for(x = 0; x < ploid_c + 1; x++){
        pep[x] = x*(1-ep)/ploid_c + ((ploid_c-x)*ep)/ploid_c;
        pepNeg[x] = 1 - pep[x];
        px[x] = binomial(x, ploid_c-x) * pow(pfreq, x) * pow(1 - pfreq, ploid_c - x);
      }
      // Run the iteration
      for(ind = 0; ind < nInd_c; ind++){
        sum = 0;
        sum2 = 0;
        a = pref[ind + nInd_c*snp];
        b = palt[ind + nInd_c*snp];
        // M-step
        sum_temp = 0;
        for(x = 0; x < ploid_c + 1; x++){
          zi[x] = (pow(pep[x], a) * pow(pepNeg[x], b) * px[x]);
          sum_temp += zi[x];
        }
        for(x = 0; x < ploid_c + 1; x++){
          zi[x] = zi[x]/sum_temp;
          sum += x * zi[x];
          sum2 += zi[x];
        }
        // E-step
        pfreq_new += sum/(ploid_c*sum2*nInd_c);
        sumA += zi[0]*a + zi[2]*b;
        sumB += zi[0]*b + zi[2]*a;
        llval += log(sum_temp);
      }
      // update the genotype probabilities
      pfreq = pfreq_new;
      ep = sumA/(sumA+sumB);
    }
    // Update the output variables
    llvec[snp] = llval;
    pvec[snp]  = pfreq;
    epvec[snp] = ep;
  }
  // Define the output variables
  double *pll, *ppara_p, *ppara_ep;
  SEXP out;
  PROTECT(out = allocVector(VECSXP, 3));
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, nSnps_c));
  SEXP para_p;
  PROTECT(para_p = allocVector(REALSXP, nSnps_c));
  SEXP para_ep;
  PROTECT(para_ep = allocVector(REALSXP, nSnps_c));
  pll = REAL(ll);
  ppara_p = REAL(para_p);
  ppara_ep = REAL(para_ep);
  // update output variables
  for(snp = 0; snp < nSnps_c; snp++){
    pll[snp] = llvec[snp];
    ppara_p[snp] = pvec[snp];
    ppara_ep[snp] = epvec[snp];
  }
  SET_VECTOR_ELT(out, 0, ll);
  SET_VECTOR_ELT(out, 1, para_p);
  SET_VECTOR_ELT(out, 2, para_ep);
  // return results
  UNPROTECT(4);
  return out;
}




