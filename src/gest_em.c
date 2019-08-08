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
SEXP gest_em_c(SEXP ginit, SEXP ep, SEXP ploid, SEXP ref, SEXP alt, SEXP nInd, SEXP nSnps,
               SEXP nThreads, SEXP EMpara){

  int nInd_c, nSnps_c, ploid_c, *pref, *palt, nThreads_c, maxThreads;
  double *pep_c, *pginit, nIter, delta;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  ploid_c = INTEGER(ploid)[0];
  // Define the pointers to the other input R variables
  pep_c = REAL(ep);
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  pginit = REAL(ginit);
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
  double llvec[nSnps_c], gvec[ploid_c * nSnps_c];


  // Perform the EM algorithm.
  int ind, snp;
#pragma omp parallel for num_threads(nThreads_c) private(ind)
  for(snp = 0; snp < nSnps_c; snp++){
    // define variables
    double llval, prellval, sum;
    double geno_temp[ploid_c + 1], geno[ploid_c + 1], geno_new[ploid_c + 1];
    int iter = 0, x, a, b;
    double pep[ploid_c + 1], pepNeg[ploid_c + 1];
    // initialize some variables
    llval = 0;
    prellval = 0;
    iter = 0;
    sum = 0;
    // Set some parameter values
    for(x = 0; x < ploid_c; x++){
      geno[x] = pginit[x + snp*ploid_c];
      sum += geno[x];
    }
    geno[ploid_c] = 1.0 - sum;
    // compute p_ep and 1-p_ep as these are unchanged for specific value of ep
    for(x = 0; x < ploid_c + 1; x++){
      pep[x] = x*(1-pep_c[snp])/ploid_c + ((ploid_c-x)*pep_c[snp])/ploid_c;
      pepNeg[x] = 1 - pep[x];
    }
    // Run the EM algorithm
    while((iter < 2) || ((iter < nIter) & ((llval - prellval) > delta))){
      // update algorithm parameters
      iter = iter + 1;
      prellval = llval;
      llval = 0;
      // Set the geno probablities for the new iteration
      for(x = 0; x < ploid_c + 1; x++){
        geno_new[x] = 0;
      }
      // Run the iteration
      for(ind = 0; ind < nInd_c; ind++){
        sum = 0;
        a = pref[ind + nInd_c*snp];
        b = palt[ind + nInd_c*snp];
        // M-step
        for(x = 0; x < ploid_c + 1; x++){
          geno_temp[x] = pow(pep[x], a) * pow(pepNeg[x], b) * geno[x];
          sum += geno_temp[x];
        }
        // E-step
        for(x = 0; x < ploid_c + 1; x++){
          geno_new[x] += geno_temp[x]/(sum*nInd_c*1.0);
        }
        llval += log(sum);
      }
      // update the genotype probabilities
      for(x = 0; x < ploid_c + 1; x++){
        geno[x] = geno_new[x];
      }
    }
    // Update the output variables
    llvec[snp] = llval;
    for(x = 0; x < ploid_c; x++){
      gvec[x + snp*ploid_c] = geno[x];
    }
  }
  // Define the output variables
  double *pll, *ppara;
  SEXP out;
  PROTECT(out = allocVector(VECSXP, 2));
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, nSnps_c));
  SEXP para;
  PROTECT(para = allocVector(REALSXP, ploid_c * nSnps_c));
  pll = REAL(ll);
  ppara = REAL(para);
  // update output variables
  int xx;
  for(snp = 0; snp < nSnps_c; snp++){
    pll[snp] = llvec[snp];
    for(xx = 0; xx < ploid_c; xx++)
      ppara[xx + snp*ploid_c] = gvec[xx + snp*ploid_c];
  }
  SET_VECTOR_ELT(out, 0, ll);
  SET_VECTOR_ELT(out, 1, para);
  // return results
  UNPROTECT(3);
  return out;
}
