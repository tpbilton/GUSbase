#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
void extractAlleleCounts(CharacterMatrix fmt, CharacterMatrix dat, NumericMatrix ref, NumericMatrix alt) {
  int nsnps = dat.nrow();    // rows = SNPs
  int nind = dat.ncol();     // columns = individuals
  for(int i=0; i<nsnps; i++) { // SNPs (dat rows)
    std::string format_str = as<std::string>(fmt(i,0));
    std::vector<std::string> fields;
    std::istringstream f(format_str);
    std::string s;
    while (getline(f, s, ':')) fields.push_back(s);
    int ad_idx = -1, ro_idx = -1, ao_idx = -1, dp4_idx = -1;
    for(size_t k=0; k<fields.size(); k++) {
      if(fields[k] == "AD") ad_idx = k;
      if(fields[k] == "RO") ro_idx = k;
      if(fields[k] == "AO") ao_idx = k;
      if(fields[k] == "DP4") dp4_idx = k;
    }
    // Error if none of the required fields are present
    if (ad_idx == -1 && dp4_idx == -1 && (ro_idx == -1 || ao_idx == -1)) {
      std::string msg = "FORMAT string in SNP row " + std::to_string(i+1) + " must contain AD, DP4, or both RO and AO fields.";
      Rcpp::stop(msg);
    }
    for(int j=0; j<nind; j++) { // individuals (dat columns)
      std::string dat_str = as<std::string>(dat(i,j)); // dat: row=SNP, col=ind
      if(dat_str == "." || dat_str == "./." || dat_str == ".|.") {
        ref(j,i) = NA_REAL;
        alt(j,i) = NA_REAL;
        continue;
      }
      std::vector<std::string> vals;
      std::istringstream d(dat_str);
      while (getline(d, s, ':')) vals.push_back(s);
      if(ad_idx != -1 && vals.size() > ad_idx) {
        std::istringstream ad(vals[ad_idx]);
        std::string ad_ref, ad_alt;
        getline(ad, ad_ref, ',');
        getline(ad, ad_alt, ',');
        ref(j,i) = ad_ref.empty() ? NA_REAL : atof(ad_ref.c_str());
        alt(j,i) = ad_alt.empty() ? NA_REAL : atof(ad_alt.c_str());
      } else if(ro_idx != -1 && ao_idx != -1 && vals.size() > std::max(ro_idx, ao_idx)) {
        ref(j,i) = vals[ro_idx].empty() ? NA_REAL : atof(vals[ro_idx].c_str());
        alt(j,i) = vals[ao_idx].empty() ? NA_REAL : atof(vals[ao_idx].c_str());
      } else if(dp4_idx != -1 && vals.size() > dp4_idx) {
        std::istringstream dp(vals[dp4_idx]);
        double rF=NA_REAL, rR=NA_REAL, aF=NA_REAL, aR=NA_REAL;
        std::string v;
        if(getline(dp, v, ',')) rF = atof(v.c_str());
        if(getline(dp, v, ',')) rR = atof(v.c_str());
        if(getline(dp, v, ',')) aF = atof(v.c_str());
        if(getline(dp, v, ',')) aR = atof(v.c_str());
        if(rF==NA_REAL || rR==NA_REAL) ref(j,i) = NA_REAL; else ref(j,i) = rF + rR;
        if(aF==NA_REAL || aR==NA_REAL) alt(j,i) = NA_REAL; else alt(j,i) = aF + aR;
      } else {
        ref(j,i) = NA_REAL;
        alt(j,i) = NA_REAL;
      }
    }
  }
}
