
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include "functions.h"


// Function for computing binomial coefficients
// taken  from this website "https://rosettacode.org/wiki/Evaluate_binomial_coefficients#C"
static unsigned long gcd_ui(unsigned long x, unsigned long y) {
  unsigned long t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  // y1 <- x0 % y0 ; x1 <- y0
  }
  return x;
}

unsigned long binomial(unsigned long a, unsigned long b) {
  unsigned long n, d, g, r = 1;
  n = a + b;
  if (a == 0) return 1;
  if (a == 1) return n;
  if (a >= n) return (a == n);
  if (a > n/2) a = n-a;
  for (d = 1; d <= a; d++) {
    if (r >= ULONG_MAX/n) {  // Possible overflow
      unsigned long nr, dr;  // reduced numerator / denominator
      g = gcd_ui(n, d);  nr = n/g;  dr = d/g;
      g = gcd_ui(r, dr);  r = r/g;  dr = dr/g;
      if (r >= ULONG_MAX/nr) return 0;  // Unavoidable overflow
      r *= nr;
      r /= dr;
      n--;
    } else {
      r *= n--;
      r /= d;
    }
  }
  return r;
}



