
/*
 * pca-utils-stat.c: source code for basic statistics.
 * Copyright (C) 2015 Bradley Worley <geekysuavo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* include the header of type and function definitions. */
#include "pca-utils.h"

double gammap_series (double s, double x) {
  /* declare required variables. */
  double gln, sum, del, ap;

  /* calculate initial values. */
  ap = s;
  gln = lgamma (s);
  del = sum = 1.0 / s;

  /* loop until convergence. */
  for (;;) {
    ++ap;
    del *= x / ap;
    sum += del;
    if (fabs (del) < fabs (sum) * DBL_EPSILON) {
      return sum * exp (-x + s * log (x) - gln);
    }
  }

  /* satisfy the compiler. */
  return NAN;
}

double gammap_fract (double s, double x) {
  /* declare required variables. */
  double an, b, c, d, del, gln, h;
  signed long i;

  /* calculate initial values. */
  gln = lgamma (s);
  b = x + 1.0 - s;
  c = 1.0 / (DBL_MIN / DBL_EPSILON);
  d = 1.0 / b;
  h = d;

  /* loop until convergence. */
  for (i = 1;; i++) {
    an = -((double) i) * (((double) i) - s);
    b += 2.0;
    d = an * d + b;
    if (fabs (d) < (DBL_MIN / DBL_EPSILON)) d = DBL_MIN / DBL_EPSILON;
    c = b + an / c;
    if (fabs (c) < (DBL_MIN / DBL_EPSILON)) c = DBL_MIN / DBL_EPSILON;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs (del - 1.0) <= DBL_EPSILON) break;
  }

  /* return the calculated value. */
  return exp (-x + s * log (x) - gln) * h;
}

/* gammap: regularized incomplete gamma function P(s,x)
 * @s: shape parameter.
 * @x: x-value to evaluate at.
 */
double gammap (double s, double x) {
  /* ensure the arguments are valid. */
  if (x < 0.0 || s <= 0.0) {
    /* output an error message and return nothing. */
    ferror ("invalid arguments to gammap(%lf,%lf)", s, x);
    return NAN;
  }

  /* decide which method to use in the calculation. */
  if (x == 0.0) {
    /* well... umm... that was easy. */
    return 0.0;
  }

  /* determine which method of calculation to use. */
  if (x < s + 1.0) {
    /* direct series. */
    return gammap_series (s, x);
  }
  else {
    /* continuing fractions. */
    return 1.0 - gammap_fract (s, x);
  }

  /* satisfy the compiler. */
  return NAN;
}

/* betacf: incomplete beta function.
 * @x: value of the beta distribution to calculate a p-value for.
 * @a: first parameter of the incomplete beta function.
 * @b: second parameter of the incomplete beta function.
 */
double betacf (double x, double a, double b) {
  double aa, c, d, del, h, qab, qam, qap;
  long m, m2;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;
  d = 1.0 - qab * x / qap;

  d = 1.0 / d;
  h = d;

  for (m = 1; m < 10000; m++) {
    m2 = 2 * m;
    aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    c = 1.0 + aa / c;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    c = 1.0 + aa / c;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs (del - 1.0) <= DBL_EPSILON) break;
  }

  return h;
}

/* ribeta: regularized incomplete beta function.
 * @x: value of the beta distribution to calculate a p-value for.
 * @a: first parameter of the incomplete beta function.
 * @b: second parameter of the incomplete beta function.
 */
double ribeta (double x, double a, double b) {
  /* declare required variables. */
  double bt;

  /* perform bounds checking on input values. */
  if (a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0)
    return NAN;

  /* these are simple cases. */
  if (x == 0.0 || x == 1.0)
    return x;

  bt = exp (lgamma (a + b) - lgamma (a) - lgamma (b)
     + a * log (x) + b * log (1.0 - x));

  if (x < (a + 1.0) / (a + b + 2.0)) {
    return bt * betacf (x, a, b) / a;
  }
  else {
    return 1.0 - bt * betacf (1.0 - x, b, a) / b;
  }

  /* this should never occur. */
  return NAN;
}

/* fcdf: cumulative distribution function of the F distribution.
 * @x: value of the F distribution to calculate a p-value for.
 * @d1: first number of degrees of freedom.
 * @d2: second number of degrees of freedom.
 */
double fcdf (double x, double d1, double d2) {
  /* return the result of the beta function. */
  return ribeta ((d1 * x) / ((d1 * x) + d2), 0.5 * d1, 0.5 * d2);
}

/* npdf: multivariate normal probability density function.
 * @x: the x vector to calculate for.
 * @mean: the mean vector to calculate for.
 * @cov: the covariance matrix to calculate for.
 */
double npdf (struct vector *x,
             struct vector *mean,
             struct matrix *cov) {
  /* declare required variables. */
  struct matrix *L;
  double k, det, pdf;

  /* ensure the inputs are valid. */
  if (x->n != mean->n || x->n != cov->m || x->n != cov->n) {
    /* output an error message and return nothing. */
    ferror ("invalid arguments to normal pdf");
    return NAN;
  }

  /* allocate the cholesky matrix. */
  L = matrix_new (cov->m, cov->n);

  /* calculate the cholesky decomposition. */
  if (!L || !matrix_chol (cov, L)) {
    /* hmm... covariance matrices should always be PSD. */
    ferror ("invalid covariance matrix");
    return NAN;
  }

  /* calculate the dimensionality. */
  k = (double) x->n;

  /* calculate the determinant. */
  det = matrix_chol_det (L);

  /* free the cholesky decomposition. */
  matrix_free (L);

  /* calculate the pdf value. */
  pdf = exp (-0.5 * diff_matrix_diff_mul (x, mean, cov));
  pdf *= pow (2.0 * M_PI, -0.5 * k) * pow (det, -0.5);

  /* return the pdf value. */
  return pdf;
}

/* ncdf: standard normal cumulative distribution function.
 * @x: the value at which to evaluate the function.
 */
double ncdf (double x) {
  /* the cdf is a transformation of the error function. */
  return 0.5 * (1.0 + erf (x / sqrt (2.0)));
}

/* gammacdf: cumulative distribution function for the gamma distribution.
 * @x: the x-value for evaluation.
 * @a: the first parameter.
 * @b: the second parameter.
 */
double gammacdf (double x, double a, double b) {
  /* return the value. */
  return gammap (a, x / b);
}

/* chi2cdf: chi-square cumulative distribution function.
 * @x: the x-value for evaluation.
 * @k: the number of degrees of freedom.
 */
double chi2cdf (double x, double k) {
  /* return the value. */
  return gammacdf (x, k / 2.0, 2.0);
}

/* mahalanobis_pvalue: calculates a p-value from mahalanobis distance.
 * @DM: the mahalanobis distance.
 * @m: the first group point count.
 * @n: the second group point count.
 * @D: the problem dimensionality.
 */
double mahalanobis_pvalue (double DM,
                           unsigned long m,
                           unsigned long n,
                           unsigned long D) {
  /* declare required variables. */
  double na, nb, p, Tsq, F;

  /* gain floating-point references to the point counts. */
  na = (double) m;
  nb = (double) n;

  /* gain a floating-point reference to the dimensionality. */
  p = (double) D;

  /* calculate the T-squared value. */
  Tsq = ((na * nb) / (na + nb)) * DM * DM;

  /* calculate the F statistic. */
  F = ((na + nb - p - 1.0) / (p * (na + nb - 2.0))) * Tsq;

  /* return the null hypothesis p-value. */
  return 1.0 - fcdf (F, p, na + nb - p - 1.0);
}

/* mixture_pvalue: calculates a p-value from a mixture bhattacharyya distance.
 * @DB: the bhattacharyya distance.
 */
double mixture_pvalue (double DB) {
  /* declare required variables. */
  double BC, k, chisq;

  /* calculate the bhattacharyya coefficient. */
  BC = exp (-DB);

  /* calculate the chi-squared statistic. */
  chisq = 4.0 * (1.0 - BC);

  /* calculate the number of degrees of freedom. */
  k = pow ((double) pca_divisions (), 2.0);

  /* return the null hypothesis p-value. */
  return chi2cdf (chisq, k);
}

/* chisq_pvalue: return the p-value for the chi-square metric.
 * @chisq: the chi-square value.
 */
double chisq_pvalue (double chisq) {
  /* declare required variables. */
  double k;

  /* calculate the number of degrees of freedom. */
  k = pow ((double) pca_divisions (), 2.0);

  /* return the value. */
  return 1.0 - chi2cdf (chisq, k);
}

