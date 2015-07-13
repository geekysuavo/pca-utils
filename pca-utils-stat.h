
/*
 * pca-utils-stat.h: header code for basic statistics.
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

#ifndef __PCA_UTILS_STAT__
#define __PCA_UTILS_STAT__

/* begin function definitions below. */

double gammap (double s, double x);

double betacf (double x, double a, double b);

double ribeta (double x, double a, double b);

double fcdf (double x, double d1, double d2);

double npdf (struct vector *x,
             struct vector *mean,
             struct matrix *cov);

double ncdf (double x);

double gammacdf (double x, double a, double b);

double chi2cdf (double x, double k);

double mahalanobis_pvalue (double DM,
                           unsigned long m,
                           unsigned long n,
                           unsigned long D);

double mixture_pvalue (double DB);

double chisq_pvalue (double chisq);

#endif

