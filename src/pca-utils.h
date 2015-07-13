
/*
 * pca-utils.h: header code for pca-utils core functions.
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

/* include libc headers, pretty boilerplate. */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

/* include pca-utils headers. */
#include "pca-utils-type.h"
#include "pca-utils-math.h"
#include "pca-utils-rand.h"
#include "pca-utils-stat.h"
#include "pca-utils-list.h"
#include "pca-utils-dist.h"
#include "pca-utils-tree.h"
#include "pca-utils-draw.h"
#include "pca-utils-getopt.h"

#ifndef __PCA_UTILS__
#define __PCA_UTILS__

/* ferror: macro definition of a printf-style function to output error
 * messages to stderr with useful debugging information. */
#define ferror(...) ferror_fn (__FILE__, __LINE__, __VA_ARGS__)

/* define the chi-squared critical value at alpha=0.05, df={2,3}. */
#define CHISQ2 5.99146454710798082032852107658982
#define CHISQ3 7.81472790325117738774451936478727

/* begin function definitions below. */

void ferror_fn (const char *f, const int l, const char *format, ...)
  __attribute__ ((format (printf, 3, 4)));

struct pca_group *pca_group_new (const char *name, unsigned long dims);

struct pca_group *pca_group_union_new (struct pca *P,
                                       unsigned long ga,
                                       unsigned long gb);

struct pca *pca_new (const char *fname);

void pca_group_free (struct pca_group *group);

void pca_free (struct pca *P);

struct pca_group *pca_group_copy (struct pca_group *group);

struct pca_group *pca_group_copy_bootstrapped (struct pca_group *group);

int pca_group_del_point (struct pca_group *group, unsigned long idx);

int pca_group_calculate_mean (struct pca_group *group);

int pca_group_calculate_cov (struct pca_group *group);

int pca_group_calculate_eig (struct pca_group *group);

int pca_group_calculate_map (struct pca_group *group,
                             struct vector *x,
                             struct vector *y,
                             double s);

int pca_calculate_means (struct pca *P);

int pca_calculate_covs (struct pca *P);

int pca_calculate_eigs (struct pca *P);

int pca_calculate_maps (struct pca *P);

int pca_calculate_var (struct pca *P, double *vx, double *vy);

struct pca_group *pca_find_group (struct pca *P, const char *name);

unsigned long pca_group_index (struct pca *P, const char *name);

int pca_add_group (struct pca *P, struct pca_group *group);

struct pca_group *pca_del_group (struct pca *P, unsigned long gdel);

struct pca *pca_copy (struct pca *P);

struct pca *pca_copy_bootstrapped (struct pca *P);

unsigned long pca_dimensions (struct pca *P);

unsigned long pca_divisions (void);

struct matrix *pca_datamatrix (struct pca *P);

#endif

