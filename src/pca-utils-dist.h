
/*
 * pca-utils-dist.h: header code for group distance metrics.
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

#ifndef __PCA_UTILS_DIST__
#define __PCA_UTILS_DIST__

/* include the pca utils header. */
#include "pca-utils.h"

/* begin function definitions below. */

double diff_matrix_diff_mul (struct vector *x1,
                             struct vector *x2,
                             struct matrix *M);

double pca_group_distance_euclidean (struct pca_group *groupA,
                                     struct pca_group *groupB);

double pca_group_distance_mahalanobis (struct pca_group *groupA,
                                       struct pca_group *groupB);

double pca_group_distance_bhattacharyya (struct pca_group *groupA,
                                         struct pca_group *groupB);

double pca_group_distance_hellinger (struct pca_group *groupA,
                                     struct pca_group *groupB);

double pca_group_distance_mixture (struct pca_group *groupA,
                                   struct pca_group *groupB);

double pca_group_distance_chisq (struct pca_group *groupA,
                                 struct pca_group *groupB);

int pca_distances (struct pca *P, struct matrix *D,
                   pca_distance_metric metric);

int pca_distances_print (struct pca *P, struct matrix *D);

double pca_distances_min (struct matrix *D,
                          unsigned long *ga,
                          unsigned long *gb);

struct pca_tree *pca_distances_union (struct pca *P, struct matrix *D,
                                      unsigned long gb, unsigned long ga);

#endif

