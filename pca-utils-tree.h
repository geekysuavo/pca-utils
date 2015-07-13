
/*
 * pca-utils-tree.h: header code for dendrogram generation.
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

#ifndef __PCA_UTILS_TREE__
#define __PCA_UTILS_TREE__

/* include the pca utils header. */
#include "pca-utils.h"

/* begin function definitions below. */

struct pca_tree *pca_tree_new (struct pca *P, struct matrix *D);

void pca_tree_free (struct pca_tree *T);

int pca_tree_is_leaf (struct pca_tree *T);

unsigned long pca_tree_num_leaves (struct pca_tree *T);

void pca_tree_bootstrap_compare (struct pca_tree *T, struct pca_tree *Tboot);

int pca_tree_bootstrap_end (struct pca_tree *T, unsigned long n);

void pca_tree_recount_points (struct pca_tree *T);

int pca_tree_assign_pvalues (struct pca_tree *T, unsigned long D,
                             pca_distance_metric metric);

void pca_tree_enum_leaves (struct pca_tree *T, unsigned long *idx);

void pca_tree_min_length (struct pca_tree *T, double *l);

void pca_tree_max_length (struct pca_tree *T, double *l);

void pca_tree_max_x (struct pca_tree *T, double *x);

void pca_tree_max_name (struct pca_tree *T, unsigned long *n);

void pca_tree_rescale (struct pca_tree *T,
                       double imin, double imax,
                       double omin, double omax);

void pca_tree_register_y (struct pca_tree *T, double dy);

void pca_tree_register_x (struct pca_tree *T, double dx);

int pca_tree_draw (struct pca *P, struct pca_tree *T, const char *fname);

int pca_tree_print (struct pca_tree *T);

#endif

