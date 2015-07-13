
/*
 * pca-utils-dist.c: source code for group distance metrics.
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

/* diff_matrix_diff_mul: multiplies a difference vector against the inverse of
 * a matrix, then multiplies the resultant vector against the difference
 * vector.
 * @x1: the first vector.
 * @x2: the second vector.
 * @M: the matrix.
 */
double diff_matrix_diff_mul (struct vector *x1,
                             struct vector *x2,
                             struct matrix *M) {
  /* declare required variables. */
  struct vector *dx, *Idx;
  struct matrix *L, *I;
  double prod;

  /* make sure the input vectors and matrix are defined. */
  if (!x1 || !x2 || !M) {
    /* output an error message and return nan. */
    ferror ("input structures to energy routine are undefined");
    return NAN;
  }

  /* make sure the sizes match. */
  if (x1->n != x2->n || M->m != M->n || M->m != x1->n) {
    /* output an error message and return nan. */
    ferror ("input structure to energy routine size mismatch");
    return NAN;
  }

  /* allocate vectors for the difference calculations. */
  dx = vector_new (x1->n);
  Idx = vector_new (x1->n);

  /* make sure we allocated vectors. */
  if (!dx || !Idx) {
    /* output an error message and return nan. */
    ferror ("failed to allocate temporary vectors for energy routine");
    return NAN;
  }

  /* allocate matrices for the inversion. */
  L = matrix_new (M->m, M->n);
  I = matrix_new (M->m, M->n);

  /* make sure we allocated matrices. */
  if (!L || !I) {
    /* output an error message and return nan. */
    ferror ("failed to allocate temporary matrices for energy routine");
    return NAN;
  }

  /* calculate the difference vector. */
  if (!vector_sum (1.0, x1, -1.0, x2, dx)) {
    /* output an error message and return nan. */
    ferror ("failed to calculate vector difference");
    return NAN;
  }

  /* calculate the cholesky decomposition of the matrix. */
  if (!matrix_chol (M, L) || !matrix_chol_inv (L, I)) {
    /* output an error message and return nan. */
    ferror ("failed to cholesky invert input matrix to energy routine");
    return NAN;
  }

  /* calculate the matrix-vector product. */
  if (!matrix_vector_mul (1.0, I, dx, Idx)) {
    /* output an error message and return nothing. */
    ferror ("failed to calculate matrix-vector product in energy routine");
    return NAN;
  }

  /* calculate the inner product. */
  prod = vector_inner (dx, Idx);

  /* free the temporary structures. */
  vector_free (dx);
  vector_free (Idx);
  matrix_free (L);
  matrix_free (I);

  /* return the product. */
  return prod;
}

/* bhattacharyya_coefficient: calculates, yeah, the bhattacharyya
 * coefficient. this is not to be used outside this header.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double bhattacharyya_coefficient (struct pca_group *groupA,
                                  struct pca_group *groupB) {
  /* declare required variables. */
  double coeff, dxIdx, det, det1, det2;
  struct matrix *P, *L, *L1, *L2;

  /* allocate memory for the cholesky decompositions. */
  P = matrix_new (groupA->cov->m, groupA->cov->n);
  L = matrix_new (groupA->cov->m, groupA->cov->n);
  L1 = matrix_new (groupA->cov->m, groupA->cov->n);
  L2 = matrix_new (groupA->cov->m, groupA->cov->n);

  /* make sure we allocated the matrices. */
  if (!P || !L || !L1 || !L2) {
    /* output an error and return nan. */
    ferror ("failed to allocate cholesky decomposition matrices");
    return NAN;
  }

  /* calculate the average covariance matrix. */
  if (!matrix_sum (0.5, groupA->cov, 0.5, groupB->cov, P)) {
    /* output an error and return nan. */
    ferror ("failed to calculate average covariance matrix");
    return NAN;
  }

  /* calculate the cholesky decompositions of the matrices. */
  if (!matrix_chol (groupA->cov, L1) ||
      !matrix_chol (groupB->cov, L2) ||
      !matrix_chol (P, L)) {
    /* output an error message and return nan. */
    ferror ("failed to cholesky decompose covariances");
    return NAN;
  }

  /* calculate the determinants. */
  det = pow (matrix_chol_det (L), 0.5);
  det1 = pow (matrix_chol_det (L1), 0.25);
  det2 = pow (matrix_chol_det (L2), 0.25);

  /* calculate the energy scalar. */
  dxIdx = diff_matrix_diff_mul (groupA->mean, groupB->mean, P);

  /* calculate the coefficient. */
  coeff = ((det1 * det2) / det) * exp (-0.125 * dxIdx);

  /* free the temporary matrices. */
  matrix_free (P);
  matrix_free (L);
  matrix_free (L1);
  matrix_free (L2);

  /* return the coefficient. */
  return coeff;
}

/* pca_group_distance_euclidean: calculate the euclidean distance between
 * two pca group structures.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double pca_group_distance_euclidean (struct pca_group *groupA,
                                     struct pca_group *groupB) {
  /* declare required variables. */
  struct vector *d;
  double dnorm;

  /* allocate memory for the difference vector. */
  d = vector_new (groupA->mean->n);

  /* make sure we allocated the difference vector. */
  if (!d) {
    /* output an error message and return nan. */
    ferror ("failed to allocate difference vector");
    return NAN;
  }

  /* calculate the difference of the two mean vectors. */
  if (!vector_sum (1.0, groupA->mean, -1.0, groupB->mean, d)) {
    /* output an error message and return nan. */
    ferror ("failed to calculate difference vector");
    return NAN;
  }

  /* calculate the norm of the difference vector. */
  dnorm = vector_norm (d, 2.0);

  /* free the allocated vector. */
  vector_free (d);

  /* return the norm of the distance vector. */
  return dnorm;
}

/* pca_group_distance_mahalanobis: calculate the mahalanobis
 * distance between two pca group structures, using a pooled variance
 * covariance matrix.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double pca_group_distance_mahalanobis (struct pca_group *groupA,
                                       struct pca_group *groupB) {
  /* declare required variables. */
  double n1, n2, dist;
  struct matrix *P;

  /* allocate memory for the pooled variance-covariance matrix. */
  P = matrix_new (groupA->cov->m, groupA->cov->n);

  /* make sure we allocated the matrix. */
  if (!P) {
    /* output an error message and return nan. */
    ferror ("failed to allocate variance-covariance matrix");
    return NAN;
  }

  /* retrieve the sample sizes of the two groups. */
  n1 = (double) groupA->n_pts;
  n2 = (double) groupB->n_pts;

  /* calculate the pooled variance-covariance matrix. */
  if (!matrix_sum ((n1 - 1.0) / (n1 + n2 - 2.0), groupA->cov,
                   (n2 - 1.0) / (n1 + n2 - 2.0), groupB->cov, P)) {
    /* output an error and return nan. */
    ferror ("failed to calculate pooled variance-covariance matrix");
    return NAN;
  }

  /* calculate the distance. */
  dist = diff_matrix_diff_mul (groupA->mean, groupB->mean, P);

  /* free the pooled variance-covariance matrix. */
  matrix_free (P);

  /* return the distance. */
  return sqrt (dist);
}

/* pca_group_distance_bhattacharyya: calculate the bhattacharyya distance
 * between two pca group structures.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double pca_group_distance_bhattacharyya (struct pca_group *groupA,
                                         struct pca_group *groupB) {
  /* return the distance. */
  if (groupA == groupB)
    return 0.0;
  else
    return -1.0 * log (bhattacharyya_coefficient (groupA, groupB));
}

/* pca_group_distance_hellinger: calculate the hellinger distance
 * between two pca group structures.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double pca_group_distance_hellinger (struct pca_group *groupA,
                                     struct pca_group *groupB) {
  /* return the distance. */
  if (groupA == groupB)
    return 0.0;
  else
    return sqrt (1.0 - bhattacharyya_coefficient (groupA, groupB));
}

/* pca_group_distance_mixture: calculate mixture model pdf overlap.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double pca_group_distance_mixture (struct pca_group *groupA,
                                   struct pca_group *groupB) {
  /* declare required variables. */
  unsigned long i, j;
  double d, dx, dy;

  /* initialize the value. */
  d = 0.0;

  /* calculate the dx,dy values. */
  dx = fabs (groupA->mx->d[1] - groupA->mx->d[0]);
  dy = fabs (groupA->my->d[1] - groupA->my->d[0]);

  /* loop through the rows. */
  for (i = 0; i < groupA->map->m; i++) {
    /* and the columns. */
    for (j = 0; j < groupA->map->n; j++) {
      /* add the overlap values. */
      d += sqrt (groupA->map->d[i][j] * groupB->map->d[i][j]);
    }
  }

  /* multiply through by the differential area. */
  d *= dx * dy;

  /* return the distance value. */
  return -log (d);
}

/* pca_group_distance_chisq: calculate mixture model chi^2 values.
 * @groupA: the first group structure.
 * @groupB: the second group structure.
 */
double pca_group_distance_chisq (struct pca_group *groupA,
                                 struct pca_group *groupB) {
  /* declare required variables. */
  unsigned long i, j;
  double chisq, nA, nB, kA, kB;

  /* initialize the value. */
  chisq = 0.0;

  /* get the group counts. */
  nA = (double) groupA->n_pts;
  nB = (double) groupB->n_pts;

  /* calculate the scaling constants. */
  kA = sqrt (nB / nA);
  kB = sqrt (nA / nB);

  /* loop through the rows. */
  for (i = 0; i < groupA->map->m; i++) {
    /* and the columns. */
    for (j = 0; j < groupA->map->n; j++) {
      /* avoid singularities. */
      if (groupA->map->d[i][j] != 0.0 || groupB->map->d[i][j] != 0.0) {
        /* add the current term. */
        chisq +=
          pow (kA * nA * groupA->map->d[i][j]
             - kB * nB * groupB->map->d[i][j], 2.0)
          / (nA * groupA->map->d[i][j] + nB * groupB->map->d[i][j]);
      }
    }
  }

  /* return the distance value. */
  return chisq;
}

/* pca_distances: calculate a distance matrix of a pca structure.
 * @P: the pca structure to operate on.
 * @D: the output distance matrix.
 * @metric: the distance metric function pointer to use.
 */
int pca_distances (struct pca *P, struct matrix *D,
                   pca_distance_metric metric) {
  /* declare required variables. */
  unsigned long ga, gb;

  /* make sure the passed structures are defined. */
  if (!P || !D) {
    /* output an error message and return failure. */
    ferror ("distance matrix input structures undefined");
    return 0;
  }

  /* make sure the distance matrix size matches the group count. */
  if (D->m != P->n_groups || D->n != P->n_groups) {
    /* output an error message. */
    ferror ("distance matrix size is invalid for %lu groups",
            P->n_groups);

    /* return failure. */
    return 0;
  }

  /* make sure the groups have pdf matrices if calculating overlaps. */
  if (metric == pca_group_distance_mixture ||
      metric == pca_group_distance_chisq) {
    /* calculate the maps. */
    if (!pca_calculate_maps (P)) {
      /* output an error message and return failure. */
      ferror ("failed to build initial maps for mixture distances");
      return 0;
    }
  }

  /* loop through the groups. */
  for (ga = 0; ga < P->n_groups; ga++) {
    /* and again, obeying commutativity. */
    for (gb = ga; gb < P->n_groups; gb++) {
      /* store the calculated value. */
      D->d[ga][gb] = metric (P->groups[ga], P->groups[gb]);

      /* are we off-diagonal? */
      if (ga != gb)
        D->d[gb][ga] = D->d[ga][gb];
    }
  }

  /* return success. */
  return 1;
}

/* pca_distances_print: prints the distance matrix to standard output.
 * @P: the pca structure to operate on.
 * @D: the distance matrix to print.
 */
int pca_distances_print (struct pca *P, struct matrix *D) {
  /* declare required variables. */
  unsigned long ga, gb, wmax, i;
  char nstr[16];

  /* make sure the passed structures are defined. */
  if (!P || !D) {
    /* output an error message and return failure. */
    ferror ("distance matrix input structures undefined");
    return 0;
  }

  /* make sure the distance matrix size matches the group count. */
  if (D->m != P->n_groups || D->n != P->n_groups) {
    /* output an error message. */
    ferror ("distance matrix size is invalid for %lu groups",
            P->n_groups);

    /* return failure. */
    return 0;
  }

  /* loop through the group names. */
  for (wmax = 0, ga = 0; ga < P->n_groups; ga++) {
    /* see if the current group name length exceeds the current maximum. */
    if (strlen (P->groups[ga]->name) > wmax) {
      /* yes. store the new width. */
      wmax = strlen (P->groups[ga]->name);
    }
  }

  /* print the first spacer. */
  for (i = 0; i < wmax + 2; i++) fprintf (stdout, " ");

  /* loop through the groups. */
  for (ga = 1; ga < P->n_groups; ga++) {
    /* print the output string. */
    snprintf (nstr, 16, "%s", P->groups[ga]->name);

    /* should we add ellipsis? */
    if (strlen (nstr) > 12) {
      /* yeah, too long. */
      nstr[11] = '.';
      nstr[10] = '.';
      nstr[9] = '.';
    }

    /* print the string. */
    fprintf (stdout, "%-14s", nstr);
  }

  /* print a newline. */
  fprintf (stdout, "\n");

  /* loop through the groups. */
  for (ga = 0; ga < P->n_groups - 1; ga++) {
    /* print the group name. */
    fprintf (stdout, "%s", P->groups[ga]->name);

    /* pad the line with spaces. */
    for (i = 0; i < wmax - strlen (P->groups[ga]->name); i++) {
      /* print a single space. */
      fprintf (stdout, " ");
    }

    /* loop through the unspoken values. */
    for (gb = 1; gb <= ga; gb++) {
      /* print a placeholder. */
      for (i = 0; i < 14; i++) fprintf (stdout, " ");
    }

    /* and again, throwing commutativity away. */
    for (gb = ga + 1; gb < P->n_groups; gb++) {
      /* print the distance matrix element value. */
      fprintf (stdout, "%14le", D->d[ga][gb]);
    }

    /* print a newline. you know, since this is supposed to be readable... */
    fprintf (stdout, "\n");
  }

  /* return success. */
  return 1;
}

/* pca_distances_min: find the indices of the minimum value in a matrix.
 * @D: the input distance matrix to search.
 * @ga: a pointer to the minimum row index.
 * @gb: a pointer to the minimum column index.
 */
double pca_distances_min (struct matrix *D,
                          unsigned long *ga,
                          unsigned long *gb) {
  /* declare required variables. */
  unsigned long icur, jcur, imin, jmin;
  double dmin;

  /* make sure the passed structure is defined. */
  if (!D) {
    /* output an error message and return failure. */
    ferror ("distance matrix input structure undefined");
    *ga = *gb = 0;
    return 0;
  }

  /* initialize the minimum value. */
  imin = jmin = 0;
  dmin = DBL_MAX;

  /* loop through the distance matrix rows. */
  for (icur = 0; icur < D->m; icur++) {
    /* and though the columns, avoiding the diagonal. */
    for (jcur = icur + 1; jcur < D->n; jcur++) {
      /* see if this value is smaller than the current minimum. */
      if (D->d[icur][jcur] < dmin) {
        /* store the new indices and value. */
        imin = icur;
        jmin = jcur;
        dmin = D->d[icur][jcur];
      }
    }
  }

  /* return the indices. */
  *ga = imin;
  *gb = jmin;

  /* return the value. */
  return dmin;
}

/* pca_distances_union: unions two groups in a distance matrix.
 * @P: the pca structure to operate on.
 * @D: the distance matrix to operate on.
 * @ga: the first input group index to use.
 * @gb: the second index group index to use.
 */
struct pca_tree *pca_distances_union (struct pca *P, struct matrix *D,
                                      unsigned long ga, unsigned long gb) {
  /* declare required variables. */
  struct pca_group *gunion;
  struct pca_tree *T;
  double d, na, nb;
  unsigned long i;

  /* ensure the input structures are defined. */
  if (!P || !D) {
    /* output an error message and return nothing. */
    ferror ("input data structures passed to union are invalid");
    return NULL;
  }

  /* ensure the group indices are sorted. */
  if (ga >= gb) {
    /* output an error message and return nothing. */
    ferror ("group indices are unsorted or equal (%lu >= %lu)", ga, gb);
    return NULL;
  }

  /* build a unioned group. */
  gunion = pca_group_union_new (P, ga, gb);

  /* make sure we got a unioned group. */
  if (!gunion) {
    /* output an error message and return nothing. */
    ferror ("failed to create union of groups %lu and %lu", ga, gb);
    return NULL;
  }

  /* add the unioned group to the pca structure. */
  if (!pca_add_group (P, gunion)) {
    /* output an error message and return nothing. */
    ferror ("failed to add %s to pca structure", gunion->name);
    return NULL;
  }

  /* resize the matrix to a larger size, temporarily. */
  if (!matrix_resize (D, D->m + 1, D->n + 1)) {
    /* output an error message and return nothing. */
    ferror ("failed to increase size of distance matrix");
    return NULL;
  }

  /* allocate memory for the new parent tree node. */
  T = (struct pca_tree *) malloc (sizeof (struct pca_tree));

  /* ensure the new tree node was allocated. */
  if (!T) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate new parent tree node");
    return NULL;
  }

  /* copy the name of the union group to the tree node. */
  T->name = strdup (gunion->name);
  T->membstr = strdup (gunion->membstr);

  /* initialize the p-value of the bifurcation at the tree node. */
  T->p = 0.0;

  /* store the union group in the tree node. */
  T->group = gunion;

  /* initialize the tree linkages for later. */
  T->parent = T->child_a = T->child_b = NULL;

  /* get the point counts of the input groups. */
  na = (double) P->groups[ga]->n_pts;
  nb = (double) P->groups[gb]->n_pts;

  /* loop through the groups in the matrix. */
  for (i = 0; i < D->m - 1; i++) {
    /* calculate the new distance value. */
    d = (na / (na + nb)) * D->d[i][ga] + (nb / (na + nb)) * D->d[i][gb];

    /* store the new distance value. */
    D->d[D->m - 1][i] = d;
    D->d[i][D->m - 1] = d;
  }

  /* store the distances between the child nodes and the new node. */
  T->dist_a = D->d[D->m - 1][ga];
  T->dist_b = D->d[D->m - 1][gb];

  /* remove the input group rows and columns. */
  if (!matrix_del_row (D, ga) ||
      !matrix_del_row (D, gb - 1) ||
      !matrix_del_column (D, ga) ||
      !matrix_del_column (D, gb - 1)) {
    /* output an error message and return nothing. */
    ferror ("failed to remove input groups from distance matrix");
    return NULL;
  }

  /* delete the groups from the pca structure. */
  if (!pca_del_group (P, ga) || !pca_del_group (P, gb - 1)) {
    /* output an error message and return nothing. */
    ferror ("failed to remove input groups from pca structure");
    return NULL;
  }

  /* return the new tree node. */
  return T;
}

