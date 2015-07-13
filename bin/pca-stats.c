
/*
 * pca-stats.c: main program source code for pca-stats.
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

/* include the pca utils header. */
#include "pca-utils.h"

/* define the help message string. */
#define HELP "\
 pca-stats: A tool to display simple statistics of PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-stats [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -i,--input FIN       Input scores list filename (required)\n\
\n\
 The pca-stats tool reads in an input PCA scores file and generates basic\n\
 statistics for each group in the file. The statistics generated include\n\
 mean, covariance, covariance eigenvalues, covariance eigenvectors and J2.\n\
 This is mostly useful as a debug tool to validate other tools in the\n\
 pca-utils software package.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  struct matrix *M0, *C0, *V0;
  struct vector *lambda0;
  unsigned long g, i;
  double det0, J2;
  struct pca *P;
  char *ifname;

  /* parse the command-line arguments. */
  if (!opts_init (argc, argv)) {
    /* output an error message and exit. */
    ferror ("failed to parse arguments");
    fprintf (stdout, HELP);
    return 1;
  }

  /* see if we should display the help message. */
  if (opts_geti (OPTS_S_HELP)) {
    /* print the help message and quit. */
    fprintf (stdout, HELP);
    return 0;
  }

  /* get the input filename. */
  ifname = opts_gets (OPTS_S_INPUT);

  /* ensure the input filename was given. */
  if (!ifname) {
    /* output an error message and exit. */
    ferror ("input file required");
    fprintf (stdout, HELP);
    return 1;
  }

  /* read in the input file. */
  P = pca_new (ifname);

  /* make sure we read in the file. */
  if (!P) {
    /* output an error message and exit. */
    ferror ("failed to read '%s'", ifname);
    return 1;
  }

  /* combine the groups' points into a single data matrix. */
  M0 = pca_datamatrix (P);

  /* make sure that worked. */
  if (!M0) {
    /* output an error message and exit. */
    ferror ("failed to build combined group scores matrix");
    return 1;
  }

  /* calculate the covariance of the data matrix. */
  C0 = matrix_cov (M0);

  /* ensure the calculated a covariance. */
  if (!C0) {
    /* output an error message and return failure. */
    ferror ("failed to calculate combined group covariance matrix");
    return 1;
  }

  /* allocate matrices for the matrix eigendecomposition. */
  V0 = matrix_new (pca_dimensions (P), pca_dimensions (P));
  lambda0 = vector_new (pca_dimensions (P));

  /* calculate the eigendecomposition of the entire matrix. */
  if (!M0 || !lambda0 || !V0 || !matrix_eig (C0, lambda0, V0)) {
    /* output an error message and exit. */
    ferror ("failed to eigendecompose combined scores matrix");
    return 1;
  }

  /* calculate the matrix determinant. all that work for a determinant? */
  det0 = matrix_eig_det (lambda0);

  /* loop through the groups. */
  for (g = 0; g < P->n_groups; g++) {
    /* print the mean vector. */
    fprintf (stdout, "mean('%s'):\n", P->groups[g]->name);
    vector_print (P->groups[g]->mean);

    /* print the covariance matrix. */
    fprintf (stdout, "\ncov('%s'):\n", P->groups[g]->name);
    matrix_print (P->groups[g]->cov);

    /* print the eigenvalue vector. */
    fprintf (stdout, "\nlambda('%s'):\n", P->groups[g]->name);
    vector_print (P->groups[g]->eig_l);

    /* print the eigenvector matrix. */
    fprintf (stdout, "\nQ('%s'):\n", P->groups[g]->name);
    matrix_print (P->groups[g]->eig_v);

    /* print the J2 statistic. */
    J2 = det0 / matrix_eig_det (P->groups[g]->eig_l);
    fprintf (stdout, "\nJ2('%s'): %lf\n", P->groups[g]->name, J2);

    /* print a newline. */
    fprintf (stdout, "\n");

    /* should we print a separator? */
    if (g < P->n_groups - 1) {
      fprintf (stdout, "  ");
      for (i = 0; i < 76; i++) fprintf (stdout, "-");
      fprintf (stdout, "\n\n");
    }
  }

  /* free the eigendecomposition structures. */
  matrix_free (M0);
  matrix_free (C0);
  matrix_free (V0);
  vector_free (lambda0);

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

