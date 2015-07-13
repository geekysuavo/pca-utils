
/*
 * pca-overlap.c: main program source code for pca-overlap.
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
 pca-overlap: A tool to build overlap p-value tables of PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-overlap [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -i,--input FIN       Input scores list filename (required)\n\
\n\
 The pca-overlap tool reads in an input PCA scores file and generates a\n\
 symmetric matrix of p-values for null hypotheses between every group,\n\
 of which the upper triangle is displayed with group labels.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  unsigned long ga, gb;
  struct matrix *D, *Pvals;
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

  /* ensure the input filename was provided. */
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

  /* allocate distance and pvalue matrices. */
  D = matrix_new (P->n_groups, P->n_groups);
  Pvals = matrix_new (D->m, D->n);

  /* make sure we allocated the matrices. */
  if (!D || !Pvals) {
    /* output an error message and exit. */
    ferror ("failed to allocate distance matrix");
    return 1;
  }

  /* calculate the distance matrix. */
  if (!pca_distances (P, D, pca_group_distance_mahalanobis)) {
    /* output an error message and exit. */
    ferror ("failed to calculate distance matrix");
    return 1;
  }

  /* loop through the groups. */
  for (ga = 0; ga < P->n_groups; ga++) {
    /* set the diagonal. */
    Pvals->d[ga][ga] = 0.0;

    /* and again through the groups. */
    for (gb = ga + 1; gb < P->n_groups; gb++) {
      /* calculate the p-value. */
      Pvals->d[ga][gb] = mahalanobis_pvalue (D->d[ga][gb],
                                             P->groups[ga]->n_pts,
                                             P->groups[gb]->n_pts,
                                             pca_dimensions (P));

      /* also store the result in the transpose-side. */
      Pvals->d[gb][ga] = Pvals->d[ga][gb];
    }
  }

  /* print the p-value matrix. */
  pca_distances_print (P, Pvals);

  /* free the distance and pvalue matrices. */
  matrix_free (D);
  matrix_free (Pvals);

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

