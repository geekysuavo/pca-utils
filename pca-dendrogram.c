
/*
 * pca-dendrogram.c: main program source code for pca-dendrogram.
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

/* incude the pca utils header. */
#include "pca-utils.h"

/* define the help message string. */
#define HELP "\
 pca-dendrogram: A tool to build statistical trees from PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-dendrogram [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -k,--key             Add key symbols to leaf nodes (optional)\n\
   -i,--input FIN       Input scores list filename (required)\n\
   -o,--output FOUT     Output postscript tree file (optional)\n\
\n\
 The pca-dendrogram tool reads in an input PCA scores file and generates\n\
 a dendrogram for the score groups, based on a multivariate normal model\n\
 of each group and a Mahalanobis distance metric, from which an F-statistic\n\
 (and corresponding p-value) may be derived for each tree node. This tool\n\
 uses the same underlying UPGMA algorithm as pca-bootstrap, but only needs\n\
 to generate a single tree.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  struct pca_tree *T;
  struct matrix *D;
  struct pca *P;
  char *ifname, *ofname;

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

  /* read the input file. */
  P = pca_new (ifname);

  /* make sure we read in the file. */
  if (!P) {
    /* output an error message and exit. */
    ferror ("failed to read '%s'", ifname);
    return 1;
  }

  /* allocate a distance matrix. */
  D = matrix_new (P->n_groups, P->n_groups);

  /* make sure we allocated a matrix. */
  if (!D) {
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

  /* build the tree structure. */
  T = pca_tree_new (P, D);

  /* make sure we successfully built a tree. */
  if (!T) {
    /* output an error message and exit. */
    ferror ("failed to generate tree structure");
    return 1;
  }

  /* re-count the inner node point counts. */
  pca_tree_recount_points (T);

  /* generate p-values on the tree nodes. */
  if (!pca_tree_assign_pvalues (T, pca_dimensions (P),
                                pca_group_distance_mahalanobis)) {
    /* output an error message and exit. */
    ferror ("failed to calculate tree node p-values");
    return 1;
  }

  /* get the output filename. */
  ofname = opts_gets (OPTS_S_OUTPUT);

  /* see if an output filename was passed. */
  if (ofname) {
    /* output the final tree to a postscript file. */
    if (!pca_tree_draw (P, T, ofname)) {
      /* output an error message and exit. */
      ferror ("failed to draw tree to '%s'", ofname);
      return 1;
    }
  }
  else {
    /* output the final tree to standard out. */
    if (!pca_tree_print (T)) {
      /* output an error message and exit. */
      ferror ("failed to print tree to stdout");
      return 1;
    }
  }

  /* free the tree structure. */
  pca_tree_free (T);

  /* free the distance matrix. */
  matrix_free (D);

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

