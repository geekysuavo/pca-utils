
/*
 * pca-bootstrap.c: main program source code for pca-bootstrap.
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
 pca-bootstrap: A tool to build dendrograms from scores using bootstrapping.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-bootstrap [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -k,--key             Add key symbols to leaf nodes (optional)\n\
   -i,--input FIN       Input scores list filename (required)\n\
   -o,--output FOUT     Output postscript tree file (optional)\n\
   -n,--count N         Number of bootstrap iterations [100]\n\
\n\
 The pca-bootstrap tool reads in an input PCA scores file and generates a\n\
 dendrogram of the groups using a bootstrapping method. Over N iterations,\n\
 the groups are 'rebuilt' by random sampling (with replacement) of their\n\
 points. Using the UPGMA algorithm, a set of N trees is created. A main\n\
 tree is chosen based on the original points in the dataset. Bootstrap\n\
 numbers based on the frequency of occurence of each node in the randomly\n\
 sampled trees are placed at each main tree node as percentages.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  unsigned long i, n;
  struct pca *P, *Pboot;
  struct pca_tree *T, *Tboot;
  struct matrix *D, *Dboot;
  char *ifname, *ofname;

  /* initialize the random number generator. */
  rand_init (time (NULL));

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

  /* get the number of bootstrap iterations. */
  n = opts_geti (OPTS_S_COUNT);

  /* make sure a bootstrap iteration count was provided. */
  if (!n) {
    /* use the default. */
    n = 100;
  }

  /* read in the input file. */
  P = pca_new (ifname);

  /* make sure we read in the file. */
  if (!P) {
    /* output an error message and exit. */
    ferror ("failed to read '%s'", ifname);
    return 1;
  }

  /* allocate the distance matrices. */
  D = matrix_new (P->n_groups, P->n_groups);
  Dboot = matrix_new (P->n_groups, P->n_groups);

  /* make sure we allocated the matrices. */
  if (!D || !Dboot) {
    /* output an error message and exit. */
    ferror ("failed to allocate distance matrices");
    return 1;
  }

  /* calculate the main distance matrix. */
  if (!pca_distances (P, D, pca_group_distance_euclidean)) {
    /* output an error message and exit. */
    ferror ("failed to calculate original distance matrix");
    return 1;
  }

  /* build the main tree. */
  T = pca_tree_new (P, D);

  /* make sure the main tree was calculated successfully. */
  if (!T) {
    /* output an error message and exit. */
    ferror ("failed to generate original tree");
    return 1;
  }

  /* loop a certain number of times. */
  for (i = 0; i < n; i++) {
    /* make a bootstrapped copy of the original pca structure. */
    Pboot = pca_copy_bootstrapped (P);

    /* make sure we copied the pca structure. */
    if (!Pboot) {
      /* output an error message and exit. */
      ferror ("failed to copy bootstrapped pca structure");
      return 1;
    }

    /* calculate the distance matrix. */
    if (!pca_distances (Pboot, Dboot, pca_group_distance_euclidean)) {
      /* output an error message and exit. */
      ferror ("failed to calculate bootstrapped distance matrix");
      return 1;
    }

    /* build the bootstrapped tree. */
    Tboot = pca_tree_new (Pboot, Dboot);

    /* make sure the tree generation was successful. */
    if (!Tboot) {
      /* output an error message and exit. */
      ferror ("failed to generate bootstrapped tree");
      return 1;
    }

    /* compare the bootstrap tree to the main tree. */
    pca_tree_bootstrap_compare (T, Tboot);

    /* free the bootstrapped tree structure. */
    pca_tree_free (Tboot);

    /* free the bootstrapped pca structure. */
    pca_free (Pboot);
  }

  /* finalize the bootstrap run. */
  pca_tree_bootstrap_end (T, n);

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

  /* free the distance matrices. */
  matrix_free (D);
  matrix_free (Dboot);

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

