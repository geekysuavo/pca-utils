
/*
 * pca-distances.c: main program source code for pca-distances.
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
 pca-distances: A tool to build distance matrices of PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-distances [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -i,--input FIN       Input scores list filename (required)\n\
   -m,--metric DM       One of the supported distance metrics [EUC]\n\
   -D,--divisions DIV   Number of divisions for MIX, CHI only [128]\n\
\n\
 Metrics:\n\
   EUC: Euclidean\n\
        This is simply the square root of the dot product of the vector\n\
        of displacements between the group means.\n\
   MAH: Mahalanobis\n\
        This is a weighted metric, where the inverse of the pooled\n\
        variance-covariance matrix of the two groups is placed into\n\
        the multiplication.\n\
   BHA: Bhattacharyya\n\
        This is another weighted metric, based on the negative natural\n\
        logarithm of the Bhattacharyya coefficient for the two groups.\n\
   HEL: Hellinger\n\
        Similar to the Bhattacharyya distance, this is the square root\n\
        of one minus the Bhattacharyya coefficient. Unlike BHA, HEL obeys\n\
        the triangle inequality.\n\
   MIX: Mixture\n\
        Builds a Bhattacharyya distance by approximating the overlap\n\
        integral between groups in 2D space using Riemann sums. The\n\
        probability density of each group is estimated by convolving\n\
        every point in the group with a Gaussian point spread function\n\
        and normalizing to unit volume.\n\
   CHI: Chi-square\n\
        Using the same convolution as MIX, this metric calculates a chi-\n\
        square value for each pair of groups. Note that this value will\n\
        not necessarily vary as a chi-square distribution, since the\n\
        vast majority of compared pixels will be zero.\n\
\n\
 The pca-distances tool reads in an input PCA scores file and generates a\n\
 symmetric matrix of intergroup distances, using the distance metric DM,\n\
 of which it prints only the upper triangle, with nice group labels.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  pca_distance_metric metric;
  char *ifname, *metstr;
  struct matrix *D;
  struct pca *P;

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

  /* set the default distance metric. */
  metric = pca_group_distance_euclidean;

  /* get the distance metric string. */
  metstr = opts_gets (OPTS_S_METRIC);

  /* see if the distance metric was provided. */
  if (metstr) {
    /* check the distance metric argument. */
    if (strncmp (metstr, "EUC", 3) == 0 ||
        strncmp (metstr, "euc", 3) == 0) {
      /* set the euclidean metric. */
      metric = pca_group_distance_euclidean;
    }
    else if (strncmp (metstr, "MAH", 3) == 0 ||
             strncmp (metstr, "mah", 3) == 0) {
      /* set the mahalanobis distance. */
      metric = pca_group_distance_mahalanobis;
    }
    else if (strncmp (metstr, "BHA", 3) == 0 ||
             strncmp (metstr, "bha", 3) == 0) {
      /* set the bhattacharyya distance. */
      metric = pca_group_distance_bhattacharyya;
    }
    else if (strncmp (metstr, "HEL", 3) == 0 ||
             strncmp (metstr, "hel", 3) == 0) {
      /* set the hellinger distance. */
      metric = pca_group_distance_hellinger;
    }
    else if (strncmp (metstr, "MIX", 3) == 0 ||
             strncmp (metstr, "mix", 3) == 0) {
      /* set the mixture distance. */
      metric = pca_group_distance_mixture;
    }
    else if (strncmp (metstr, "CHI", 3) == 0 ||
             strncmp (metstr, "chi", 3) == 0) {
      /* set the mixture distance. */
      metric = pca_group_distance_chisq;
    }
    else {
      /* output an error message and exit. */
      ferror ("unrecognized distance metric '%s'", metstr);
      fprintf (stdout, HELP);
      return 1;
    }
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
  if (!pca_distances (P, D, metric)) {
    /* output an error message and exit. */
    ferror ("failed to calculate distance matrix");
    return 1;
  }

  /* write out the distance matrix. */
  pca_distances_print (P, D);

  /* free the distance matrix. */
  matrix_free (D);

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

