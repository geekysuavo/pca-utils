
/*
 * pca-rand.c: main program source code for pca-rand.
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
 pca-rand: A tool to generate bivariate random data in list format.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-rand [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -n,--count N         Number of generated points (required)\n\
   -L,--label LBL       Group label string (required)\n\
   -H,--header          Print a list file header string (optional)\n\
   -u,--mean (U1,U2)    Group mean vector [0.0,0.0]\n\
   -v,--var (V1,V2)     Group variances [1.0,1.0]\n\
   -r,--rotation ROT    Group rotation in degrees [0.0]\n\
\n\
 The pca-rand tool generates N bivariate normally distributed data points\n\
 with unit variance using a pseudo-random number generator. This is a small\n\
 utility that is useful is generating test datasets for analysis using the\n\
 other tools in the pca-utils software package.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  struct vector *x, *y, *z, *u;
  struct matrix *V, *R, *C;
  double theta;
  char *smean, *svar, *label;
  unsigned long n;
  int i;

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

  /* get the number of points to produce. */
  n = opts_geti (OPTS_S_COUNT);

  /* ensure the point count is valid. */
  if (!n) {
    /* output an error message and exit. */
    ferror ("invalid point count (%lu)", n);
    fprintf (stdout, HELP);
    return 1;
  }

  /* get the points label. */
  label = opts_gets (OPTS_S_LABEL);

  /* ensure the points label is valid. */
  if (!label) {
    /* output an error message and exit. */
    ferror ("invalid label '%s'", label);
    fprintf (stdout, HELP);
    return 1;
  }

  /* allocate the math structures. */
  x = vector_new (2);
  y = vector_new (2);
  z = vector_new (2);
  u = vector_new (2);
  V = matrix_new (2, 2);
  R = matrix_new (2, 2);
  C = matrix_new (2, 2);

  /* ensure we allocated the math structures. */
  if (!x || !y || !z || !u || !V || !R || !C) {
    /* output an error message and exit. */
    ferror ("failed to allocate math structures");
    return 1;
  }

  /* get the mean string. */
  smean = opts_gets (OPTS_S_MEAN);

  /* did we get a mean string? */
  if (smean) {
    /* try to parse the string. */
    if (sscanf (smean, "(%lf,%lf)", &u->d[0], &u->d[1]) != 2) {
      /* output an error message and exit. */
      ferror ("failed to parse mean string '%s'", smean);
      fprintf (stdout, HELP);
      return 1;
    }
  }
  else {
    /* initialize to the origin. */
    u->d[0] = u->d[1] = 0.0;
  }

  /* get the variance string. */
  svar = opts_gets (OPTS_S_VAR);

  /* did we get a variance string? */
  if (svar) {
    /* try to parse the string. */
    if (sscanf (svar, "(%lf,%lf)", &V->d[0][0], &V->d[1][1]) != 2) {
      /* output an error message and exit. */
      ferror ("failed to parse variance string '%s'", svar);
      fprintf (stdout, HELP);
      return 1;
    }
  }
  else {
    /* initialize to unit variances. */
    V->d[0][0] = V->d[1][1] = 1.0;
  }

  /* get the rotation value. */
  theta = opts_getf (OPTS_S_ROT) / 180.0 * M_PI;

  /* build the rotation matrix. */
  R->d[0][0] = cos (theta);
  R->d[0][1] = -sin (theta);
  R->d[1][0] = sin (theta);
  R->d[1][1] = cos (theta);

  /* calculate the covariance matrix. */
  if (!matrix_matrix_mul (1.0, R, V, C)) {
    /* output an error message and exit. */
    ferror ("failed to calculate covariance matrix");
    return 1;
  }

  /* see if we should produce a header. */
  if (opts_geti (OPTS_S_HEADER)) {
    /* yep. print one. */
    fprintf (stdout,
      "Obs ID (Primary)\tObs ID (Obs. Sec. ID:1)\tM1.t[2]\tM1.t[1]\n");
  }

  /* loop a set number of times. */
  for (i = 0; i < n; i++) {
    /* calculate the two values. */
    x->d[0] = randnormal ();
    x->d[1] = randnormal ();

    /* scale the values. */
    if (!matrix_vector_mul (1.0, C, x, y)) {
      /* output an error message and exit. */
      ferror ("failed to scale random variate %d", i + 1);
      return 1;
    }

    /* translate the values. */
    if (!vector_sum (1.0, y, 1.0, u, z)) {
      /* output an error message and exit. */
      ferror ("failed to translate random variate %d", i + 1);
      return 1;
    }

    /* print the values. */
    fprintf (stdout, "%d\t%s\t%lf\t%lf\n", i + 1, label, z->d[0], z->d[1]);
  }

  /* free the math structures. */
  vector_free (x);
  vector_free (y);
  vector_free (z);
  vector_free (u);
  matrix_free (V);
  matrix_free (R);
  matrix_free (C);

  /* return success. */
  return 0;
}

