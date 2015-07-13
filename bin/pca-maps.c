
/*
 * pca-maps.c: main program source code for pca-maps.
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
 pca-maps: A tool to calculate mixture maps from PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-maps [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -i,--input FIN       Input scores list filename (required)\n\
   -D,--divisions DIV   Number of mixture map divisions [128]\n\
\n\
 The pca-maps tool reads in an input PCA scores file and generates\n\
 a table of estimated probability density functions for each group.\n\
 Honestly, it's really just useful as a debug tool.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  unsigned long g, i, j;
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

  /* read the input file. */
  P = pca_new (ifname);

  /* make sure we read in the file. */
  if (!P) {
    /* output an error message and exit. */
    ferror ("failed to read '%s'", ifname);
    return 1;
  }

  /* calculate the map. */
  if (!pca_calculate_maps (P)) {
    /* output an error message and exit. */
    ferror ("failed to calculate density maps");
    return 1;
  }

  /* loop through the groups. */
  for (g = 0; g < P->n_groups; g++) {
    /* loop through the map rows. */
    for (i = 0; i < P->groups[g]->map->m; i++) {
      /* loop through the map columns. */
      for (j = 0; j < P->groups[g]->map->n; j++) {
        /* output the value. */
        fprintf (stdout, "%lu %le %le %le\n", g,
          P->groups[g]->mx->d[i],
          P->groups[g]->my->d[j],
          P->groups[g]->map->d[i][j]);
      }
    }
  }

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

