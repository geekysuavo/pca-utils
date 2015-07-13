
/*
 * pca-ellipsoids.c: main program source code for pca-ellipsoids.
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
 pca-ellipsoids: A tool to build confidence ellipsoids of PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-ellipses [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -i,--input FIN       Input scores list filename (required)\n\
   -1,--component1 PC1  Percent contribution for PC1 (optional)\n\
   -2,--component2 PC2  Percent contribution for PC2 (optional)\n\
   -3,--component3 PC3  Percent contribution for PC3 (optional)\n\
\n\
 The pca-ellipsoids tool reads in an input PCA scores file and generates\n\
 a confidence ellipse for each group based on estimated two-sigma variance\n\
 along its pricipal axes. The pca-ellipsoids tool generates a ready to\n\
 use gnuplot script for generating 3D plots of scores data.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  double pc1, pc2, pc3, red, green, blue;
  unsigned long g, i;
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

  /* make sure we're on a drawable set of dimensions. */
  if (P->n_indexes != 3) {
    /* output an error message and exit. */
    ferror ("ellipsoid calculation is only supported for three dimensions");
    fprintf (stdout, HELP);
    return 1;
  }

  /* get the first component contribution. */
  pc1 = opts_getf (OPTS_S_PC1);
  if (pc1)
    printf ("set xlabel 'PC1 (%.2lf%%)'\n", pc1);
  else
    printf ("set xlabel 'x'\n");

  /* get the second component contribution. */
  pc2 = opts_getf (OPTS_S_PC2);
  if (pc2)
    printf ("set ylabel 'PC2 (%.2lf%%)'\n", pc2);
  else
    printf ("set ylabel 'y'\n");

  /* get the third component contribution. */
  pc3 = opts_getf (OPTS_S_PC3);
  if (pc3)
    printf ("set zlabel 'PC3 (%.2lf%%)'\n", pc3);
  else
    printf ("set zlabel 'z'\n");

  /* set required parameters for plotting ellipsoids. */
  printf ("set view 65, 30\n");
  printf ("set parametric\n");
  printf ("set urange [0 : 2 * pi]\n");
  printf ("set vrange [0 : 2 * pi]\n");

  /* define the functions used to plot ellipsoids. */
  printf ("a(r,u,v) = r * cos(u) * cos(v)\n");
  printf ("b(r,u,v) = r * sin(u) * cos(v)\n");
  printf ("c(r,u,v) = r * sin(v)\n");
  printf ("R(mu,k1,k2,k3,ra,rb,rc,u,v) = ");
  printf ("mu + k1 * a(ra,u,v) + k2 * b(rb,u,v) + k3 * c(rc,u,v)\n");

  /* start plotting stuff. */
  printf ("splot \\\n");

  /* loop through the groups in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* get the group color. */
    draw_rgb (g, P->n_groups, &red, &green, &blue);

    /* output the ellipsoid for this group. */
    printf ("  R(%le,%le,%le,%le,%le,%le,%le,u,v), \\\n",
      P->groups[g]->mean->d[0],
      P->groups[g]->eig_v->d[0][0],
      P->groups[g]->eig_v->d[0][1],
      P->groups[g]->eig_v->d[0][2],
      sqrt (P->groups[g]->eig_l->d[0] * CHISQ3),
      sqrt (P->groups[g]->eig_l->d[1] * CHISQ3),
      sqrt (P->groups[g]->eig_l->d[2] * CHISQ3));
    printf ("  R(%le,%le,%le,%le,%le,%le,%le,u,v), \\\n",
      P->groups[g]->mean->d[1],
      P->groups[g]->eig_v->d[1][0],
      P->groups[g]->eig_v->d[1][1],
      P->groups[g]->eig_v->d[1][2],
      sqrt (P->groups[g]->eig_l->d[0] * CHISQ3),
      sqrt (P->groups[g]->eig_l->d[1] * CHISQ3),
      sqrt (P->groups[g]->eig_l->d[2] * CHISQ3));
    printf ("  R(%le,%le,%le,%le,%le,%le,%le,u,v) notitle \\\n",
      P->groups[g]->mean->d[2],
      P->groups[g]->eig_v->d[2][0],
      P->groups[g]->eig_v->d[2][1],
      P->groups[g]->eig_v->d[2][2],
      sqrt (P->groups[g]->eig_l->d[0] * CHISQ3),
      sqrt (P->groups[g]->eig_l->d[1] * CHISQ3),
      sqrt (P->groups[g]->eig_l->d[2] * CHISQ3));

    /* finish the group plot statement. */
    printf ("  linetype rgb '#%02x%02x%02x', \\\n",
      (unsigned int) (255 * red),
      (unsigned int) (255 * green),
      (unsigned int) (255 * blue));
  }

  /* loop again through the groups. */
  for (g = 0; g < P->n_groups; g++) {
    /* get the group color. */
    draw_rgb (g, P->n_groups, &red, &green, &blue);

    /* print the plot statement for the group points. */
    printf ("  '-' with points title '%s' ",
      P->groups[g]->name);

    /* finish the plot statement. */
    printf ("pt 7 ps 3 lc rgb '#%02x%02x%02x'%s\n",
      (unsigned int) (200 * red),
      (unsigned int) (200 * green),
      (unsigned int) (200 * blue),
      g == P->n_groups - 1 ? "" : ", \\");
  }

  /* loop one final time through the groups. */
  for (g = 0; g < P->n_groups; g++) {
    /* loop through the points. */
    for (i = 0; i < P->groups[g]->n_pts; i++) {
      /* print the point. */
      printf ("  %le %le %le\n",
        P->groups[g]->pts[i]->d[0],
        P->groups[g]->pts[i]->d[1],
        P->groups[g]->pts[i]->d[2]);
    }

    /* end the point set. */
    printf ("end\n");
  }

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

