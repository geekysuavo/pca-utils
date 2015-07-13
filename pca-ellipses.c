
/*
 * pca-ellipses.c: main program source code for pca-ellipses.
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

/* define constants for the plot layout. */
#define PLOT_W 800.0
#define PLOT_H 540.0
#define PLOT_PAD 25.0

/* define macros to transform scores space into plot space. */
#define XCOORD(a,max) \
  (3.0 * PLOT_PAD + ((PLOT_W - 6.0 * PLOT_PAD) / (2.0 * max)) * (a + max))
#define YCOORD(a,max) \
  (3.0 * PLOT_PAD + ((PLOT_H - 6.0 * PLOT_PAD) / (2.0 * max)) * (a + max))

/* define the help message string. */
#define HELP "\
 pca-ellipses: A tool to build confidence ellipses of PCA group scores.\n\
 Copyright (C) 2015 Bradley Worley. Released under the GNU GPL 3.0.\n\
\n\
 Usage:\n\
   pca-ellipses [arguments]\n\
\n\
 Arguments:\n\
   -h,--help            Display this help message\n\
   -k,--key             Add a key in the upper right (optional)\n\
   -i,--input FIN       Input scores list filename (required)\n\
   -o,--output FOUT     Output postscript tree file (optional)\n\
   -1,--component1 PC1  Percent contribution for PC1 (optional)\n\
   -2,--component2 PC2  Percent contribution for PC2 (optional)\n\
\n\
 The pca-ellipses tool reads in an input PCA scores file and generates\n\
 a confidence ellipse for each group based on estimated two-sigma variance\n\
 along its pricipal axes. When no output file is specified, pca-ellipses\n\
 outputs confidence ellipse data that is directly readable by the gnuplot\n\
 program. If an output filename is specified, pca-ellipses generates a\n\
 postscript plot of the scores.\n\
\n\
"

/* main: application entry point.
 * @argc: the number of commandline arguments.
 * @argv: the commandline argument string array.
 */
int main (int argc, char **argv) {
  /* declare required variables. */
  double mu0, mu1, l0, l1, v00, v01, v10, v11, r0, r1;
  double t, x, y, xmax, ymax, xtic, ytic, ex, ey, ssx, ssy, pc1, pc2;
  double red, green, blue;
  struct matrix *ellipse, *total;
  struct list *ellipses, *iter;
  char *ifname, *ofname, nstr[64];
  unsigned long g, i;
  struct pca *P;
  FILE *fh;

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
  if (P->n_indexes != 2) {
    /* output an error message and exit. */
    ferror ("ellipse calculation is only supported for two dimensions");
    fprintf (stdout, HELP);
    return 1;
  }

  /* initialize the plot boundaries. */
  xmax = ymax = 0.0;

  /* calculate the variances of the entire dataset. */
  pca_calculate_var (P, &ssx, &ssy);

  /* allocate an ellipse matrix for the total ellipse. */
  total = matrix_new (1024, 2);

  /* make sure we allocated an ellipse. */
  if (!total) {
    /* output an error message and exit. */
    ferror ("failed to allocate matrix for total ellipse");
    return 1;
  }

  /* add points to the total ellipse. */
  for (i = 0, t = 0.0; i < 1024; i++, t += M_PI / 512.0) {
    /* add the point into the ellipse. */
    ex = total->d[i][0] = sqrt (ssx * CHISQ2) * cos (t);
    ey = total->d[i][1] = sqrt (ssy * CHISQ2) * sin (t);

    /* check for new maxima. */
    if (fabs (ex) > xmax) xmax = fabs (ex);
    if (fabs (ey) > ymax) ymax = fabs (ey);
  }

  /* allocate memory for the ellipse list. */
  ellipses = list_new (NULL);

  /* make sure we allocated a list. */
  if (!ellipses) {
    /* output an error message and exit. */
    ferror ("failed to allocate ellipse storage list");
    return 1;
  }

  /* loop through the groups in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* create temporary storage variables for mathematical values. */
    mu0 = P->groups[g]->mean->d[0];
    mu1 = P->groups[g]->mean->d[1];
    l0 = P->groups[g]->eig_l->d[0];
    l1 = P->groups[g]->eig_l->d[1];
    v00 = P->groups[g]->eig_v->d[0][0];
    v01 = P->groups[g]->eig_v->d[0][1];
    v10 = P->groups[g]->eig_v->d[1][0];
    v11 = P->groups[g]->eig_v->d[1][1];

    /* calculate the radii of the ellipse. */
    r0 = sqrt (l0 * CHISQ2);
    r1 = sqrt (l1 * CHISQ2);

    /* allocate the next ellipse. */
    ellipse = matrix_new (1024, 2);

    /* make sure we allocated the ellipse. */
    if (!ellipse) {
      /* output an error message and exit. */
      ferror ("failed to allocate matrix for ellipse %lu", g);
      return 1;
    }

    /* loop over the parameter to draw the ellipse. */
    for (i = 0, t = 0.0; i < 1024; i++, t += M_PI / 512.0) {
      /* calculate the x,y pair. */
      ex = mu0 + r0 * v00 * cos (t) + r1 * v01 * sin (t);
      ey = mu1 + r0 * v10 * cos (t) + r1 * v11 * sin (t);

      /* store the x,y pair. */
      ellipse->d[i][0] = ex;
      ellipse->d[i][1] = ey;

      /* check for new maxima. */
      if (fabs (ex) > xmax) xmax = fabs (ex);
      if (fabs (ey) > ymax) ymax = fabs (ey);
    }

    /* store the ellipse into the list. */
    list_append (ellipses, ellipse);
  }

  /* get the output filename. */
  ofname = opts_gets (OPTS_S_OUTPUT);

  /* was an output filename specified? */
  if (ofname) {
    /* open the output postscript file. */
    fh = draw_document_open (ofname, PLOT_W, PLOT_H);

    /* ensure we opened the document. */
    if (!fh) {
      /* output an error message and exit. */
      ferror ("failed to open '%s' for drawing", ofname);
      return 1;
    }

    /* draw the plot boundaries. */
    draw_newpath (fh);
    draw_moveto (fh, 3.0 * PLOT_PAD, 3.0 * PLOT_PAD);
    draw_lineto (fh, 3.0 * PLOT_PAD, PLOT_H - 3.0 * PLOT_PAD);
    draw_lineto (fh, PLOT_W - 3.0 * PLOT_PAD, PLOT_H - 3.0 * PLOT_PAD);
    draw_lineto (fh, PLOT_W - 3.0 * PLOT_PAD, 3.0 * PLOT_PAD);
    draw_lineto (fh, 3.0 * PLOT_PAD, 3.0 * PLOT_PAD);
    draw_closepath (fh);
    draw_stroke (fh);

    /* draw the plot crosshairs. */
    draw_line (fh,
      3.0 * PLOT_PAD, PLOT_H / 2.0,
      PLOT_W - 3.0 * PLOT_PAD, PLOT_H / 2.0);
    draw_line (fh,
      PLOT_W / 2.0, 3.0 * PLOT_PAD,
      PLOT_W / 2.0, PLOT_H - 3.0 * PLOT_PAD);

    /* get a good tick spacing. */
    for (xtic = 1.0; xmax / xtic > 5.0; xtic += 1.0);
    for (ytic = 1.0; ymax / ytic > 5.0; ytic += 1.0);

    /* loop through the x tick marks. */
    for (x = xtic; x <= xmax; x += xtic) {
      /* draw the positive tick. */
      draw_line (fh,
        XCOORD (x, xmax), 3.0 * PLOT_PAD,
        XCOORD (x, xmax), 3.0 * PLOT_PAD - 6.0);

      /* print the positive tick label. */
      snprintf (nstr, 64, "%.1lf", x);
      draw_moveto (fh,
        XCOORD (x, xmax) - 4.0 * strlen (nstr),
        3.0 * PLOT_PAD - 16.0);
      draw_text (fh, 12, nstr);

      /* draw the negative tick. */
      draw_line (fh,
        XCOORD (-x, xmax), 3.0 * PLOT_PAD,
        XCOORD (-x, xmax), 3.0 * PLOT_PAD - 6.0);

      /* print the negative tick label. */
      snprintf (nstr, 64, "%.1lf", -x);
      draw_moveto (fh,
        XCOORD (-x, xmax) - 4.0 * strlen (nstr),
        3.0 * PLOT_PAD - 16.0);
      draw_text (fh, 12, nstr);
    }

    /* draw the x zero tick. */
    draw_line (fh,
      XCOORD (0.0, xmax), 3.0 * PLOT_PAD,
      XCOORD (0.0, xmax), 3.0 * PLOT_PAD - 6.0);

    /* print the x zero tick label. */
    snprintf (nstr, 64, "%.1lf", 0.0);
    draw_moveto (fh,
      XCOORD (0.0, xmax) - 4.0 * strlen (nstr),
      3.0 * PLOT_PAD - 16.0);
    draw_text (fh, 12, nstr);

    /* loop through the y tick marks. */
    for (y = ytic; y <= ymax; y += ytic) {
      /* draw the positive tick. */
      draw_line (fh,
        3.0 * PLOT_PAD, YCOORD (y, ymax),
        3.0 * PLOT_PAD - 6.0, YCOORD (y, ymax));

      /* print the positive tick label. */
      snprintf (nstr, 64, "%.1lf", y);
      draw_moveto (fh,
        3.0 * PLOT_PAD - 18.0 - 5.0 * strlen (nstr),
        YCOORD (y, ymax) - 4.0);
      draw_text (fh, 12, nstr);

      /* draw the negative tick. */
      draw_line (fh,
        3.0 * PLOT_PAD, YCOORD (-y, ymax),
        3.0 * PLOT_PAD - 6.0, YCOORD (-y, ymax));

      /* print the negative tick label. */
      snprintf (nstr, 64, "%.1lf", -y);
      draw_moveto (fh,
        3.0 * PLOT_PAD - 18.0 - 5.0 * strlen (nstr),
        YCOORD (-y, ymax) - 4.0);
      draw_text (fh, 12, nstr);
    }

    /* draw the y zero tick. */
    draw_line (fh,
      3.0 * PLOT_PAD, YCOORD (0.0, ymax),
      3.0 * PLOT_PAD - 6.0, YCOORD (0.0, ymax));

    /* print the y zero tick label. */
    snprintf (nstr, 64, "%.1lf", 0.0);
    draw_moveto (fh,
      3.0 * PLOT_PAD - 18.0 - 5.0 * strlen (nstr),
      YCOORD (0.0, ymax) - 4.0);
    draw_text (fh, 12, nstr);

    /* get the component-1 contribution value. */
    pc1 = opts_getf (OPTS_S_PC1);

    /* did we get a value? */
    if (pc1) {
      /* cool, print a valid axis label. */
      snprintf (nstr, 64, "PC1 (%.1lf%%)", pc1);
    }
    else {
      /* damn. just do something stupid. */
      snprintf (nstr, 64, "x");
    }

    /* print the x axis label. */
    draw_moveto (fh,
      XCOORD (0.0, xmax) - 4.0 * strlen (nstr),
      2.0 * PLOT_PAD - 16.0);
    draw_text (fh, 16, nstr);

    /* get the component-2 contribution value. */
    pc2 = opts_getf (OPTS_S_PC2);

    /* did we get a value? */
    if (pc2) {
      /* cool, print a valid axis label. */
      snprintf (nstr, 64, "PC2 (%.1lf%%)", pc2);
    }
    else {
      /* damn. just do something stupid. */
      snprintf (nstr, 64, "y");
    }

    /* print the y axis label. */
    draw_moveto (fh,
      2.0 * PLOT_PAD - 30.0,
      YCOORD (0.0, ymax) - 4.0 * strlen (nstr));
    fprintf (fh, "90 rotate\n");
    draw_text (fh, 16, nstr);
    fprintf (fh, "-90 rotate\n");

    /* loop through the ellipses. */
    for (iter = ellipses->next; iter; iter = iter->next) {
      /* gain a reference to the ellipse. */
      ellipse = (struct matrix *) iter->d;

      /* start a new path for the ellipse. */
      draw_newpath (fh);
      draw_moveto (fh,
        XCOORD (ellipse->d[0][0], xmax),
        YCOORD (ellipse->d[0][1], ymax));

      /* loop through the points of the ellipse. */
      for (i = 1; i < ellipse->m; i++) {
        /* draw the next point. */
        draw_lineto (fh,
          XCOORD (ellipse->d[i][0], xmax),
          YCOORD (ellipse->d[i][1], ymax));
      }

      /* close the ellipse path. */
      draw_closepath (fh);
      draw_stroke (fh);
    }

    /* open a new path for the total ellipse. */
    draw_newpath (fh);
    draw_moveto (fh,
      XCOORD (total->d[0][0], xmax),
      YCOORD (total->d[0][1], ymax));

    /* loop through the points of the ellipse. */
    for (i = 1; i < total->m; i++) {
      /* draw the next point. */
      draw_lineto (fh,
        XCOORD (total->d[i][0], xmax),
        YCOORD (total->d[i][1], ymax));
    }

    /* close the total ellipse path. */
    draw_closepath (fh);
    draw_stroke (fh);

    /* loop through the groups. */
    for (g = 0; g < P->n_groups; g++) {
      /* does the user want to draw the key? */
      if (opts_geti (OPTS_S_KEY)) {
        /* print a key symbol. */
        draw_symbol (fh, g, P->n_groups,
          XCOORD (xmax, xmax) - 10.0,
          YCOORD (ymax, ymax) - 10.0 * (g + 1));

        /* print the group name. */
        draw_moveto (fh,
          XCOORD (xmax, xmax) - 60.0,
          YCOORD (ymax, ymax) - 10.0 * g - 12.0);
        draw_text (fh, 10, P->groups[g]->name);
      }

      /* loop through the points. */
      for (i = 0; i < P->groups[g]->n_pts; i++) {
        /* draw the currently indexed point. */
        draw_symbol (fh, g, P->n_groups,
          XCOORD (P->groups[g]->pts[i]->d[0], xmax),
          YCOORD (P->groups[g]->pts[i]->d[1], ymax));
      }
    }

    /* close the postscript file. */
    draw_document_close (fh);
  }
  else {
    /* print the first component contribution. */
    pc1 = opts_getf (OPTS_S_PC1);
    if (pc1)
      printf ("set xlabel 'PC1 (%.2lf%%)'\n", pc1);
    else
      printf ("set xlabel 'x'\n");

    /* print the second component contribution. */
    pc2 = opts_getf (OPTS_S_PC2);
    if (pc2)
      printf ("set ylabel 'PC2 (%.2lf%%)'\n", pc2);
    else
      printf ("set ylabel 'y'\n");

    /* set required parameters for plotting ellipses. */
    printf ("set parametric\n");
    printf ("set trange [0 : 2 * pi]\n");

    /* define the functions used to plot ellipses. */
    printf ("c(r,t) = r * cos(t)\n");
    printf ("s(r,t) = r * sin(t)\n");
    printf ("R(mu,k1,k2,ra,rb,t) = ");
    printf ("mu + k1 * c(ra,t) + k2 * s(rb,t)\n");

    /* start plotting stuff. */
    printf ("plot \\\n");

    /* output the ellipse for the main ellipse. */
    printf ("  R(%le,%le,%le,%le,%le,t), \\\n",
      0.0, 1.0, 0.0, sqrt (ssx * CHISQ2), sqrt (ssy * CHISQ2));
    printf ("  R(%le,%le,%le,%le,%le,t) notitle \\\n",
      0.0, 0.0, 1.0, sqrt (ssx * CHISQ2), sqrt (ssy * CHISQ2));

    /* finish the main ellipse plot statement. */
    printf ("  linetype rgb '#000000', \\\n");

    /* loop through the groups in the pca structure. */
    for (g = 0; g < P->n_groups; g++) {
      /* get the group color. */
      draw_rgb (g, P->n_groups, &red, &green, &blue);

      /* output the ellipse for this group. */
      printf ("  R(%le,%le,%le,%le,%le,t), \\\n",
        P->groups[g]->mean->d[0],
        P->groups[g]->eig_v->d[0][0],
        P->groups[g]->eig_v->d[0][1],
        sqrt (P->groups[g]->eig_l->d[0] * CHISQ2),
        sqrt (P->groups[g]->eig_l->d[1] * CHISQ2));
      printf ("  R(%le,%le,%le,%le,%le,t) notitle \\\n",
        P->groups[g]->mean->d[1],
        P->groups[g]->eig_v->d[1][0],
        P->groups[g]->eig_v->d[1][1],
        sqrt (P->groups[g]->eig_l->d[0] * CHISQ2),
        sqrt (P->groups[g]->eig_l->d[1] * CHISQ2));

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
      printf ("pt 7 ps 2 lc rgb '#%02x%02x%02x'%s\n",
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
        printf ("  %le %le\n",
          P->groups[g]->pts[i]->d[0],
          P->groups[g]->pts[i]->d[1]);
      }

      /* end the point set. */
      printf ("end\n");
    }
  }

  /* free the total ellipse. */
  matrix_free (total);

  /* free the ellipse list. */
  list_foreach (ellipses, (void (*) (void *)) &matrix_free);
  list_free (ellipses);

  /* free the pca structure. */
  pca_free (P);

  /* return success. */
  return 0;
}

