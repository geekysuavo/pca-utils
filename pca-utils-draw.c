
/*
 * pca-utils-draw.c: source code for postscript drawing.
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

/* define the postscript header to all rendered documents. */
#define DRAW_HEADER "\
%%!PS-Adobe-2.0\n\
%%%%Title: pca-utils generated %04d-%02d-%02d %02d:%02d:%02d\n\
%%%%DocumentFonts: (atend)\n\
%%%%Pages: 1 1\n\
%%%%BoundingBox: 0 0 %lu %lu\n\
%%%%EndComments\n\
/l {newpath moveto lineto stroke} def\n\
%%%%EndProlog\n\
%%%%\n\
%%%%Page: 1 1\n\
%%%%PageBoundingBox: 0 0 %lu %lu\n\
%%%%PageFonts: (atend)\n\
%%%%BeginPageSetup\n\
"

/* define the postscript footer to all rendered documents. */
#define DRAW_FOOTER "\
showpage\n\
%%%%PageTrailer\n\
%%%%PageFonts: Courier\n\
%%%%Trailer\n\
%%%%DocumentFonts: Courier\n\
"

/* draw_document_open: prints a postscript header to a file handle.
 * @fname: the output filename.
 * @w: the document width.
 * @h: the document height.
 */
FILE *draw_document_open (const char *fname,
                          unsigned long w,
                          unsigned long h) {
  /* declare required variables. */
  struct tm *timestamp;
  time_t tsraw;
  FILE *fh;

  /* try to open the output file. */
  fh = fopen (fname, "wb");

  /* make sure we opened the file. */
  if (!fh) {
    /* output an error and return nothing. */
    ferror ("failed to open '%s' for writing", fname);
    return NULL;
  }

  /* grab the timestamp information. */
  time (&tsraw);
  timestamp = localtime (&tsraw);

  /* output the document header. */
  fprintf (fh, DRAW_HEADER,
    timestamp->tm_year + 1900,
    timestamp->tm_mon + 1,
    timestamp->tm_mday,
    timestamp->tm_hour,
    timestamp->tm_min,
    timestamp->tm_sec,
    w + 30, h + 30,
    w, h);

  /* return the output file handle. */
  return fh;
}

/* draw_document_close: prints a postscript footer to a file handle.
 * @fh: the file handle to print to and close.
 */
void draw_document_close (FILE *fh) {
  /* output the document footer. */
  fprintf (fh, DRAW_FOOTER);

  /* close the output file. */
  fclose (fh);
}

/* draw_newpath: opens a new path in the postscript file.
 * @fh: the file handle to print to.
 */
void draw_newpath (FILE *fh) {
  /* print the newpath statement. */
  fprintf (fh, "newpath\n");
}

/* draw_closepath: closes a path in the postscript file.
 * @fh: the file handle to print to.
 */
void draw_closepath (FILE *fh) {
  /* print the closepath statement. */
  fprintf (fh, "closepath\n");
}

/* draw_setcolor: sets a color.
 * @fh: the file handle to print to.
 * @r: the red component to set.
 * @g: the green component to set.
 * @b: the blue component to set.
 */
void draw_setcolor (FILE *fh, double r, double g, double b) {
  /* print the setrgbcolor statement. */
  fprintf (fh, "%lf %lf %lf setrgbcolor\n", r, g, b);
}

/* draw_fill: sets a color and fills.
 * @fh: the file handle to print to.
 * @r: the red component to set.
 * @g: the green component to set.
 * @b: the blue component to set.
 */
void draw_fill (FILE *fh, double r, double g, double b) {
  /* print the setrgbcolor and fill statements. */
  fprintf (fh, "%lf %lf %lf setrgbcolor\nfill\n", r, g, b);
}

/* draw_stroke: closes a path in the postscript file.
 * @fh: the file handle to print to.
 */
void draw_stroke (FILE *fh) {
  /* print the stroke statement. */
  fprintf (fh, "stroke\n");
}

/* draw_moveto: prints a moveto command to a file handle.
 * @fh: the file handle to print to.
 * @x: the x-axis position to move to.
 * @y: the y-axis position to move to.
 */
void draw_moveto (FILE *fh, double x, double y) {
  /* print the moveto statement. */
  fprintf (fh, "%lf %lf moveto\n", x, y);
}

/* draw_lineto: prints a lineto command to a file handle.
 * @fh: the file handle to print to.
 * @x: the x-axis position to line to.
 * @y: the y-axis position to line to.
 */
void draw_lineto (FILE *fh, double x, double y) {
  /* print the moveto statement. */
  fprintf (fh, "%lf %lf lineto\n", x, y);
}

/* draw_line: prints a full line segment to a file handle.
 * @fh: the file handle to print to.
 * @x1: the start x position of the line.
 * @y1: the start y position of the line.
 * @x2: the end x position of the line.
 * @y2: the end y position of the line.
 */
void draw_line (FILE *fh, double x1, double y1, double x2, double y2) {
  /* print the custom-defined line statement. */
  fprintf (fh, "%lf %lf %lf %lf l\n", x1, y1, x2, y2);
}

/* draw_text: prints a block of postscript text to a file handle.
 * @fh: the file handle to print to.
 * @size: the font size of the text.
 * @str: the text string to print.
 */
void draw_text (FILE *fh, unsigned long size, const char *str) {
  /* draw the text to the document. */
  fprintf (fh, "/Courier %lu selectfont\n(%s) show\n", size, str);
}

/* draw_state_save: saves the current path.
 * @fh: the file handle to print to.
 */
void draw_state_save (FILE *fh) {
  /* print the gsave statement. */
  fprintf (fh, "gsave\n");
}

/* draw_state_restore: restores the saved path.
 * @fh: the file handle to print to.
 */
void draw_state_restore (FILE *fh) {
  /* print the grestore statement. */
  fprintf (fh, "grestore\n");
}

/* draw_rgb: builds a unique color for each g in [0,n-1].
 * @g: the index to get a color for.
 * @n: the number of possible indices.
 * @R: pointer to the output red component.
 * @G: pointer to the output green component.
 * @B: pointer to the output blue component.
 */
void draw_rgb (unsigned long g, unsigned long n,
               double *R, double *G, double *B) {
  /* declare required variables. */
  double h, hx;

  /* calculate the desired hue of the symbol. */
  h = 5.0 * ((double) g) / ((double) n);

  /* calculate the intermediate value. */
  hx = 1.0 - fabs (fmod (h, 2.0) - 1.0);

  /* determine what the rgb values are. */
  if (h < 1.0) {
    /* set the values. */
    *R = 1.0;
    *G = hx;
    *B = 0.0;
  }
  else if (h < 2.0) {
    /* set the values. */
    *R = hx;
    *G = 1.0;
    *B = 0.0;
  }
  else if (h < 3.0) {
    /* set the values. */
    *R = 0.0;
    *G = 1.0;
    *B = hx;
  }
  else if (h < 4.0) {
    /* set the values. */
    *R = 0.0;
    *G = hx;
    *B = 1.0;
  }
  else if (h < 5.0) {
    /* set the values. */
    *R = hx;
    *G = 0.0;
    *B = 1.0;
  }
  else {
    /* set the values. */
    *R = 1.0;
    *G = 0.0;
    *B = hx;
  }
}

/* draw_symbol: draws a point in a pca/opls scores plot.
 * @fh: the file handle to print to.
 * @g: the group index, used for shape and color.
 * @n: the number of groups, used for color.
 * @x: the x-coordinate of the symbol.
 * @y: the y-coordinate of the symbol.
 */
void draw_symbol (FILE *fh, unsigned long g, unsigned long n,
                  double x, double y) {
  /* declare required variables. */
  double red, green, blue;

  /* get the group color. */
  draw_rgb (g, n, &red, &green, &blue);

  /* see which shape we should draw. */
  switch (g % 4) {
    /* square. */
    case 0:
      draw_newpath (fh);
      draw_moveto (fh, x - 2.0, y - 2.0);
      draw_lineto (fh, x - 2.0, y + 2.0);
      draw_lineto (fh, x + 2.0, y + 2.0);
      draw_lineto (fh, x + 2.0, y - 2.0);
      draw_lineto (fh, x - 2.0, y - 2.0);
      draw_closepath (fh);
      draw_state_save (fh);
      draw_fill (fh, red, green, blue);
      draw_state_restore (fh);
      draw_setcolor (fh, 0.0, 0.0, 0.0);
      draw_stroke (fh);
      break;

    /* diamond. */
    case 1:
      draw_newpath (fh);
      draw_moveto (fh, x, y - 2.8);
      draw_lineto (fh, x - 2.8, y);
      draw_lineto (fh, x, y + 2.8);
      draw_lineto (fh, x + 2.8, y);
      draw_lineto (fh, x, y - 2.8);
      draw_closepath (fh);
      draw_state_save (fh);
      draw_fill (fh, red, green, blue);
      draw_state_restore (fh);
      draw_setcolor (fh, 0.0, 0.0, 0.0);
      draw_stroke (fh);
      break;

    /* triangle. */
    case 2:
      draw_newpath (fh);
      draw_moveto (fh, x, y + 2.5);
      draw_lineto (fh, x + 2.5, y - 2.5);
      draw_lineto (fh, x - 2.5, y - 2.5);
      draw_lineto (fh, x, y + 2.5);
      draw_closepath (fh);
      draw_state_save (fh);
      draw_fill (fh, red, green, blue);
      draw_state_restore (fh);
      draw_setcolor (fh, 0.0, 0.0, 0.0);
      draw_stroke (fh);
      break;

    /* circle. */
    case 3:
      draw_newpath (fh);
      fprintf (fh, "%lf %lf 2.0 0.0 360.0 arc\n", x, y);
      draw_closepath (fh);
      draw_state_save (fh);
      draw_fill (fh, red, green, blue);
      draw_state_restore (fh);
      draw_setcolor (fh, 0.0, 0.0, 0.0);
      draw_stroke (fh);
      break;
  }
}

