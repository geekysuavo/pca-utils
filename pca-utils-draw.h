
/*
 * pca-utils-draw.h: header code for postscript drawing.
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

#ifndef __PCA_UTILS_DRAW__
#define __PCA_UTILS_DRAW__

/* begin function definitions below. */

FILE *draw_document_open (const char *fname,
                          unsigned long w,
                          unsigned long h);

void draw_document_close (FILE *fh);

void draw_newpath (FILE *fh);

void draw_closepath (FILE *fh);

void draw_setcolor (FILE *fh, double r, double g, double b);

void draw_fill (FILE *fh, double r, double g, double b);

void draw_stroke (FILE *fh);

void draw_moveto (FILE *fh, double x, double y);

void draw_lineto (FILE *fh, double x, double y);

void draw_line (FILE *fh, double x1, double y1, double x2, double y2);

void draw_text (FILE *fh, unsigned long size, const char *str);

void draw_state_save (FILE *fh);

void draw_state_restore (FILE *fh);

void draw_rgb (unsigned long g, unsigned long n,
               double *R, double *G, double *B);

void draw_symbol (FILE *fh, unsigned long g, unsigned long n,
                  double x, double y);

#endif

