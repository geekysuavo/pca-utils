
/*
 * pca-utils-math.h: header code for linear algebra.
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

#ifndef __PCA_UTILS_MATH__
#define __PCA_UTILS_MATH__

/* begin function definitions below. */

struct vector *vector_new (unsigned long n);

struct vector *vector_copy (struct vector *vi);

void vector_free (struct vector *v);

double vector_get (struct vector *v, unsigned long i);

int vector_set (struct vector *v, unsigned long i, double x);

void vector_print (struct vector *v);

int vector_sum (double da, struct vector *va,
                double db, struct vector *vb,
                struct vector *vo);

double vector_inner (struct vector *va, struct vector *vb);

double vector_norm (struct vector *v, double p);

struct matrix *matrix_new (unsigned long m, unsigned long n);

struct matrix *matrix_copy (struct matrix *Mi);

void matrix_free (struct matrix *M);

double matrix_get (struct matrix *M, unsigned long i, unsigned long j);

int matrix_set (struct matrix *M, unsigned long i, unsigned long j,
                double x);

void matrix_print (struct matrix *M);

int matrix_sum (double da, struct matrix *Ma,
                double db, struct matrix *Mb,
                struct matrix *Mo);

int matrix_vector_mul (double d,
                       struct matrix *M,
                       struct vector *v,
                       struct vector *vo);

int matrix_matrix_mul (double d,
                       struct matrix *Ma,
                       struct matrix *Mb,
                       struct matrix *Mo);

struct vector *matrix_mean (struct matrix *M);

struct matrix *matrix_cov (struct matrix *M);

int matrix_eig (struct matrix *M,
                struct vector *lambda,
                struct matrix *Q);

double matrix_eig_det (struct vector *lambda);

int matrix_chol (struct matrix *M, struct matrix *L);

int matrix_chol_inv (struct matrix *L, struct matrix *Minv);

double matrix_chol_det (struct matrix *L);

int matrix_resize (struct matrix *M, unsigned long m, unsigned long n);

int matrix_del_row (struct matrix *M, unsigned long idel);

int matrix_del_column (struct matrix *M, unsigned long jdel);

void matrix_rotx (struct matrix *R, double theta);

void matrix_roty (struct matrix *R, double theta);

void matrix_rotz (struct matrix *R, double theta);

struct matrix *matrix_rot (double alpha, double beta, double gamma);

#endif

