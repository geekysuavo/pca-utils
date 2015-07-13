
/*
 * pca-utils-math.c: source code for linear algebra.
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

/* vector_new: creates a pointer to a newly allocated vector struct.
 * @n: the desired dimensionality of the new vector.
 */
struct vector *vector_new (unsigned long n) {
  /* declare required variables. */
  struct vector *v;
  unsigned long i;

  /* ensure the requested dimensionality makes sense. */
  if (!n) {
    /* output an error message and return nothing. */
    ferror ("nonzero vector size required");
    return NULL;
  }

  /* allocate a pointer to a vector structure. */
  v = (struct vector *) malloc (sizeof (struct vector));

  /* ensure the vector pointer was allocated. */
  if (!v) {
    /* outout an error message and return nothing. */
    ferror ("failed to allocate pca_vector pointer");
    return NULL;
  }

  /* store the dimensionality in the structure. */
  v->n = n;

  /* allocate memory fot the array of values. */
  v->d = (double *) malloc (sizeof (double) * v->n);

  /* ensure the array memory was allocated. */
  if (!v->d) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate data array");
    return NULL;
  }

  /* loop through the values, zeroing them out. */
  for (i = 0; i < v->n; i++)
    v->d[i] = 0.0;

  /* return the newly created and initialized vector struct pointer. */
  return v;
}

/* vector_copy: copies a vector.
 * @vi: the input vector.
 */
struct vector *vector_copy (struct vector *vi) {
  /* declare required variables. */
  struct vector *vo;
  unsigned long i;

  /* make sure the input vector is defined. */
  if (!vi) {
    /* output an error message and return nothing. */
    ferror ("input vector to copy is undefined");
    return NULL;
  }

  /* make a new vector. */
  vo = vector_new (vi->n);

  /* make sure we allocated a new vector. */
  if (!vo) {
    /* output an error message and return nothng. */
    ferror ("failed to create copy vector of size %lu", vi->n);
    return NULL;
  }

  /* copy the values. */
  for (i = 0; i < vi->n; i++)
    vo->d[i] = vi->d[i];

  /* return the new vector. */
  return vo;
}

/* vector_free: frees a vector structure.
 * @v: the vector to free.
 */
void vector_free (struct vector *v) {
  /* don't free a null structure. */
  if (!v) return;

  /* free the values array. */
  free (v->d);
  v->d = NULL;

  /* set the number of values to zero. */
  v->n = 0;

  /* free the structure pointer. */
  free (v);
  v = NULL;
}

/* vector_get: get the value of a vector element.
 * @v: the vector to operate on.
 * @i: the element index.
 */
double vector_get (struct vector *v, unsigned long i) {
  /* make sure the vector is defined. */
  if (!v) {
    /* output an error message and return failure. */
    ferror ("vector undefined");
    return 0.0;
  }

  /* make sure the element index is in bounds. */
  if (i >= v->n) {
    /* output an error message and return failure. */
    ferror ("vector index %lu out of bounds", i);
    return 0.0;
  }

  /* return the vector element. */
  return v->d[i];
}

/* vector_set: set the value of a vector element.
 * @v: the vector to operate on.
 * @i: the vector element index to set.
 * @x: the new value of the vector element.
 */
int vector_set (struct vector *v, unsigned long i, double x) {
  /* make sure the vector is defined. */
  if (!v) {
    /* output an error message and return failure. */
    ferror ("vector undefined");
    return 0;
  }

  /* make sure the element index is in bounds. */
  if (i >= v->n) {
    /* output an error message and return failure. */
    ferror ("vector index %lu out of bounds", i);
    return 0;
  }

  /* set the vector element. */
  v->d[i] = x;

  /* return success. */
  return 1;
}

/* vector_print: prints a vector to standard output.
 * @v: the vector to output.
 */
void vector_print (struct vector *v) {
  /* declare required variables. */
  unsigned long i;

  /* loop through the elements. */
  for (i = 0; i < v->n; i++)
    fprintf (stdout, "  %16le\n", v->d[i]);
}

/* vector_sum: sums two scaled vectors.
 * @da: the first scale factor.
 * @va: the first input vector.
 * @db: the second scale factor.
 * @vb: the second input vector.
 * @vo: the output vector.
 */
int vector_sum (double da, struct vector *va,
                double db, struct vector *vb,
                struct vector *vo) {
  /* declare required variables. */
  unsigned long i;

  /* make sure the two vectors are defined. */
  if (!va || !vb || !vo) {
    /* output an error message and return failure. */
    ferror ("input vectors to sum not defined");
    return 0;
  }

  /* make sure the two vectors are equal size. */
  if (va->n != vb->n) {
    /* output an error message and return failure. */
    ferror ("input vector to sum size mismatch");
    return 0;
  }

  /* make sure the output vector is the correct size. */
  if (vo->n != va->n) {
    /* output an error message and return failure. */
    ferror ("output vector of sum size mismatch");
    return 0;
  }

  /* loop through the elements. */
  for (i = 0; i < va->n; i++)
    vo->d[i] = da * va->d[i] + db * vb->d[i];

  /* return success. */
  return 1;
}

/* vector_inner: calculate the inner (dot) product of two vectors.
 * @va: the first input vector.
 * @vb: the second input vector.
 */
double vector_inner (struct vector *va, struct vector *vb) {
  /* declare required variables. */
  unsigned long i;
  double d;

  /* make sure the vectors are defined. */
  if (!va || !vb) {
    /* output an error message and return nan. */
    ferror ("inner product input vectors undefined");
    return NAN;
  }

  /* make sure the vector sizes match. */
  if (va->n != vb->n) {
    /* output an error message and return nan. */
    ferror ("inner product input vector size mismatch");
    return NAN;
  }

  /* loop through the elements, calculating the product. */
  for (d = 0.0, i = 0; i < va->n; i++)
    d += va->d[i] * vb->d[i];

  /* return the product. */
  return d;
}

/* vector_norm: calculate the p-norm of a vector.
 * @v: the vector of which to calculate the norm.
 * @p: the type of norm to calculate.
 */
double vector_norm (struct vector *v, double p) {
  /* declare required variables. */
  unsigned long i;
  double norm;

  /* make sure the input vector is defined. */
  if (!v) {
    /* output an error and return nan. */
    ferror ("input vector to norm not defined");
    return NAN;
  }

  /* make sure the norm is greated than or equal to one. */
  if (p < 1.0) {
    /* output an error message and return nan. */
    ferror ("%.2lf-norm of vector is invalid", p);
    return NAN;
  }

  /* calculate the radicand of the norm. */
  for (norm = 0.0, i = 0; i < v->n; i++)
    norm += pow (v->d[i], p);

  /* return the norm. */
  return pow (norm, 1.0 / p);
}

/* matrix_new: creates a pointer to a newly allocated matrix struct.
 * @m: the row-wise dimensionality of the matrix.
 * @n: the column-wise dimensionality of the matrix.
 */
struct matrix *matrix_new (unsigned long m, unsigned long n) {
  /* declare required variables. */
  struct matrix *M;
  unsigned long i, j;

  /* ensure the requested dimensionalities make sense. */
  if (!m || !n) {
    /* output an error message and return nothing. */
    ferror ("nonzero matrix sizes required (%lux%lu invalid)", m, n);
    return NULL;
  }

  /* allocate a pointer to the matrix struct. */
  M = (struct matrix *) malloc (sizeof (struct matrix));

  /* ensure the matrix struct pointer was allocated. */
  if (!M) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate pca_matrix pointer");
    return NULL;
  }

  /* store the dimensionalities in the structure. */
  M->m = m;
  M->n = n;

  /* allocate memory for the table of values. */
  M->d = (double **) malloc (sizeof (double *) * M->m);

  /* ensure the matrix row memory was allocated. */
  if (!M->d) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate data array");
    return NULL;
  }

  /* loop through the matrix rows. */
  for (i = 0; i < M->m; i++) {
    /* allocate the column values in the row. */
    M->d[i] = (double *) malloc (sizeof (double) * M->n);

    /* ensure the column array was allocated. */
    if (!M->d[i]) {
      /* output an error message and return nothing. */
      ferror ("failed to allocate data array row '%lu'", i);
      return NULL;
    }

    /* loop through the column values, zeroing them out. */
    for (j = 0; j < M->n; j++)
      M->d[i][j] = 0.0;
  }

  /* return the newly created and initialized matrix. */
  return M;
}

/* matrix_copy: copies a matrix structure.
 * @Mi: the input matrix.
 */
struct matrix *matrix_copy (struct matrix *Mi) {
  /* declare required variables. */
  unsigned long i, j;
  struct matrix *Mo;

  /* make sure the input matrix is defined. */
  if (!Mi) {
    /* output an error message and return nothing. */
    ferror ("input matrix to copy is undefined");
    return NULL;
  }

  /* make a new matrix. */
  Mo = matrix_new (Mi->m, Mi->n);

  /* make sure we made the matrix. */
  if (!Mo) {
    /* output an error message and return nothing. */
    ferror ("failed to create copy matrix");
    return NULL;
  }

  /* copy the matrix values. */
  for (i = 0; i < Mi->m; i++)
    for (j = 0; j < Mi->n; j++)
      Mo->d[i][j] = Mi->d[i][j];

  /* return the copied matrix. */
  return Mo;
}

/* matrix_free: frees a matrix structure.
 * @M: the matrix to free.
 */
void matrix_free (struct matrix *M) {
  /* declare required variables. */
  unsigned long i;

  /* don't free a null structure. */
  if (!M) return;

  /* loop through the rows. */
  for (i = 0; i < M->m; i++) {
    /* free the currently indexes column array. */
    free (M->d[i]);
    M->d[i] = NULL;
  }

  /* free the values array. */
  free (M->d);
  M->d = NULL;

  /* set the number of rows and columns to zero. */
  M->m = 0;
  M->n = 0;

  /* free the structure pointer. */
  free (M);
  M = NULL;
}

/* matrix_get: get the value of a matrix element.
 * @M: the matrix to operate on.
 * @i: the row index to get.
 * @j: the column index to get.
 */
double matrix_get (struct matrix *M,
                   unsigned long i,
                   unsigned long j) {
  /* ensure the matrix is defined. */
  if (!M) {
    /* output an error message and return failure. */
    ferror ("matrix undefined");
    return 0.0;
  }

  /* ensure the indices are in bounds. */
  if (i >= M->m || j >= M->n) {
    /* output an error message and return failure. */
    ferror ("matrix indices %lu,%lu out of bounds", i, j);
    return 0.0;
  }

  /* return the matrix element. */
  return M->d[i][j];
}

/* matrix_set: sets the value of a matrix element.
 * @M: the matrix to operate on.
 * @i: the row index to set.
 * @j: the column index to set.
 * @x: the value to set.
 */
int matrix_set (struct matrix *M,
                unsigned long i,
                unsigned long j,
                double x) {
  /* ensure the matrix is defined. */
  if (!M) {
    /* output an error message and return failure. */
    ferror ("matrix undefined");
    return 0;
  }

  /* ensure the indices are in bounds. */
  if (i >= M->m || j >= M->n) {
    /* output an error message and return failure. */
    ferror ("matrix indices %lu,%lu out of bounds", i, j);
    return 0;
  }

  /* set the matrix element. */
  M->d[i][j] = x;

  /* return success. */
  return 1;
}

/* matrix_print: prints a matrix to standard output.
 * @M: the matrix to output.
 */
void matrix_print (struct matrix *M) {
  /* declare required variables. */
  unsigned long i, j;

  /* loop through the matrix rows. */
  for (i = 0; i < M->m; i++) {
    /* and through the matrix columns. */
    for (j = 0; j < M->n; j++)
      fprintf (stdout, " %16le ", M->d[i][j]);

    /* print a newline. */
    fprintf (stdout, "\n");
  }
}

/* matrix_sum: calculate the sum of two scaled matrices.
 * @da: the first scale factor.
 * @Ma: the first input matrix.
 * @db: the second scale factor.
 * @Mb: the second input matrix.
 * @Mo: the output matrix.
 */
int matrix_sum (double da, struct matrix *Ma,
                double db, struct matrix *Mb,
                struct matrix *Mo) {
  /* declare required variables. */
  unsigned long i, j;

  /* make sure the two matrices are defined. */
  if (!Ma || !Mb || !Mo) {
    /* output an error message and return failure. */
    ferror ("input matrices to sum not defined");
    return 0;
  }

  /* make sure the two matrices are equal size. */
  if (Ma->m != Mb->m || Ma->n != Mb->n) {
    /* output an error message and return failure. */
    ferror ("input matrix to sum size mismatch");
    return 0;
  }

  /* make sure the output matrix is the correct size. */
  if (Mo->m != Ma->m || Mo->n != Ma->n) {
    /* output an error message and return failure. */
    ferror ("output matrix of sum size mismatch");
    return 0;
  }

  /* sum the values in the matrices. */
  for (i = 0; i < Ma->m; i++)
    for (j = 0; j < Ma->n; j++)
      Mo->d[i][j] = da * Ma->d[i][j] + db * Mb->d[i][j];

  /* return success. */
  return 1;
}

/* matrix_vector_mul: right-multiplies a matrix by a column vector.
 * @d: the scale value.
 * @M: the input matrix.
 * @v: the input vector.
 * @vo: the output vector.
 */
int matrix_vector_mul (double d,
                       struct matrix *M,
                       struct vector *v,
                       struct vector *vo) {
  /* declare required variables. */
  unsigned long i, j;

  /* make sure the input structures are defined. */
  if (!M || !v || !vo) {
    /* output an error message and return failure. */
    ferror ("input structures to matrix-vector multiply undefined");
    return 0;
  }

  /* make sure the sizes match. */
  if (M->m != vo->n || M->n != v->n) {
    /* output an error message and return failure. */
    ferror ("matrix-vector input structure size mismatch");
    return 0;
  }

  /* calculate the matrix-vector product. */
  for (i = 0; i < M->m; i++)
    for (vo->d[i] = 0.0, j = 0; j < M->n; j++)
      vo->d[i] += d * M->d[i][j] * v->d[j];

  /* return success. */
  return 1;
}

/* matrix_matrix_mul: multiplies two matrices.
 * @d: the scale value.
 * @Ma: the first matrix.
 * @Mb: the second matrix.
 * @Mo: the output matrix.
 */
int matrix_matrix_mul (double d,
                       struct matrix *Ma,
                       struct matrix *Mb,
                       struct matrix *Mo) {
  /* declare required variables. */
  unsigned long i, j, k;

  /* make sure the input structures are defined. */
  if (!Ma || !Mb || !Mo) {
    /* output an error message and return failure. */
    ferror ("input structures to matrix-matrix multiply undefined");
    return 0;
  }

  /* make sure the sizes match. */
  if (Ma->n != Mb->m || Ma->m != Mo->m || Mb->n != Mo->n) {
    /* output an error message and return failure. */
    ferror ("matrix-matrix input structure size mismatch");
    return 0;
  }

  /* calculate the matrix-matrix product. */
  for (i = 0; i < Ma->m; i++)
    for (j = 0; j < Mb->n; j++)
      for (Mo->d[i][j] = 0.0, k = 0; k < Ma->n; k++)
        Mo->d[i][j] += d * Ma->d[i][k] * Mb->d[k][j];

  /* return success. */
  return 1;
}

/* matrix_mean: calculate the mean vector of a data matrix.
 * @M: the data matrix to operate on.
 */
struct vector *matrix_mean (struct matrix *M) {
  /* declare required variables. */
  unsigned long i, j;
  struct vector *mu;

  /* allocate a mean vector. */
  mu = vector_new (M->n);

  /* ensure we got a vector. */
  if (!mu) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate mean vector");
    return NULL;
  }

  /* loop through the dimensions. */
  for (j = 0; j < M->n; j++) {
    /* initialize the mean value. */
    mu->d[j] = 0.0;

    /* loop through the rows of the data matrix. */
    for (i = 0; i < M->m; i++) {
      /* sum into the mean vector. */
      mu->d[j] += M->d[i][j];
    }

    /* normalize the value to a mean. */
    mu->d[j] /= (double) M->m;
  }

  /* return the mean vector. */
  return mu;
}

/* matrix_cov: calculate the covariance matrix of a data matrix.
 * @M: the data matrix to operate on.
 */
struct matrix *matrix_cov (struct matrix *M) {
  /* declare required variables. */
  unsigned long i, j, k;
  struct vector *mu;
  struct matrix *C;

  /* ensure the matrix contains sufficient rows to calculate covariance. */
  if (M->m < 2) {
    /* output an error message and return nothing. */
    ferror ("matrix has too few rows, cannot calculate covariance");
    return NULL;
  }

  /* calculate the matrix column means. */
  mu = matrix_mean (M);

  /* ensure we got a mean vector. */
  if (!mu) {
    /* output an error message and return nothing. */
    return NULL;
  }

  /* allocate a covariance matrix. */
  C = matrix_new (M->n, M->n);

  /* ensure we allocated a matrix. */
  if (!C) {
    /* output an error message and return nothing. */
    return NULL;
  }

  /* loop through the dimensions once. */
  for (j = 0; j < M->n; j++) {
    /* and loop through the dimensions again. */
    for (k = j; k < M->n; k++) {
      /* initialize the matrix element. */
      C->d[j][k] = 0.0;

      /* loop through the rows in the current dimension pair. */
      for (i = 0; i < M->m; i++) {
        /* sum up into the matrix element. */
        C->d[j][k] += (M->d[i][j] - mu->d[j]) * (M->d[i][k] - mu->d[k]);
      }

      /* normalize the matrix element to a variance. */
      C->d[j][k] /= (double) (M->m - 1);

      /* are we off-diagonal? */
      if (j != k) C->d[k][j] = C->d[j][k];
    }
  }

  /* free the mean vector. */
  vector_free (mu);

  /* return the covariance matrix. */
  return C;
}
/* matrix_eig_sort: sorts the eigenvalues and eigenvectors extracted from
 * a matrix into descending order by eigenvalue. you have no need to run
 * this in isolation; it's a helper function to matrix_eig() and assumes
 * the checks on memory sanity there have already been performed.
 * @lambda: the input/output eigenvalue vector.
 * @Q: the input/output eigenvector matrix.
 */
void matrix_eig_sort (struct vector *lambda, struct matrix *Q) {
  /* declare required variables. */
  unsigned long n, i, j;
  double tmp;
  int done;

  /* get the size of the vector. */
  n = lambda->n;
  done = 0;

  /* loop through the vector. */
  while (!done) {
    /* we can be optimistic, right? */
    done = 1;

    /* look for swappable indices. */
    for (i = 1; i < n; i++) {
      /* are we swappable? */
      if (lambda->d[i - 1] < lambda->d[i]) {
        /* swap the values. */
        tmp = lambda->d[i - 1];
        lambda->d[i - 1] = lambda->d[i];
        lambda->d[i] = tmp;

        /* swap the columns. */
        for (j = 0; j < n; j++) {
          /* swap the values. */
          tmp = Q->d[j][i - 1];
          Q->d[j][i - 1] = Q->d[j][i];
          Q->d[j][i] = tmp;
        }

        /* not done yet! */
        done = 0;
      }
    }
  }
}

/* matrix_eig: eigendecompose a (modestly sized) matrix using jacobi
 * rotations.
 * @M: the matrix to eigendecompose.
 * @lambda: the vector of output eigenvalues.
 * @Q: the matrix of output eigenvectors.
 */
int matrix_eig (struct matrix *M,
                struct vector *lambda,
                struct matrix *Q) {
  /* declare required variables. */
  double sum, thresh, g, h, rg, rh, t, theta, tau, c, s;
  unsigned long n, i, j, ip, iq;
  struct vector *b, *z;
  struct matrix *A;

  /* ensure the input matrix is defined and square. */
  if (!M || M->m != M->n) {
    /* output an error message and return failure. */
    ferror ("input matrix undefined or not square");
    return 0;
  }

  /* ensure the output vector is defined and matches in size. */
  if (!lambda || lambda->n != M->m) {
    /* output an error message and return failure. */
    ferror ("eigenvalue vector undefined or of invalid size");
    return 0;
  }

  /* ensure the output matrix is defined and matches in size. */
  if (!Q || Q->m != Q->n || Q->m != M->m) {
    /* output an error message and return failure. */
    ferror ("eigenvector matrix undefined or of invalid size");
    return 0;
  }

  /* store the dimensionality of the problem in an easy variable. */
  n = M->m;

  /* make a copy of the input matrix. */
  A = matrix_copy (M);

  /* make sure we made a copy. */
  if (!A) {
    /* output an error message and return failure. */
    ferror ("failed to copy input matrix to eigendecomposition");
    return 0;
  }

  /* make two temporary vectors. */
  b = vector_new (n);
  z = vector_new (n);

  /* make sure we allocated the temporary vectors. */
  if (!b || !z) {
    /* output an error message and return failure. */
    ferror ("failed to allocate temporary vectors");
    return 0;
  }

  /* initialize the eigenvector matrix to the identity. */
  for (ip = 0; ip < n; ip++)
    for (iq = 0; iq < n; iq++)
      Q->d[ip][iq] = (ip == iq ? 1.0 : 0.0);

  /* initialize vectors to the input matrix diagonal. */
  for (ip = 0; ip < n; ip++) {
    /* set the values. */
    b->d[ip] = lambda->d[ip] = A->d[ip][ip];
    z->d[ip] = 0.0;
  }

  /* iterate over the matrix a maximum number of times. */
  for (i = 1; i <= 50; i++) {
    /* sum the off-diagonal elements. */
    for (sum = 0.0, ip = 0; ip < n - 1; ip++)
      for (iq = ip + 1; iq < n; iq++)
        sum += fabs (A->d[ip][iq]);

    /* see if we've converged. */
    if (sum == 0.0) {
      /* sort the eigenvalues and eigenvectors. */
      matrix_eig_sort (lambda, Q);

      /* free the temporary structures. */
      matrix_free (A);
      vector_free (b);
      vector_free (z);

      /* return success. */
      return 1;
    }

    /* calculate the threshold. */
    if (i < 4)
      thresh = 0.2 * sum / ((double) n * (double) n);
    else
      thresh = 0.0;

    /* loop again through the off-diagonals. */
    for (ip = 0; ip < n - 1; ip++) {
      for (iq = ip + 1; iq < n; iq++) {
        /* get the value of the current off-diagonal element. */
        g = 100.0 * fabs (A->d[ip][iq]);

        /* skip rotation in the case of late small element values. */
        if (i > 4 &&
            g <= DBL_EPSILON * fabs (lambda->d[ip]) &&
            g <= DBL_EPSILON * fabs (lambda->d[iq])) {
          /* skip the rotation. */
          A->d[ip][iq] = 0.0;
        }
        else if (fabs (A->d[ip][iq]) > thresh) {
          /* get the difference in the eigenvalues. */
          h = lambda->d[iq] - lambda->d[ip];

          /* calculate the t parameter. */
          if (g <= DBL_EPSILON * fabs (h)) {
            /* t = 1 / (2*theta) */
            t = A->d[ip][iq] / h;
          }
          else {
            /* otherwise... */
            theta = 0.5 * h / A->d[ip][iq];
            t = 1.0 / (fabs (theta) + sqrt (1.0 + theta * theta));
            if (theta < 0.0) t = -t;
          }

          /* calculate values for the rotations. */
          c = 1.0 / sqrt (1 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * A->d[ip][iq];

          /* perform the rotations. */
          z->d[ip] -= h;
          z->d[iq] += h;
          lambda->d[ip] -= h;
          lambda->d[iq] += h;

          /* zero the matrix element. */
          A->d[ip][iq] = 0.0;

          /* case of rotations 0 <= j < p */
          for (j = 0; j < ip; j++) {
            /* perform the rotation operation. */
            rg = A->d[j][ip];
            rh = A->d[j][iq];
            A->d[j][ip] = rg - s * (rh + rg * tau);
            A->d[j][iq] = rh + s * (rg - rh * tau);
          }

          /* case of rotations p < j < q */
          for (j = ip + 1; j < iq; j++) {
            /* perform the rotation operation. */
            rg = A->d[ip][j];
            rh = A->d[j][iq];
            A->d[ip][j] = rg - s * (rh + rg * tau);
            A->d[j][iq] = rh + s * (rg - rh * tau);
          }

          /* case of rotations q < j < n */
          for (j = iq + 1; j < n; j++) {
            /* perform the rotation operation. */
            rg = A->d[ip][j];
            rh = A->d[iq][j];
            A->d[ip][j] = rg - s * (rh + rg * tau);
            A->d[iq][j] = rh + s * (rg - rh * tau);
          }

          /* final case. */
          for (j = 0; j < n; j++) {
            /* perform the rotation operation. */
            rg = Q->d[j][ip];
            rh = Q->d[j][iq];
            Q->d[j][ip] = rg - s * (rh + rg * tau);
            Q->d[j][iq] = rh + s * (rg - rh * tau);
          }
        }
      }
    }

    /* update the eigenvalue vector. */
    for (ip = 0; ip < n; ip++) {
      /* set the values. */
      b->d[ip] += z->d[ip];
      lambda->d[ip] = b->d[ip];
      z->d[ip] = 0.0;
    }
  }

  /* free the temporary structures. */
  matrix_free (A);
  vector_free (b);
  vector_free (z);

  /* output an error message. */
  ferror ("exceeded maximum iteration count in eigendecomposition");

  /* return success. */
  return 0;
}

/* matrix_eig_det: calculate the determinant of an eigendecomposed matrix.
 * @lambda: the vector of eigenvalues of the matrix to analyze.
 */
double matrix_eig_det (struct vector *lambda) {
  /* declare required variables. */
  double det = 1.0;
  unsigned long i;

  /* loop through the eigenvalues. */
  for (i = 0; i < lambda->n; i++) {
    /* multiply by the current eigenvalue. */
    det *= lambda->d[i];
  }

  /* return the determinant. */
  return det;
}

/* matrix_chol: cholesky decompose a matrix.
 * @M: the input matrix.
 * @L: the lower triangular output matrix.
 */
int matrix_chol (struct matrix *M, struct matrix *L) {
  /* declare required variables. */
  signed long i, j, k;
  double tmp;

  /* ensure the input matrix is defined and square. */
  if (!M || M->m != M->n) {
    /* output an error message and return failure. */
    ferror ("input matrix undefined or not square");
    return 0;
  }

  /* ensure the output matrix is defined and matches in size. */
  if (!L || L->m != L->n || L->m != M->m) {
    /* output an error message and return failure. */
    ferror ("lower triangular matrix undefined or of invalid size");
    return 0;
  }

  /* copy the input matrix into the output matrix. */
  for (i = 0; i < M->m; i++)
    for (j = 0; j < M->n; j++)
      L->d[i][j] = M->d[i][j];

  /* loop through the matrix rows. */
  for (i = 0; i < L->m; i++) {
    /* and through the upper triangle. */
    for (j = i; j < L->n; j++) {
      /* sum up the values in the column. */
      for (tmp = L->d[i][j], k = i - 1; k >= 0; k--)
        tmp -= L->d[i][k] * L->d[j][k];

      /* act based on sitting on the diagonal. */
      if (i == j) {
        /* check for a zero pivot. */
        if (tmp <= 0.0) {
          /* return failure. */
          return 0;
        }

        /* store the result. */
        L->d[i][i] = sqrt (tmp);
      }
      else {
        /* store the result. */
        L->d[j][i] = tmp / L->d[i][i];
      }
    }
  }

  /* zero out the upper triangle. */
  for (i = 0; i < L->m; i++)
    for (j = 0; j < i; j++)
      L->d[j][i] = 0.0;

  /* return success. */
  return 1;
}

/* matrix_chol_inv: computes the inverse of a cholesky-decomposed matrix.
 * @L: the lower triangular decomposition of M from matrix_chol().
 * @Minv: the output inverse matrix of M.
 */
int matrix_chol_inv (struct matrix *L, struct matrix *Minv) {
  /* declare required variables. */
  signed long i, j, k;
  double tmp;

  /* ensure the input matrix is defined and square. */
  if (!L || L->m != L->n) {
    /* output an error message and return failure. */
    ferror ("lower triangular matrix undefined or not square");
    return 0;
  }

  /* ensure the output matrix is defined and matches in size. */
  if (!Minv || Minv->m != Minv->n || Minv->m != L->m) {
    /* output an error message and return failure. */
    ferror ("inverse matrix undefined or of invalid size");
    return 0;
  }

  /* loop through the rows of the decomposition. */
  for (i = 0; i < L->m; i++) {
    /* and through the lower triangle. */
    for (j = 0; j <= i; j++) {
      /* initialize the temporary parameter. */
      tmp = (i == j ? 1.0 : 0.0);

      /* sum over the rows. */
      for (k = i - 1; k >= j; k--)
        tmp -= L->d[i][k] * Minv->d[j][k];

      /* store the value. */
      Minv->d[j][i] = tmp / L->d[i][i];
    }
  }

  /* loop backwards through the rows. */
  for (i = L->m - 1; i >= 0; i--) {
    /* and through the lower triangle. */
    for (j = 0; j <= i; j++) {
      /* initialize the temporary parameter. */
      tmp = (i < j ? 0.0 : Minv->d[j][i]);

      /* sum over the columns. */
      for (k = i + 1; k < L->m; k++)
        tmp -= L->d[k][i] * Minv->d[j][k];

      /* store the value. */
      Minv->d[i][j] = Minv->d[j][i] = tmp / L->d[i][i];
    }
  }

  /* return success. */
  return 1;
}

/* matrix_chol_det: calculate the determinant of a cholesky factored matrix.
 * @L: the cholesky factorization of the original matrix M.
 */
double matrix_chol_det (struct matrix *L) {
  /* declare required variables. */
  unsigned long i;
  double d;

  /* make sure the matrix is defined. */
  if (!L) {
    /* output an error message and return nan. */
    ferror ("input matrix to cholesky determinant is undefined");
    return NAN;
  }

  /* calculate the determinant. */
  for (d = 1.0, i = 0; i < L->m; i++)
    d *= L->d[i][i];

  /* return the determinant. */
  return d;
}

/* matrix_resize: resizes a matrix.
 * @M: the matrix to resize.
 * @m: the new number of rows.
 * @n: the new number of columns.
 */
int matrix_resize (struct matrix *M, unsigned long m, unsigned long n) {
  /* declare required variables. */
  unsigned long i;

  /* reallocate the matrix row array. */
  M->d = (double **) realloc (M->d, sizeof (double *) * m);

  /* make sure the reallocation was successful. */
  if (!M->d) {
    /* output an error message and return failure. */
    ferror ("failed to reallocate matrix to %lu rows", m);
    return 0;
  }

  /* loop through the matrix rows. */
  for (i = 0; i < m; i++) {
    /* does the row need fresh allocation? */
    if (i >= M->m) {
      /* yes. allocate the matrix row. */
      M->d[i] = (double *) calloc (n, sizeof (double));
    }
    else {
      /* no. reallocate the matrix row. */
      M->d[i] = (double *) realloc (M->d[i], sizeof (double) * n);
    }

    /* make sure the row reallocation was successful. */
    if (!M->d[i]) {
      /* output an error message and return failure. */
      ferror ("failed to reallocate matrix row %lu to %lu columns", i, n);
      return 0;
    }
  }

  /* store the new sizes. */
  M->m = m;
  M->n = n;

  /* return success. */
  return 1;
}

/* matrix_del_row: deletes a row of data from a matrix.
 * @M: the matrix to modify.
 * @idel: the row index to remove.
 */
int matrix_del_row (struct matrix *M, unsigned long idel) {
  /* declare required variables. */
  unsigned long i;

  /* free the row to remove. */
  free (M->d[idel]);
  M->d[idel] = NULL;

  /* loop through the matrix rows. */
  for (i = idel + 1; i < M->m; i++) {
    /* move the row up. */
    M->d[i - 1] = M->d[i];
  }

  /* return success if the matrix resizes. */
  return matrix_resize (M, M->m - 1, M->n);
}

/* matrix_del_column: deletes a column of data from a matrix.
 * @M: the matrix to modify.
 * @jdel: the column index to remove.
 */
int matrix_del_column (struct matrix *M, unsigned long jdel) {
  /* declare required variables. */
  unsigned long i, j;

  /* loop through the matrix rows. */
  for (i = 0; i < M->m; i++) {
    /* loop through the row values. */
    for (j = jdel + 1; j < M->n; j++) {
      /* move the value left. */
      M->d[i][j - 1] = M->d[i][j];
    }
  }

  /* return success if the matrix resizes. */
  return matrix_resize (M, M->m, M->n - 1);
}

/* matrix_rotx: create an x-rotation matrix.
 * @R: the matrix to set values in.
 * @theta: the rotation angle.
 */
void matrix_rotx (struct matrix *R, double theta) {
  /* set the values. */
  R->d[0][0] = 1.0;
  R->d[0][1] = 0.0;
  R->d[0][2] = 0.0;
  R->d[1][0] = 0.0;
  R->d[1][1] = cos (theta);
  R->d[1][2] = -sin (theta);
  R->d[2][0] = 0.0;
  R->d[2][1] = sin (theta);
  R->d[2][2] = cos (theta);
}

/* matrix_roty: create a y-rotation matrix.
 * @R: the matrix to set values in.
 * @theta: the rotation angle.
 */
void matrix_roty (struct matrix *R, double theta) {
  /* set the values. */
  R->d[0][0] = cos (theta);
  R->d[0][1] = 0.0;
  R->d[0][2] = sin (theta);
  R->d[1][0] = 0.0;
  R->d[1][1] = 1.0;
  R->d[1][2] = 0.0;
  R->d[2][0] = -sin (theta);
  R->d[2][1] = 0.0;
  R->d[2][2] = cos (theta);
}

/* matrix_rotz: create a z-rotation matrix.
 * @R: the matrix to set values in.
 * @theta: the rotation angle.
 */
void matrix_rotz (struct matrix *R, double theta) {
  /* set the values. */
  R->d[0][0] = cos (theta);
  R->d[0][1] = -sin (theta);
  R->d[0][2] = 0.0;
  R->d[1][0] = sin (theta);
  R->d[1][1] = cos (theta);
  R->d[1][2] = 0.0;
  R->d[2][0] = 0.0;
  R->d[2][1] = 0.0;
  R->d[2][2] = 1.0;
}

/* matrix_rot: create a rotation matrix from the three proper
 * euler angles: alpha, beta and gamma, according to the
 * following: Rzz (gamma) * Rx(beta) * Rz(alpha)
 * @alpha: the first euler angle.
 * @beta: the second euler angle.
 * @gamma: the third euler angle.
 */
struct matrix *matrix_rot (double alpha, double beta, double gamma) {
  /* declare required variables. */
  struct matrix *R, *R1, *R2, *R3, *Ri;

  /* allocate the matrices. */
  R = matrix_new (3, 3);
  R1 = matrix_new (3, 3);
  R2 = matrix_new (3, 3);
  R3 = matrix_new (3, 3);
  Ri = matrix_new (3, 3);

  /* ensure the matrices were allocated. */
  if (!R || !R1 || !R2 || !R3 || !Ri) {
    /* output an error and return nothing. */
    ferror ("failed to allocate euler rotation matrices");
    return NULL;
  }

  /* create the individual rotation matrices. */
  matrix_rotz (R1, alpha);
  matrix_rotx (R2, beta);
  matrix_rotz (R3, gamma);

  /* multiply the first two matrices together. */
  matrix_matrix_mul (1.0, R3, R2, Ri);

  /* multiply the intermediate matrix with the third. */
  matrix_matrix_mul (1.0, Ri, R1, R);

  /* free the individual rotation matrices. */
  matrix_free (R1);
  matrix_free (R2);
  matrix_free (R3);
  matrix_free (Ri);

  /* return the final matrix. */
  return R;
}

