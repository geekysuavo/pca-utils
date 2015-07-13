
/*
 * pca-utils.c: source code for pca-utils core functions.
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

/* use the tab character for delimiter and a buffer size of 2k chars. */
#define DELIM "\t"
#define N_BUF 2048

/* ferror_fn: function used by ferror() to output descriptive error messages.
 * @f: the filename string from which the error emanated.
 * @l: the line of the file from which the error emanated.
 * @format: the printf-style format string defining the message.
 * @...: the arguments required based on @format.
 */
void ferror_fn (const char *f, const int l, const char *format, ...) {
  /* declare required variables. */
  va_list vl;
  char *str;

  /* begin reading in the variable argument list. */
  va_start (vl, format);

  /* allocate memory for the output string. */
  str = (char *) malloc (sizeof (char) * N_BUF);

  /* ensure the string memory was allocated. */
  if (!str) {
    /* print a failure message and exit the function. */
    fprintf (stderr, "wow, sad... ferror_fn failed to print a failure!\n");
    return;
  }

  /* write out the variable arguments list to the output string. */
  vsnprintf (str, N_BUF, format, vl);

  /* print the output string and the extra information to stderr. */
  fprintf (stderr, "%s:%d: %s!\n", f, l, str);

  /* free the output string and end processing the arguments list. */
  free (str);
  va_end (vl);
}

/* pca_group_new: creates a pointer to a newly allocated group struct.
 * @name: the name of the new group.
 * @dims: the number of dimensions of the group.
 */
struct pca_group *pca_group_new (const char *name, unsigned long dims) {
  /* declare required variables. */
  struct pca_group *group;

  /* ensure the requested variables make sense. */
  if (!dims || !name || !strlen (name)) {
    /* output an error message and return nothing. */
    ferror ("invalid group parameters");
    return NULL;
  }

  /* allocate a pointer to the group struct. */
  group = (struct pca_group *) malloc (sizeof (struct pca_group));

  /* ensure the pointer was allocated. */
  if (!group) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate pca_group pointer");
    return NULL;
  }

  /* copy the group name string. */
  group->name = strdup (name);

  /* make sure the group name string was copied. */
  if (!group->name) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate group name");
    return NULL;
  }

  /* initialize the points array. */
  group->pts = NULL;
  group->n_pts = 0;

  /* allocate a mean vector. */
  group->mean = vector_new (dims);

  /* make sure we allocated a mean vector. */
  if (!group->mean) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate group mean vector");
    return NULL;
  }

  /* allocate a covariance matrix. */
  group->cov = matrix_new (dims, dims);

  /* make sure we allocated a covariance matrix. */
  if (!group->cov) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate group covariance matrix");
    return NULL;
  }

  /* allocate an eigenvalue vector. */
  group->eig_l = vector_new (dims);

  /* make sure we allocated an eigenvalue vector. */
  if (!group->eig_l) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate group eigenvalue vector");
    return NULL;
  }

  /* allocate an eigenvector matrix. */
  group->eig_v = matrix_new (dims, dims);

  /* make sure we allocated an eigenvector matrix. */
  if (!group->eig_v) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate group eigenvector matrix");
    return NULL;
  }

  /* initialize the pdf matrix. */
  group->map = NULL;
  group->mx = NULL;
  group->my = NULL;

  /* initialize the parenthood status. */
  group->is_parent = 0;
  group->child_a = 0;
  group->child_b = 0;

  /* initialize the membership info. */
  group->membstr = NULL;
  group->memb = NULL;
  group->n_memb = 0;

  /* return the allocated structure pointer. */
  return group;
}

/* forward declaration. */
void pca_group_heapsort_heapify (char **memb, unsigned long n,
                                 unsigned long i);

/* pca_group_heapsort_flatten: flattens a string array to a string.
 * @memb: the string array.
 * @n: the string array size.
 */
char *pca_group_heapsort_flatten (char **memb, unsigned long n) {
  /* declare required variables. */
  unsigned long i, len;
  char *flat;

  /* determine the flat length. */
  for (len = 2, i = 0; i < n; i++) len += 4 + strlen (memb[i]);

  /* allocate the flat string. */
  flat = (char *) malloc (sizeof (char) * len);
  strcpy (flat, "");

  /* append each string. */
  for (i = 0; i < n; i++) {
    /* append the string in single quotes. */
    strcat (flat, "'");
    strcat (flat, memb[i]);
    strcat (flat, "'");

    /* add a comma if necessary. */
    if (i < n - 1) strcat (flat, ",");
  }

  /* return the flat string. */
  return flat;
}

/* pca_group_heapsort_deroot: pops the top element off a heap.
 * @memb: the string array.
 * @n: pointer to the array size.
 */
char *pca_group_heapsort_deroot (char **memb, unsigned long *n) {
  /* declare required variables. */
  char *sroot;

  /* get the value at the root. */
  sroot = memb[0];

  /* set the new root value to that of the smallest leaf. */
  memb[0] = memb[*n - 1];
  (*n)--;

  /* re-heapify. */
  pca_group_heapsort_heapify (memb, *n, 0);

  /* return the root value. */
  return sroot;
}

/* pca_group_heapsort_heapify: heapifies a portion of an array.
 * @memb: the string array.
 * @n: the string array size.
 */
void pca_group_heapsort_heapify (char **memb, unsigned long n,
                                 unsigned long i) {
  /* declare required variables. */
  unsigned long left, right, max;
  char *sleft, *sright, *smax, *swap;

  /* set up the indices. */
  left = (2 * i) + 1;
  right = (2 * i) + 2;
  max = i;
  swap = NULL;

  /* get the values at those indices. */
  sleft = left < n ? memb[left] : NULL;
  sright = right < n ? memb[right] : NULL;
  smax = max < n ? memb[max] : NULL;

  /* de we have a larger values on the left node? */
  if (left < n && strcmp (sleft, smax) > 0) {
    /* swap indices. */
    max = left;
    smax = sleft;
  }

  /* do we have a larger values on the right node? */
  if (right < n && strcmp (sright, smax) > 0) {
    /* swap indices. */
    max = right;
    smax = sright;
  }

  /* did we perform a swap? */
  if (max != i) {
    /* actually swap the values. */
    swap = memb[i];
    memb[i] = smax;
    memb[max] = swap;

    /* heapify the child subtree. */
    pca_group_heapsort_heapify (memb, n, max);
  }
}

/* pca_group_heapsort_build: heapifies an array of strings.
 * @memb: the string array.
 * @n: the string array size.
 */
void pca_group_heapsort_build (char **memb, unsigned long n) {
  /* declare required variables. */
  signed long i;

  /* heapify each leaf node, working back up to the root. */
  for (i = (signed long) floor ((double) n / 2.0) - 1; i >= 0; i--)
    pca_group_heapsort_heapify (memb, n, (unsigned long) i);
}

/* pca_group_union_new: creates a union of two groups as a new group.
 * @P: the pca structure to read from.
 * @ga: the first group index to use in the union.
 * @gb: the second group index to use in the union.
 */
struct pca_group *pca_group_union_new (struct pca *P,
                                       unsigned long ga,
                                       unsigned long gb) {
  /* declare required variables. */
  unsigned long i, j, u_ga, u_gb;
  struct pca_group *group;
  double fa, fb;

  /* allocate a pointer to the group struct. */
  group = (struct pca_group *) malloc (sizeof (struct pca_group));

  /* ensure the pointer was allocated. */
  if (!group) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate pca_group union pointer");
    return NULL;
  }

  /* allocate memory for the group name. */
  group->name = (char *) malloc (sizeof (char) *
    (strlen (P->groups[ga]->name) + strlen (P->groups[gb]->name) + 16));

  /* make sure the string was allocated. */
  if (!group->name) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate pca group union name string");
    return NULL;
  }

  /* determine the order in which to place the child strings. */
  if (strcmp (P->groups[ga]->name, P->groups[gb]->name) > 0) {
    /* regular. */
    u_ga = ga;
    u_gb = gb;
  }
  else {
    /* reverse. */
    u_ga = gb;
    u_gb = ga;
  }

  /* build the name string. */
  sprintf (group->name, "union(%s%s%s,%s%s%s)",
    P->groups[u_ga]->is_parent ? "" : "'",
    P->groups[u_ga]->name,
    P->groups[u_ga]->is_parent ? "" : "'",
    P->groups[u_gb]->is_parent ? "" : "'",
    P->groups[u_gb]->name,
    P->groups[u_gb]->is_parent ? "" : "'");

  /* initialize the points array. */
  group->n_pts = P->groups[ga]->n_pts + P->groups[gb]->n_pts;
  group->pts = NULL;

  /* initialize the group statistics. */
  group->mean = NULL;
  group->cov = NULL;

  /* initialize the group eigenvalues and eigenvectors. */
  group->eig_l = NULL;
  group->eig_v = NULL;

  /* calculate the fractions for unweighted sums. */
  fa = ((double) P->groups[ga]->n_pts) / ((double) group->n_pts);
  fb = ((double) P->groups[gb]->n_pts) / ((double) group->n_pts);

  /* initialize the group pdf matrix. */
  if (P->groups[ga]->map && P->groups[gb]->map) {
    /* allocate the matrix. */
    group->map = matrix_new (P->groups[ga]->map->m,
                             P->groups[ga]->map->n);

    /* loop over the matrix rows. */
    for (i = 0; i < group->map->m; i++) {
      /* and over the matrix columns. */
      for (j = 0; j < group->map->n; j++) {
        /* calculate the average pdf value for the sum. */
        group->map->d[i][j] =
          fa * P->groups[ga]->map->d[i][j] +
          fb * P->groups[gb]->map->d[i][j];
      }
    }

    /* copy the x,y values. */
    group->mx = vector_copy (P->groups[ga]->mx);
    group->my = vector_copy (P->groups[ga]->my);
  }
  else {
    /* initialize the pdf matrix to null. */
    group->map = NULL;
    group->mx = NULL;
    group->my = NULL;
  }

  /* initialize the group index information. */
  group->is_parent = 1;
  group->child_a = ga;
  group->child_b = gb;

  /* initialize the group membership information. */
  group->n_memb =
    (P->groups[ga]->n_memb ? P->groups[ga]->n_memb : 1) +
    (P->groups[gb]->n_memb ? P->groups[gb]->n_memb : 1);
  group->memb = (char **) malloc (sizeof (char *) * group->n_memb);

  /* ensure we allocated group member strings. */
  if (!group->memb) {
    /* ouptut an error message and return nothing. */
    ferror ("failed to allocate union membership array");
    return NULL;
  }

  /* init the index counter. */
  j = 0;

  /* add the first group members. */
  if (P->groups[ga]->n_memb) {
    /* add the strings from the first group's membership array. */
    for (i = 0; i < P->groups[ga]->n_memb; i++, j++)
      group->memb[j] = strdup (P->groups[ga]->memb[i]);
  }
  else {
    /* add the first group name. */
    group->memb[j] = strdup (P->groups[ga]->name);
    j++;
  }

  /* add the second group members. */
  if (P->groups[gb]->n_memb) {
    /* add the strings from the second group's membership array. */
    for (i = 0; i < P->groups[gb]->n_memb; i++, j++)
      group->memb[j] = strdup (P->groups[gb]->memb[i]);
  }
  else {
    /* add the first group name. */
    group->memb[j] = strdup (P->groups[gb]->name);
    j++;
  }

  /* heapify the string array. */
  pca_group_heapsort_build (group->memb, group->n_memb);

  /* sort the string array. */
  j = group->n_memb;
  for (i = 0; i < group->n_memb; i++) {
    /* pop the next name off the heap. */
    group->memb[group->n_memb - i - 1] =
      pca_group_heapsort_deroot (group->memb, &j);
  }

  /* flatten the string array. */
  group->membstr = pca_group_heapsort_flatten (group->memb, group->n_memb);

  /* return the group. */
  return group;
}

/* pca_new: parses a SIMCA+ pca scores file into a pca structure.
 * @fname: the input text filename to open.
 */
struct pca *pca_new (const char *fname) {
  /* declare required variables. */
  char buf[N_BUF], *pbuf, *pbuf_r;
  unsigned long i, j, g;
  struct pca *P;
  int uniq;

  /* ensure the filename is valid. */
  if (!fname || !strlen (fname)) {
    /* output an error message and return nothing. */
    ferror ("invalid pca filename");
    return NULL;
  }

  /* allocate a pca structure pointer. */
  P = (struct pca *) malloc (sizeof (struct pca));

  /* make sure we allocated the structure pointer. */
  if (!P) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate pca pointer");
    return NULL;
  }

  /* store the filename string. */
  P->fname = strdup (fname);

  /* make sure we stored the filename. */
  if (!P->fname) {
    /* output an error message and return nothing. */
    ferror ("failed to copy pca filename");
    return NULL;
  }

  /* open the input text file. */
  P->fh = fopen (P->fname, "rb");

  /* make sure we opened the input text file. */
  if (!P->fh) {
    /* output an error message and return nothing. */
    ferror ("failed to open file '%s'", P->fname);
    return NULL;
  }

  /* read the first line of text into the buffer. */
  if (!fgets (buf, N_BUF, P->fh)) {
    /* output an error message and return nothing. */
    ferror ("failed to read file header");
    return NULL;
  }

  /* strip newline characters off the end of the buffer. */
  while (buf[strlen (buf) - 1] == '\n' || buf[strlen (buf) - 1] == '\r')
    buf[strlen (buf) - 1] = '\0';

  /* initialize the header string array. */
  P->headers = NULL;
  P->n_headers = 0;

  /* get the first string token from the header line. */
  pbuf = strtok_r (buf, DELIM, &pbuf_r);

  /* loop through the string tokens. */
  while (pbuf) {
    /* increment the header count. */
    P->n_headers++;

    /* resize the header array accordingly. */
    if (P->headers) {
      /* reallocate an existing header array. */
      P->headers = (char **)
        realloc (P->headers, sizeof (char *) * P->n_headers);
    }
    else {
      /* allocate a new header array. */
      P->headers = (char **)
        malloc (sizeof (char *) * P->n_headers);
    }

    /* store the header string in the array. */
    P->headers[P->n_headers - 1] = strdup (pbuf);

    /* get the next string token. */
    pbuf = strtok_r (NULL, DELIM, &pbuf_r);
  }

  /* make sure that we parsed enough header strings, and that the first
   * two header strings are integer and string identifiers of the points.
   */
  if (P->n_headers < 2 ||
      !strstr (P->headers[0], "Obs") ||
      !strstr (P->headers[1], "Obs")) {
    /* output an error message and return nothing. */
    ferror ("primary and secondary IDs required");
    return NULL;
  }

  /* loop through the remaining headers. */
  for (i = 2; i < P->n_headers; i++) {
    /* and all other headers. */
    for (j = i + 1; j < P->n_headers; j++) {
      /* make sure the header strings are unique. sometimes users will
       * export PC1 vs PC1, and wonder why their scores collapsed to a
       * line. this will give them a hint before they fail too much.
       */
      if (strcmp (P->headers[i], P->headers[j]) == 0) {
        /* output an error message. */
        ferror ("data columns must be unique! ('%s'[%lu] == '%s'[%lu])",
          P->headers[i], i + 1, P->headers[j], j + 1);

        /* return nothing. */
        return NULL;
      }
    }
  }

  /* store the number of indexes to create. */
  P->n_indexes = P->n_headers - 2;

  /* make sure we have enough indices. */
  if (P->n_indexes < 2) {
    /* output an error message. */
    ferror ("at least two columns of data required (found %lu)",
            P->n_indexes);

    /* return nothing. */
    return NULL;
  }

  /* allocate memory for the index array. */
  P->indexes = (unsigned long *)
    malloc (sizeof (unsigned long) * P->n_indexes);

  /* make sure we allocated the index array memory. */
  if (!P->indexes) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate dimension index array");
    return NULL;
  }

  /* store the column index of the given data dimension. */
  for (i = 0; i < P->n_indexes; i++) P->indexes[i] = i + 2;

  /* initialize the group array. */
  P->groups = NULL;
  P->n_groups = 0;

  /* loop until we've read the entire input file. */
  while (!feof (P->fh)) {
    /* get the next line of the input file. */
    if (fgets (buf, N_BUF, P->fh)) {
      /* strip newline characters off the end of the buffer string. */
      while (buf[strlen (buf) - 1] == '\n' || buf[strlen (buf) - 1] == '\r')
        buf[strlen (buf) - 1] = '\0';

      /* discard the first column value, assumed to be the integer index. */
      pbuf = strtok_r (buf, DELIM, &pbuf_r);
      if (pbuf) pbuf = strtok_r (NULL, DELIM, &pbuf_r);

      /* if we parsed a buffer token. */
      if (pbuf) {
        /* see if the current group name (buffer token) has been seen yet. */
        for (uniq = 1, g = 0; g < P->n_groups; g++) {
          /* is the currently indexed group name the same as the token? */
          if (strcmp (P->groups[g]->name, pbuf) == 0) {
            /* the token is not unique. break the loop. */
            uniq = 0;
            break;
          }
        }

        /* are we sitting on a unique new group name? */
        if (uniq) {
          /* increment the group count. */
          P->n_groups++;

          /* resize the group array accordingly. */
          if (P->groups) {
            /* reallocate the currently existing group array. */
            P->groups = (struct pca_group **)
              realloc (P->groups, sizeof (struct pca_group *) * P->n_groups);
          }
          else {
            /* allocate a new group array. */
            P->groups = (struct pca_group **)
              malloc (sizeof (struct pca_group *) * P->n_groups);
          }

          /* store the new group at the new index. */
          P->groups[P->n_groups - 1] = pca_group_new (pbuf, P->n_indexes);

          /* make sure we created the new group. */
          if (!P->groups[P->n_groups - 1]) {
            /* output an error message and return nothing. */
            ferror ("failed to allocate new group '%s'", pbuf);
            return NULL;
          }
        }
      }
    }
  }

  /* loop through each unique group we created. */
  for (g = 0; g < P->n_groups; g++) {
    /* initialize the group points array. */
    P->groups[g]->pts = NULL;
    P->groups[g]->n_pts = 0;

    /* rewind to the beginning of the input file. */
    fseek (P->fh, 0L, SEEK_SET);

    /* loop through the entire input file. */
    while (!feof (P->fh)) {
      if (fgets (buf, N_BUF, P->fh)) {
        /* strip newline characters from the end of the buffer. */
        while (buf[strlen (buf) - 1] == '\n' || buf[strlen (buf) - 1] == '\r')
          buf[strlen (buf) - 1] = '\0';

        /* discard the first column value, assumed to be the integer index. */
        pbuf = strtok_r (buf, DELIM, &pbuf_r);
        if (pbuf) pbuf = strtok_r (NULL, DELIM, &pbuf_r);

        /* do we have a point in a group that matches the current one? */
        if (pbuf && strcmp (P->groups[g]->name, pbuf) == 0) {
          /* increment the current group points count. */
          P->groups[g]->n_pts++;

          /* resize the points array accordingly. */
          if (P->groups[g]->pts) {
            /* reallocate an already existing points array. */
            P->groups[g]->pts = (struct vector **)
              realloc (P->groups[g]->pts,
                sizeof (struct vector *) * P->groups[g]->n_pts);
          }
          else {
            /* allocate a new points array. */
            P->groups[g]->pts = (struct vector **)
              malloc (sizeof (struct vector *) * P->groups[g]->n_pts);
          }

          /* allocate a new vector in the new points array. */
          P->groups[g]->pts[P->groups[g]->n_pts - 1] =
            vector_new (P->n_indexes);

          /* make sure we allocated the new vector. */
          if (!P->groups[g]->pts[P->groups[g]->n_pts - 1]) {
            /* output an error message. */
            ferror ("failed to allocate vector %lu in group '%s'",
                    P->groups[g]->n_pts - 1, P->groups[g]->name);

            /* return nothing. */
            return NULL;
          }

          /* get the next buffer token (column). */
          pbuf = strtok_r (NULL, DELIM, &pbuf_r);
          i = 0;

          /* loop through the remaining columns. */
          while (pbuf) {
            /* make sure the column index is in bounds. */
            if (i >= P->n_indexes) {
              /* output an error message. */
              ferror ("data column index out of bounds! (%lu > %lu)",
                      i + 1, P->n_indexes);

              /* return nothing. */
              return NULL;
            }

            /* get the value of the current column into the vector. */
            P->groups[g]->pts[P->groups[g]->n_pts - 1]->d[i] = atof (pbuf);

            /* get the next buffer token (column). */
            pbuf = strtok_r (NULL, DELIM, &pbuf_r);
            i++;
          }
        }
      }
    }
  }

  /* close the input text file handle. */
  fclose (P->fh);

  /* loop through the groups. */
  for (g = 0; g < P->n_groups; g++) {
    /* make sure the group has more points than dimensions. */
    if (P->groups[g]->n_pts <= pca_dimensions (P)) {
      /* output an error message. */
      ferror ("too few points (%lu <= %lu) in group '%s'",
        P->groups[g]->n_pts, pca_dimensions (P),
        P->groups[g]->name);

      /* return nothing. */
      return NULL;
    }
  }

  /* try to calculate the group means. */
  if (!pca_calculate_means (P)) {
    /* output an error message and return nothing. */
    ferror ("failed to calculate means");
    return NULL;
  }

  /* try to calculate the group covariances. */
  if (!pca_calculate_covs (P)) {
    /* output an error message and return nothing. */
    ferror ("failed to calculate covariances");
    return NULL;
  }

  /* try to calculate the group eigenvalues/eigenvectors. */
  if (!pca_calculate_eigs (P)) {
    /* output an error message and return nothing. */
    ferror ("failed to calculate eigendecomposition");
    return NULL;
  }

  /* return the new pca structure. */
  return P;
}

/* pca_group_free: free a pca group structure.
 * @group: the pca group structure to free.
 */
void pca_group_free (struct pca_group *group) {
  /* declare required variables. */
  unsigned long i, p;

  /* don't free a null structure. */
  if (!group) return;

  /* free the group name. */
  free (group->name);
  group->name = NULL;

  /* is the points array allocated? */
  if (group->pts) {
    /* free the points in the group. */
    for (p = 0; p < group->n_pts; p++)
      vector_free (group->pts[p]);

    /* free the points array. */
    free (group->pts);
    group->pts = NULL;
  }

  /* zero the points count. */
  group->n_pts = 0;

  /* free the mean and covariance structures. */
  vector_free (group->mean);
  matrix_free (group->cov);

  /* free the eigendecomposition structures. */
  vector_free (group->eig_l);
  matrix_free (group->eig_v);

  /* free the pdf matrix. */
  if (group->map) {
    /* free all the related structs. */
    matrix_free (group->map);
    vector_free (group->mx);
    vector_free (group->my);
  }

  /* free the membership string array. */
  if (group->memb) {
    /* free the strings in the array. */
    for (i = 0; i < group->n_memb; i++) {
      /* free the string. */
      free (group->memb[i]);
      group->memb[i] = NULL;
    }

    /* free the array. */
    free (group->memb);
    group->memb = NULL;
    group->n_memb = 0;
  }

  /* free the membership string. */
  free (group->membstr);
  group->membstr = NULL;

  /* free the group structure pointer. */
  free (group);
  group = NULL;
}

/* pca_free: free a pca structure.
 * @P: the pca structure to free.
 */
void pca_free (struct pca *P) {
  /* declare required variables. */
  unsigned long g, h;

  /* don't free a null structure. */
  if (!P) return;

  /* loop through the headers. */
  for (h = 0; h < P->n_headers; h++) {
    /* free the currently indexed header. */
    free (P->headers[h]);
    P->headers[h] = NULL;
  }

  /* free the header array. */
  free (P->headers);
  P->headers = NULL;

  /* free the index array. */
  free (P->indexes);
  P->indexes = NULL;

  /* free the groups in the group array. */
  for (g = 0; g < P->n_groups; g++)
    pca_group_free (P->groups[g]);

  /* free the groups array. */
  free (P->groups);
  P->groups = NULL;
  P->n_groups = 0;

  /* null out the filename and file handle. */
  free (P->fname);
  P->fname = NULL;
  P->fh = NULL;

  /* free the pca structure pointer. */
  free (P);
  P = NULL;
}

/* pca_group_copy: copies an exact replica of a group.
 * @group: the group to copy.
 */
struct pca_group *pca_group_copy (struct pca_group *group) {
  /* declare required variables. */
  struct pca_group *gcopy;
  unsigned long i;

  /* make sure the input group is defined. */
  if (!group) {
    /* output an error and return nothing. */
    ferror ("input group structure to copy is undefined");
    return NULL;
  }

  /* build a new group structure. */
  gcopy = pca_group_new (group->name, group->pts[0]->n);

  /* make sure we got a new group. */
  if (!gcopy) {
    /* output an error and return nothing. */
    ferror ("failed to initialize group copy");
    return NULL;
  }

  /* allocate memory for the points array. */
  gcopy->n_pts = group->n_pts;
  gcopy->pts = (struct vector **)
    malloc (sizeof (struct vector *) * gcopy->n_pts);

  /* make sure the points array was allocated. */
  if (!gcopy->pts) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate copied points array");
    return NULL;
  }

  /* loop through the input points array. */
  for (i = 0; i < group->n_pts; i++) {
    /* copy the point over. */
    gcopy->pts[i] = vector_copy (group->pts[i]);

    /* make sure the point copied ok. */
    if (!gcopy->pts[i]) {
      /* output an error message and return nothing. */
      ferror ("failed to copy point %lu from group '%s'", i, group->name);
      return NULL;
    }
  }

  /* copy the statistics. no sense in recalculating them. */
  gcopy->mean = vector_copy (group->mean);
  gcopy->cov = matrix_copy (group->cov);

  /* copy the eigendecomposition. */
  gcopy->eig_l = vector_copy (group->eig_l);
  gcopy->eig_v = matrix_copy (group->eig_v);

  /* copy the pdf matrix. */
  if (group->map) {
    /* copy all related structs. */
    gcopy->map = matrix_copy (group->map);
    gcopy->mx = vector_copy (group->mx);
    gcopy->my = vector_copy (group->my);
  }

  /* make sure everything copied ok. */
  if (!gcopy->mean || !gcopy->cov ||
      !gcopy->eig_l || !gcopy->eig_v) {
    /* output an error message and return nothing. */
    ferror ("failed to copy statistics for '%s' copy", group->name);
    return NULL;
  }

  /* copy the flat member string. */
  if (group->membstr) gcopy->membstr = strdup (group->membstr);

  /* copy the member strings. */
  if (group->memb && group->n_memb) {
    /* copy all related structs. */
    gcopy->n_memb = group->n_memb;
    gcopy->memb = (char **) malloc (sizeof (char *) * gcopy->n_memb);

    /* loop through the strings. */
    for (i = 0; i < gcopy->n_memb; i++)
      gcopy->memb[i] = strdup (group->memb[i]);
  }

  /* return the copied group. */
  return gcopy;
}

/* pca_group_copy_bootstrapped: copies a bootstrapped version of a group.
 * @group: the group to bootstrap-copy.
 */
struct pca_group *pca_group_copy_bootstrapped (struct pca_group *group) {
  /* declare required variables. */
  struct pca_group *gboot;
  unsigned long i, pidx;

  /* make sure the input group is defined. */
  if (!group) {
    /* output an error and return nothing. */
    ferror ("input group structure to bootstrap copy is undefined");
    return NULL;
  }

  /* build a new group structure. */
  gboot = pca_group_new (group->name, group->pts[0]->n);

  /* make sure we got a new group. */
  if (!gboot) {
    /* output an error and return nothing. */
    ferror ("failed to initialize bootstrapped group copy");
    return NULL;
  }

  /* allocate memory for the new point set. */
  gboot->n_pts = group->n_pts;
  gboot->pts = (struct vector **)
    malloc (sizeof (struct vector *) * gboot->n_pts);

  /* make sure we allocated a new array. */
  if (!gboot->pts) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate points array for bootstrapped group copy");
    return NULL;
  }

  /* build the same number of points as the original group. */
  for (i = 0; i < gboot->n_pts; i++) {
    /* get the point index to copy. */
    pidx = (randull ()) % group->n_pts;

    /* copy the vector into the new point slot. */
    gboot->pts[i] = vector_copy (group->pts[pidx]);

    /* make sure we copied the vector. */
    if (!gboot->pts[i]) {
      /* output an error message. */
      ferror ("failed to copy point %lu from group '%s'", pidx, group->name);
      return NULL;
    }
  }

  /* calculate the statistics again. */
  if (!pca_group_calculate_mean (gboot) ||
      !pca_group_calculate_cov (gboot) ||
      !pca_group_calculate_eig (gboot)) {
    /* output an error message. */
    ferror ("failed to recalculate bootstrapped statistics for '%s' copy",
      group->name);

    /* return nothing. */
    return NULL;
  }

  /* don't worry about the pdf matrix. */
  gboot->map = NULL;
  gboot->mx = NULL;
  gboot->my = NULL;

  /* return the new group. */
  return gboot;
}

/* pca_group_del_point: deletes a point from a given pca group structure.
 * @group: the pca group structure to operate on.
  @idx: the point index to delete.
 */
int pca_group_del_point (struct pca_group *group, unsigned long idx) {
  /* declare required variables. */
  unsigned long i;

  /* ensure the group index is in bounds. */
  if (idx >= group->n_pts) {
    /* output an error message and return failure. */
    ferror ("group index out of bounds (%lu >= %lu)", idx, group->n_pts);
    return 0;
  }

  /* free the point. */
  vector_free (group->pts[idx]);

  /* loop through the remaining points. */
  for (i = idx + 1; i < group->n_pts; i++) {
    /* shift the point up. */
    group->pts[i - 1] = group->pts[i];
  }

  /* resize the points array. */
  group->n_pts--;
  group->pts = (struct vector **) realloc (group->pts,
    sizeof (struct vector *) * group->n_pts);

  /* return success. */
  return 1;
}

/* pca_group_calculate_mean: calculate the mean vector of a given pca
 * group structure.
 * @group: the pca group structure to operate on.
 */
int pca_group_calculate_mean (struct pca_group *group) {
  /* declare required variables. */
  unsigned long i, d, D;

  /* ensure the group contains sufficient points to calculate a mean. */
  if (!group->n_pts) {
    /* output an error message and return failure. */
    ferror ("group '%s' has no points, cannot calculate mean", group->name);
    return 0;
  }

  /* get the dimensionality of the points. */
  D = group->pts[0]->n;

  /* loop through the dimensions of the points. */
  for (d = 0; d < D; d++) {
    /* initialize the vector element. */
    group->mean->d[d] = 0.0;

    /* loop through the points, summing into the vector element. */
    for (i = 0; i < group->n_pts; i++)
      group->mean->d[d] += group->pts[i]->d[d];

    /* normalize to yield a mean. */
    group->mean->d[d] /= (double) group->n_pts;
  }

  /* return success. */
  return 1;
}

/* pca_group_calculate_cov: calculate the covariance matrix of a given
 * pca group structure.
 * @group: the pca group structure to operate on.
 */
int pca_group_calculate_cov (struct pca_group *group) {
  /* declare required variables. */
  unsigned long i, j, k, D;

  /* ensure the group contains sufficient points to calculate covariance. */
  if (group->n_pts < 2) {
    /* output an error message. */
    ferror ("group '%s' has too few points, cannot calculate covariance",
            group->name);

    /* return failure. */
    return 0;
  }

  /* get the dimensionality of the points. */
  D = group->pts[0]->n;

  /* loop through the dimensions once. */
  for (j = 0; j < D; j++) {
    /* and loop through the dimensions again. */
    for (k = j; k < D; k++) {
      /* initialize the matrix element. */
      group->cov->d[j][k] = 0.0;

      /* loop through the points in the current dimension pair. */
      for (i = 0; i < group->n_pts; i++) {
        /* sum up into the matrix element. */
        group->cov->d[j][k] += (group->pts[i]->d[j] - group->mean->d[j])
                             * (group->pts[i]->d[k] - group->mean->d[k]);
      }

      /* normalize the matrix element to a variance. */
      group->cov->d[j][k] /= (double) (group->n_pts - 1);

      /* are we off-diagonal? */
      if (j != k) group->cov->d[k][j] = group->cov->d[j][k];
    }
  }

  /* return success. */
  return 1;
}

/* pca_group_calculate_eig: eigendecompose the covariance matrix of a
 * given pca group structure.
 * @group: the pca group structure to operate on.
 */
int pca_group_calculate_eig (struct pca_group *group) {
  /* calculate the matrix eigendecomposition of the group covariance
   * matrix, storing the values in the appropriate group structures.
   */
  if (!matrix_eig (group->cov, group->eig_l, group->eig_v)) {
    /* output an error message. */
    ferror ("failed to calculate eigendecomposition of cov['%s']",
            group->name);

    /* return failure. */
    return 0;
  }

  /* return success. */
  return 1;
}

/* pca_group_calculate_map: builds a gaussian mixture model pdf matrix.
 * @group: the pca group structure to operate on.
 * @x: a vector of x-values to calculate the map for.
 * @y: a vector of y-values to calculate the map for.
 * @s: the desired standard deviation of the smoothing kernel.
 */
int pca_group_calculate_map (struct pca_group *group,
                             struct vector *x,
                             struct vector *y,
                             double s) {
  /* declare required variables. */
  unsigned long i, j, k;

  /* allocate the group pdf matrix. */
  group->map = matrix_new (x->n, y->n);

  /* ensure the pdf matrix was allocated. */
  if (!group->map) {
    /* output an error message and return failure. */
    ferror ("failed to allocate pdf matrix for '%s'", group->name);
    return 0;
  }

  /* copy the map x,y values. */
  group->mx = vector_copy (x);
  group->my = vector_copy (y);

  /* make a loop through the x values. */
  for (j = 0; j < x->n; j++) {
    /* and the y values. */
    for (k = 0; k < y->n; k++) {
      /* loop through the points in the group. */
      for (group->map->d[j][k] = 0.0, i = 0; i < group->n_pts; i++) {
        /* add the current point's contribution. */
        group->map->d[j][k] +=
          (1.0 / (2.0 * M_PI * s * s * ((double) group->n_pts))) *
          exp (-0.5 * (pow ((group->pts[i]->d[0] - x->d[j]) / s, 2.0)
                     + pow ((group->pts[i]->d[1] - y->d[k]) / s, 2.0)));
      }
    }
  }

  /* return success. */
  return 1;
}

/* pca_calculate_means: calculate the mean vectors of the groups in the
 * pca structure P.
 * @P: the pca structure to operate on.
 */
int pca_calculate_means (struct pca *P) {
  /* declare required variables. */
  unsigned long g;

  /* loop through each group in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* calculate the current group mean. */
    if (!pca_group_calculate_mean (P->groups[g])) {
      /* output an error message and return failure. */
      ferror ("failed to calculate pca group %lu mean", g);
      return 0;
    }
  }

  /* return success. */
  return 1;
}

/* pca_calculate_covs: calculate the covariance matrices of the groups
 * in the pca structure P.
 * @P: the pca structure to operate on.
 */
int pca_calculate_covs (struct pca *P) {
  /* declare required variables. */
  unsigned long g;

  /* loop through each group in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* calculate the current group covariance. */
    if (!pca_group_calculate_cov (P->groups[g])) {
      /* output an error message and return failure. */
      ferror ("failed to calculate pca group %lu covariance", g);
      return 0;
    }
  }

  /* return success. */
  return 1;
}

/* pca_calculate_eigs: eigendecompose the covariance matrices of the groups
 * in the pca structure P.
 * @P: the pca structure to operate on.
 */
int pca_calculate_eigs (struct pca *P) {
  /* declare required variables. */
  unsigned long g;

  /* loop through each group in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* calculate the current group eigendecomposition. */
    if (!pca_group_calculate_eig (P->groups[g])) {
      /* output an error message and return failure. */
      ferror ("failed to calculate pca group %lu eigendecomposition", g);
      return 0;
    }
  }

  /* return success. */
  return 1;
}

/* pca_calculate_maps: calculates the (2D) pdf matrices of a pca dataset.
 * @P: the pca structure to operate on.
 */
int pca_calculate_maps (struct pca *P) {
  /* declare required variables. */
  struct vector *x, *y;
  unsigned long ndiv, i, g;
  double ss, ssx, ssy;

  /* get the number of divisions. */
  ndiv = pca_divisions ();

  /* calculate the boundaries. */
  pca_calculate_var (P, &ssx, &ssy);
  ss = sqrt (ssx + ssy);

  /* scale the boundaries. */
  ssx = sqrt (2.0 * ssx * CHISQ2);
  ssy = sqrt (2.0 * ssy * CHISQ2);

  /* allocate the value vectors. */
  x = vector_new (ndiv);
  y = vector_new (ndiv);

  /* ensure the vectors were allocated. */
  if (!x || !y) {
    /* output an error message and return failure. */
    ferror ("failed to allocate boundary vectors");
    return 0;
  }

  /* loop through the tics. */
  for (i = 0; i < ndiv; i++) {
    /* calculate the values. */
    x->d[i] = 2.0 * ssx * (((double) i) / ((double) ndiv)) - ssx;
    y->d[i] = 2.0 * ssy * (((double) i) / ((double) ndiv)) - ssy;
  }

  /* loop through the groups. */
  for (g = 0; g < P->n_groups; g++) {
    /* calculate the current group pdf matrix. */
    if (!pca_group_calculate_map (P->groups[g], x, y, 0.25 * ss)) {
      /* output an error message and return failure. */
      ferror ("failed to calculate pdf matrix for group %lu", g);
      return 0;
    }
  }

  /* free the vectors. */
  vector_free (x);
  vector_free (y);

  /* return success. */
  return 1;
}

/* pca_calculate_var: calculates the (2D) variances of the entire dataset.
 * @P: the pca structure to operate on.
 * @vx: a pointer to the output x-axis variance.
 * @vy: a pointer to the output y-axis variance.
 */
int pca_calculate_var (struct pca *P, double *vx, double *vy) {
  /* declare required variables. */
  unsigned long n, g, i;

  /* initialize the variances and the counter. */
  *vx = 0.0;
  *vy = 0.0;
  n = 0;

  /* loop through the groups in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* and through the points of the group. */
    for (i = 0; i < P->groups[g]->n_pts; i++) {
      /* sum up the squared distances from the origin. */
      *vx += pow (P->groups[g]->pts[i]->d[0], 2.0);
      *vy += pow (P->groups[g]->pts[i]->d[1], 2.0);

      /* increment the counter. */
      n++;
    }
  }

  /* finalize the variances. */
  *vx /= ((double) (n - 1));
  *vy /= ((double) (n - 1));

  /* return success. */
  return 1;
}

/* pca_find_group: returns the group structure having a given name in a
 * given pca structure.
 * @P: the pca structure to search.
 * @name: the name of the group to find.
 */
struct pca_group *pca_find_group (struct pca *P, const char *name) {
  /* declare required variables. */
  unsigned long g;

  /* loop through the groups in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* compare the current group name to that sought. */
    if (strcmp (P->groups[g]->name, name) == 0) {
      /* they match. return the group structure. */
      return P->groups[g];
    }
  }

  /* output an error message detailing our woeful failure. */
  ferror ("failed to find pca group with name '%s'", name);

  /* return nothing. */
  return NULL;
}

/* pca_find_group: returns the group structure having a given name in a
 * given pca structure.
 * @P: the pca structure to search.
 * @name: the name of the group to find.
 */
unsigned long pca_group_index (struct pca *P, const char *name) {
  /* declare required variables. */
  unsigned long g;

  /* loop through the groups in the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* compare the current group name to that sought. */
    if (strcmp (P->groups[g]->name, name) == 0) {
      /* they match. return the group structure. */
      return g;
    }
  }

  /* output an error message detailing our woeful failure. */
  ferror ("failed to find pca group with name '%s'", name);

  /* return nothing. */
  return 0;
}

/* pca_add_group: add a pca group structure to a pca structure.
 * @P: the pca structure to add the group to.
 * @group: the group to add.
 */
int pca_add_group (struct pca *P, struct pca_group *group) {
  /* make sure the two pointers are valid. */
  if (!P || !group) {
    /* output an error message and return failure. */
    ferror ("pca or pca_group structures are invalid");
    return 0;
  }

  /* increment the group counter. */
  P->n_groups++;

  /* resize the groups array accordingly. */
  if (P->groups) {
    /* reallocate the existing groups array. */
    P->groups = (struct pca_group **)
      realloc (P->groups, sizeof (struct pca_group *) * P->n_groups);
  }
  else {
    /* allocate a new groups array. */
    P->groups = (struct pca_group **)
      malloc (sizeof (struct pca_group *) * P->n_groups);
  }

  /* store the new group in the array. */
  P->groups[P->n_groups - 1] = group;

  /* return success. */
  return 1;
}

/* pca_del_group: deletes a group from a pca structure.
 * @P: the pca structure to modify.
 * @gdel: the group index to remove.
 */
struct pca_group *pca_del_group (struct pca *P, unsigned long gdel) {
  /* declare required variables. */
  struct pca_group *group;
  unsigned long g;

  /* make sure the input structure is valid. */
  if (!P) {
    /* output an error message and return nothing. */
    ferror ("input pca structure is undefined");
    return NULL;
  }

  /* make sure the group index is in range. */
  if (gdel >= P->n_groups) {
    /* output an error message and return nothing. */
    ferror ("group index out of bounds (%lu >= %lu)", gdel, P->n_groups);
    return NULL;
  }

  /* keep a pointer to the group. */
  group = P->groups[gdel];

  /* loop through the groups. */
  for (g = gdel + 1; g < P->n_groups; g++) {
    /* move the group up in the list. */
    P->groups[g - 1] = P->groups[g];
  }

  /* resize the group array. */
  P->n_groups--;
  P->groups = (struct pca_group **)
    realloc (P->groups, sizeof (struct pca_group *) * P->n_groups);

  /* ensure we successfully reallocated. */
  if (!P->groups) {
    /* output an error message and return nothing. */
    ferror ("failed to resize pca groups array to %lu groups", P->n_groups);
    return NULL;
  }

  /* return the removed group. */
  return group;
}

/* pca_copy: copies an exact replica of a pca structure.
 * @P: the pca structure to copy.
 */
struct pca *pca_copy (struct pca *P) {
  /* declare required variables. */
  struct pca_group *gcopy;
  struct pca *Pcopy;
  unsigned long g;

  /* make sure the input structure is valid. */
  if (!P) {
    /* output an error message and return nothing. */
    ferror ("input pca structure to copy is undefined");
    return NULL;
  }

  /* allocate a new structure. */
  Pcopy = (struct pca *) malloc (sizeof (struct pca));

  /* make sure we allocated a new structure. */
  if (!Pcopy) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate copy of pca structure");
    return NULL;
  }

  /* initialize the information that can be pulled from the original. */
  Pcopy->headers = NULL;
  Pcopy->n_headers = 0;
  Pcopy->indexes = NULL;
  Pcopy->n_indexes = 0;
  Pcopy->fname = NULL;
  Pcopy->fh = NULL;

  /* initialize the groups array. */
  Pcopy->groups = NULL;
  Pcopy->n_groups = 0;

  /* loop through the groups in the original pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* make a copy of the group. */
    gcopy = pca_group_copy (P->groups[g]);

    /* make sure we got a new copy. */
    if (!gcopy) {
      /* output an error message and return nothing. */
      ferror ("failed to copy group '%s'", P->groups[g]->name);
      return NULL;
    }

    /* add a bootstrapped copy to the new struct. */
    if (!pca_add_group (Pcopy, gcopy)) {
      /* output an error message and return nothing. */
      ferror ("failed to add copied group '%s'", P->groups[g]->name);
      return NULL;
    }
  }

  /* return the new structure. */
  return Pcopy;
}

/* pca_copy_bootstrapped: copies a bootstrapped version of a pca structure.
 * @P: the pca structure to bootstrap-copy.
 */
struct pca *pca_copy_bootstrapped (struct pca *P) {
  /* declare required variables. */
  struct pca_group *gboot;
  struct pca *Pboot;
  unsigned long g;

  /* make sure the input structure is valid. */
  if (!P) {
    /* output an error message and return nothing. */
    ferror ("input pca structure to bootstrap copy is undefined");
    return NULL;
  }

  /* allocate a new structure. */
  Pboot = (struct pca *) malloc (sizeof (struct pca));

  /* make sure we allocated a new structure. */
  if (!Pboot) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate bootstrap copy of pca structure");
    return NULL;
  }

  /* initialize the information that can be pulled from the original. */
  Pboot->headers = NULL;
  Pboot->n_headers = 0;
  Pboot->indexes = NULL;
  Pboot->n_indexes = 0;
  Pboot->fname = NULL;
  Pboot->fh = NULL;

  /* initialize the groups array. */
  Pboot->groups = NULL;
  Pboot->n_groups = 0;

  /* loop through the groups in the original pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* make a bootstrapped group copy. */
    gboot = pca_group_copy_bootstrapped (P->groups[g]);

    /* make sure we got a new copy. */
    if (!gboot) {
      /* output an error message and return nothing. */
      ferror ("failed to copy a bootstrapped group '%s'", P->groups[g]->name);
      return NULL;
    }

    /* add a bootstrapped copy to the new struct. */
    if (!pca_add_group (Pboot, gboot)) {
      /* output an error message and return nothing. */
      ferror ("failed to add bootstrapped group '%s'", P->groups[g]->name);
      return NULL;
    }
  }

  /* return the new pca structure. */
  return Pboot;
}

/* pca_dimensions: returns the number of dimensions of a pca structure.
 * @P: the pca structure to read.
 */
unsigned long pca_dimensions (struct pca *P) {
  /* make sure the pca structure is defined. */
  if (!P) {
    /* output an error message and return nothing. */
    ferror ("input pca structure is undefined");
    return 0;
  }

  /* make sure the pca structure has groups. */
  if (!P->n_groups) {
    /* output an error message and return nothing. */
    ferror ("input pca structure has no groups");
    return 0;
  }

  /* make sure the first group has points. */
  if (!P->groups[0]->n_pts) {
    /* output an error message and return nothing. */
    ferror ("input pca structure group has no points");
    return 0;
  }

  /* return the number of dimensions. */
  return P->groups[0]->pts[0]->n;
}

/* pca_divisions: get the number of map divisions via getopt.
 */
unsigned long pca_divisions (void) {
  /* declare required variables. */
  unsigned long ndiv = 0;

  /* get the number of divisions. */
  ndiv = opts_geti (OPTS_S_DIVS);

  /* ensure a valid count. */
  if (!ndiv) {
    /* set the default. */
    ndiv = 128;
  }

  /* return the count. */
  return ndiv;
}

/* pca_datamatrix: combine groups' points into a single matrix.
 * @P: the pca structure to operate on.
 */
struct matrix *pca_datamatrix (struct pca *P) {
  /* declare required variables. */
  unsigned long m, n, g, i, j, k;
  struct matrix *M;

  /* get the number of columns. (dimension count) */
  n = pca_dimensions (P);

  /* get the number of rows. (total point count) */
  for (m = 0, g = 0; g < P->n_groups; g++) {
    /* add to the total count. */
    m += P->groups[g]->n_pts;
  }

  /* allocate memory for the matrix. */
  M = matrix_new (m, n);

  /* ensure we allocated a matrix. */
  if (!M) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate data matrix");
    return NULL;
  }

  /* loop once again through the groups. */
  for (k = 0, g = 0; g < P->n_groups; g++) {
    /* loop through the points in the group. */
    for (i = 0; i < P->groups[g]->n_pts; i++) {
      /* loop through the dimensions. */
      for (j = 0; j < n; j++) {
        M->d[k][j] = P->groups[g]->pts[i]->d[j];
      }

      /* increment the row index. */
      k++;
    }
  }

  /* return the new matrix. */
  return M;
}

