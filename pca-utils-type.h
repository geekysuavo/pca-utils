
/*
 * pca-utils-type.h: header code for custom structure types.
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

#ifndef __PCA_UTILS_TYPE__
#define __PCA_UTILS_TYPE__

/* macros to calculate maximum and minimum. */
#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)

/* list_t: a type definition required for recursion in the list structure.
 */
typedef struct list list_t;

/* list: linked list containing void pointers as data elements.
 * @n: the number of elements/links in the list.
 * @prev: the previous member of the list.
 * @next: the next member of the list.
 * @last: the last member of the list.
 * @cmpfn: the comparison function for list data.
 * @d: the list node data.
 */
struct list {
  unsigned long n;

  list_t *prev;
  list_t *next;
  list_t *last;

  int (*cmpfn) (void *d1, void *d2);

  void *d;
};

/* pca_group: structure containing all group-related information.
 * @name: the textual (short) description of the group.
 * @memb: the array of member group names, or NULL for leaf group.
 * @membstr: the flat string version of the string array.
 * @n_memb: the number of member group names, or zero for leaf group.
 * @pts: the points in scores space belonging to the group.
 * @n_pts: the number of points belonging to the group.
 * @mean: the mean value (centroid) of the points in the group.
 * @cov: the covariance matrix of the points in the group.
 * @eig_l: eigenvalues of the group covariance matrix.
 * @eig_v: eigenvectors of the group covariance matrix.
 * @map: group mixture probability density matrix (2D scores only).
 * @mx: the x-values in the pdf matrix.
 * @my: the y-values in the pdf matrix.
 * @is_parent: whether or not the group is composed of child groups.
 * @child_a: the first child group index from the pca structure.
 * @child_b: the second child group index from the pca structure.
 */
struct pca_group {
  char *name;

  char **memb, *membstr;
  unsigned long n_memb;

  struct vector **pts;
  unsigned long n_pts;

  struct vector *mean;
  struct matrix *cov;

  struct vector *eig_l;
  struct matrix *eig_v;

  struct matrix *map;
  struct vector *mx, *my;

  int is_parent;
  unsigned long child_a;
  unsigned long child_b;
};

/* pca_tree_t: a type definition required for recursion in the tree structure.
 */
typedef struct pca_tree pca_tree_t;

/* pca_tree: structure containing a hierarchical clustering of groups.
 * @name: the name of the current node in the tree.
 * @membstr: the flat member string for the node.
 * @p: the p-value of the bifurcation at the current node.
 * @dist_a: the distance to the first child of the current node.
 * @dist_b: the distance to the second child of the current node.
 * @na: the number of points in the first child.
 * @nb: the number of points in the second child.
 * @group: the pca group of the current node.
 * @parent: the parent tree node, or NULL at the root.
 * @child_a: the first child tree node.
 * @child_b: the second child tree node.
 * @x: the x-location of the node on a plot.
 * @y: the y-location of the node on a plot.
 * @dy: the height of a (non-leaf) node on a plot.
 * @idx: the index of the node. leaves only.
 */
struct pca_tree {
  char *name, *membstr;

  double p;
  double dist_a, dist_b;
  unsigned long na, nb;

  struct pca_group *group;

  pca_tree_t *parent;
  pca_tree_t *child_a;
  pca_tree_t *child_b;

  double x, y, dy;
  unsigned long idx;
};

/* pca: structure containing all groups of a given PCA analysis.
 * @headers: column heading strings from the input text file.
 * @n_headers: number of column heading strings.
 * @indexes: column indices of the data dimensions.
 * @n_indexes: number of data dimensions.
 * @groups: the groups of the multivariate analysis.
 * @n_groups: the number of groups in the analysis.
 * @fname: the filename of the file used to read/write the struct.
 * @fh: filesystem handle of the file used to read/write the struct.
 */
struct pca {
  char **headers;
  unsigned long n_headers;

  unsigned long *indexes;
  unsigned long n_indexes;

  struct pca_group **groups;
  unsigned long n_groups;

  char *fname;
  FILE *fh;
};

/* vector: structure that contains a n-dimensional vector of scalar values.
 * @n: dimensionality of the vector.
 * @d: array of double-precision floating point scalar values.
 */
struct vector {
  unsigned long n;
  double *d;
};

/* matrix: structure that contains an m-by-n dimensional matrix of reals.
 * @m: row-wise dimensionality of the matrix.
 * @n: column-wise dimensionality of the matrix.
 * @d: 2d array of double-precision floating point scalar values.
 */
struct matrix {
  unsigned long m, n;
  double **d;
};

/* define the pca group-to-group distance metric function prototype. */
typedef double (*pca_distance_metric) (struct pca_group *groupA,
                                       struct pca_group *groupB);

#endif

