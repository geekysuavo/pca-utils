
/*
 * pca-utils-tree.c: source code for dendrogram generation.
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

/* forest_new: create a linked list of trees from a pca structure.
 * @P: the pca structure to use.
 */
struct list *forest_new (struct pca *P) {
  /* declare required variables. */
  struct list *forest;
  struct pca_tree *T;
  unsigned long g;

  /* allocate the new list. */
  forest = list_new (NULL);

  /* ensure we successfully made the list. */
  if (!forest) {
    /* output an error message and return nothing. */
    ferror ("failed to create new list for tree storage");
    return NULL;
  }

  /* loop through the groups of the pca structure. */
  for (g = 0; g < P->n_groups; g++) {
    /* create a new tree node. */
    T = (struct pca_tree *) malloc (sizeof (struct pca_tree));

    /* make sure we created a new tree node. */
    if (!T) {
      /* output an error and return nothing. */
      ferror ("failed to allocate new tree node for storage");
      return NULL;
    }

    /* store the name of the tree node. */
    T->name = strdup (P->groups[g]->name);
    T->membstr = NULL;

    /* make sure we stored the name. */
    if (!T->name) {
      /* output an error message and return nothing. */
      ferror ("failed to copy group name to tree node");
      return NULL;
    }

    /* set the numeric values of the tree node. */
    T->p = T->dist_a = T->dist_b = 0.0;

    /* store a pointer to the group. */
    T->group = P->groups[g];

    /* set the linkage values of the tree node. */
    T->parent = T->child_a = T->child_b = NULL;

    /* add the tree to the new list. */
    if (!list_append (forest, T)) {
      /* output an error message and return nothing. */
      ferror ("failed to add tree '%s' to storage list", T->name);
      return NULL;
    }
  }

  /* return the new list. */
  return forest;
}

/* forest_search: searches the forest structure for a group node.
 * @forest: the forest structure to search.
 * @group: the pca group to search for a node of.
 */
struct pca_tree *forest_search (struct list *forest, struct pca_group *group) {
  /* declare required variables. */
  struct pca_tree *Titer;
  struct list *iter;

  /* loop through the list. */
  for (iter = forest->next; iter; iter = iter->next) {
    /* get a reference to the current tree. */
    Titer = (struct pca_tree *) iter->d;

    /* make sure we have a tree. */
    if (!Titer) {
      /* output an error message and return nothing. */
      ferror ("forest list node contains invalid tree pointer");
      return NULL;
    }

    /* see if the group of the current tree matches our search group. */
    if (strcmp (Titer->group->name, group->name) == 0) {
      /* we found it. return the tree structure. */
      return Titer;
    }
  }

  /* failure. */
  ferror ("failed to find matching tree for '%s' in forest", group->name);
  return NULL;
}

/* forest_union: joins two trees together in a forest structure.
 * @forest: the forest structure to manipulate.
 * @T: the parent tree node of the union.
 * @Ta: the first tree of the union.
 * @Tb: the second tree of the union.
 */
int forest_union (struct list *forest, struct pca_tree *T,
                  struct pca_tree *Ta, struct pca_tree *Tb) {
  /* delete the first list node. */
  if (!list_delete (forest, Ta)) {
    /* output an error message and return failure. */
    ferror ("failed to delete forest list node for '%s'", Ta->name);
    return 0;
  }

  /* delete the second list node. */
  if (!list_delete (forest, Tb)) {
    /* output an error message and return failure. */
    ferror ("failed to delete forest list node for '%s'", Tb->name);
    return 0;
  }

  /* link up the child nodes to the parent node. */
  T->child_a = Ta;
  T->child_b = Tb;

  /* link up the parent node to the child nodes. */
  Ta->parent = T;
  Tb->parent = T;

  /* store the new parent node in the forest. */
  list_append (forest, T);

  /* return success. */
  return 1;
}

/* pca_tree_new: generate a tree from a pca group distance matrix.
 * @P: the pca structure to use.
 * @D: the distance matrix to use.
 */
struct pca_tree *pca_tree_new (struct pca *P, struct matrix *D) {
  /* declare required variables. */
  struct pca_group *groupA, *groupB;
  unsigned long i, j;
  struct list *forest;
  struct pca_tree *T, *Ta, *Tb;
  struct matrix *Dop;
  struct pca *Pop;

  /* initialize the output tree. */
  T = NULL;

  /* make sure the input data structures are valid. */
  if (!P || !D) {
    /* output an error message and return nothing. */
    ferror ("input data structures to pca tree generation are invalid");
    return NULL;
  }

  /* copy the pca structure. */
  Pop = pca_copy (P);

  /* make sure we copied the pca structure. */
  if (!Pop) {
    /* output an error message and return nothing. */
    ferror ("failed to copy input pca structure");
    return NULL;
  }

  /* copy the distance matrix. */
  Dop = matrix_copy (D);

  /* make sure we copied the distance matrix. */
  if (!Dop) {
    /* output an error message and return nothing. */
    ferror ("failed to copy input distance matrix");
    return NULL;
  }

  /* build the list of singleton trees to start from. */
  forest = forest_new (P);

  /* make sure we allocated the forest. */
  if (!forest) {
    /* output an error message and return nothing. */
    ferror ("failed to build initial tree storage list");
    return NULL;
  }

  /* loop until the tree has as many leaves as the number of input groups. */
  while (Dop->m > 1) {
    /* find the closest pair of groups based on Dop. */
    pca_distances_min (Dop, &i, &j);

    /* store the group pointers before they're lost. */
    groupA = Pop->groups[i];
    groupB = Pop->groups[j];

    /* get the tree nodes from the forest. */
    Ta = forest_search (forest, groupA);
    Tb = forest_search (forest, groupB);

    /* get the new parent node and reduce the distance matrix. */
    T = pca_distances_union (Pop, Dop, i, j);

    /* store the group points counts in the tree node. */
    T->na = groupA->n_pts;
    T->nb = groupB->n_pts;

    /* ensure the matrix reduction step worked. */
    if (!T) {
      /* output an error message and return nothing. */
      ferror ("failed to build distance matrix union");
      return NULL;
    }

    /* reduce the tree list based on the new tree node. */
    if (!forest_union (forest, T, Ta, Tb)) {
      /* output an error message and return nothing. */
      ferror ("failed to build forest union");
      return NULL;
    }
  }

  /* free the forest structure. */
  list_free (forest);

  /* free the copied distance matrix. */
  matrix_free (Dop);

  /* free the copied pca structure. */
  pca_free (Pop);

  /* return the resultant tree. */
  return T;
}

/* pca_tree_free: frees a generated pca tree structure.
 * @T: the pca tree to free.
 */
void pca_tree_free (struct pca_tree *T) {
  /* don't free a null structure. */
  if (!T) return;

  /* free the node name. */
  free (T->name);
  T->name = NULL;

  /* free the node membership string. */
  free (T->membstr);
  T->membstr = NULL;

  /* free the child nodes. */
  pca_tree_free (T->child_a);
  pca_tree_free (T->child_b);

  /* free the current node. */
  free (T);
  T = NULL;
}

/* pca_tree_is_leaf: returns whether the given tree node is a leaf.
 * @T: the pca tree node to check.
 */
int pca_tree_is_leaf (struct pca_tree *T) {
  /* if the tree node has children, it's not a leaf. */
  return ((T->child_a || T->child_b) ? 0 : 1);
}

/* pca_tree_num_leaves: counts the number of leaf nodes in a tree.
 * @T: the pca tree structure to count leaves in.
 */
unsigned long pca_tree_num_leaves (struct pca_tree *T) {
  /* declare required variables. */
  unsigned long n;

  /* return nothing if the current node is undefined. */
  if (!T) return 0;

  /* initialize the counter at the current level. */
  n = pca_tree_is_leaf (T);

  /* add the number of leaves from each child node. */
  n += pca_tree_num_leaves (T->child_a);
  n += pca_tree_num_leaves (T->child_b);

  /* return the counted value. */
  return n;
}

/* pca_tree_bootstrap_compare_sub: runs a bootstrap comparison of two trees.
 * @T: the main tree to add bootstrap numbers to.
 * @Tboot: the bootstrap iteration tree to search.
 */
void pca_tree_bootstrap_compare_sub (struct pca_tree *T,
                                     struct pca_tree *Tboot) {
  /* stop running when we've hit a leaf. */
  if (!T || !Tboot || pca_tree_is_leaf (T) || pca_tree_is_leaf (Tboot))
    return;

  /* see if the current main tree node is found at this level of
   * the bootstrap tree with any given substructure.
   */
  if (strcmp (Tboot->membstr, T->membstr) == 0) {
    /* a match was found. increment the bootstrap number. */
    T->p += 1.0;
  }

  /* traverse the bootstrap tree child nodes. */
  pca_tree_bootstrap_compare_sub (T, Tboot->child_a);
  pca_tree_bootstrap_compare_sub (T, Tboot->child_b);
}

/* pca_tree_bootstrap_compare: runs a bootstrap comparison of two trees.
 * @T: the main tree to add bootstrap numbers to.
 * @Tboot: the bootstrap iteration tree to search.
 */
void pca_tree_bootstrap_compare (struct pca_tree *T, struct pca_tree *Tboot) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return;

  /* search the bootstrap tree for the current main tree node. */
  pca_tree_bootstrap_compare_sub (T, Tboot);

  /* traverse the child nodes. */
  pca_tree_bootstrap_compare (T->child_a, Tboot);
  pca_tree_bootstrap_compare (T->child_b, Tboot);
}

/* pca_tree_bootstrap_end: finishes a bootstrap run of a tree.
 * @T: the main tree to finalize bootstrap numbers of.
 * @n: the number of bootstrap iterations run.
 */
int pca_tree_bootstrap_end (struct pca_tree *T, unsigned long n) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return 1;

  /* divide the bootstrap count by the bootstrap number. */
  T->p = ceil (100.0 * T->p / ((double) n));

  /* traverse the child nodes. */
  return pca_tree_bootstrap_end (T->child_a, n) &&
         pca_tree_bootstrap_end (T->child_b, n);
}

/* pca_tree_recount_points: re-counts the points in inner tree nodes.
 * @T: the main tree to re-count points inside of.
 */
void pca_tree_recount_points (struct pca_tree *T) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return;

  /* traverse the tree first. */
  pca_tree_recount_points (T->child_a);
  pca_tree_recount_points (T->child_b);

  /* only change inner nodes. */
  if (!pca_tree_is_leaf (T->child_a)) {
    /* use the average of the left child node's child point counts. */
    T->na = (T->child_a->na + T->child_a->nb) / 2.0;
  }

  /* again, only change inner nodes. */
  if (!pca_tree_is_leaf (T->child_b)) {
    /* use the average of the right child node's child point counts. */
    T->nb = (T->child_b->na + T->child_b->nb) / 2.0;
  }
}

/* pca_tree_assign_pvalues: assigns probabilities to inner nodes of a tree.
 * @T: the main tree to calculate p-values inside of.
 * @D: the number of dimensions in the input scores.
 * @metric: the metric used to generate distances.
 */
int pca_tree_assign_pvalues (struct pca_tree *T,
                             unsigned long D,
                             pca_distance_metric metric) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return 1;

  /* store the p-value on the tree node. */
  if (metric == pca_group_distance_mahalanobis) {
    /* use the T2 -> F -> p-value method. */
    T->p = mahalanobis_pvalue (T->dist_a + T->dist_b, T->na, T->nb, D);
  }
  else if (metric == pca_group_distance_mixture) {
    /* use the indirect chi-square method. */
    T->p = mixture_pvalue (T->dist_a + T->dist_b);
  }
  else if (metric == pca_group_distance_chisq) {
    /* use the direct chi-square method. */
    T->p = chisq_pvalue (T->dist_a + T->dist_b);
  }
  else {
    /* output an error message and return failure. */
    ferror ("invalid distance metric");
    return 0;
  }

  /* traverse the child nodes. */
  return pca_tree_assign_pvalues (T->child_a, D, metric) &&
         pca_tree_assign_pvalues (T->child_b, D, metric);
}

/* pca_tree_enum_leaves: assigns indices to leaf nodes of a tree.
 * @T: the tree to operate on.
 * @idx: always run this function with this equal to NULL.
 */
void pca_tree_enum_leaves (struct pca_tree *T, unsigned long *idx) {
  /* stop running when we've hit a null. */
  if (!T) return;

  /* traverse the tree first. */
  pca_tree_enum_leaves (T->child_a, idx);
  pca_tree_enum_leaves (T->child_b, idx);

  /* is this node a leaf? */
  if (pca_tree_is_leaf (T)) {
    /* use the current value of idx and increment it for the next leaf. */
    T->idx = *idx;
    (*idx)++;
  }
}

/* pca_tree_min_length: returns the smallest distance in the tree.
 * @T: the tree to operate on.
 * @l: a pointer to the output value.
 */
void pca_tree_min_length (struct pca_tree *T, double *l) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return;

  /* compare the distance to the left child node length. */
  if (T->dist_a < *l) *l = T->dist_a;

  /* compare the distance to the right child node length. */
  if (T->dist_b < *l) *l = T->dist_b;

  /* traverse the child nodes. */
  pca_tree_min_length (T->child_a, l);
  pca_tree_min_length (T->child_b, l);
}

/* pca_tree_max_length: returns the largest distance in the tree.
 * @T: the tree to operate on.
 * @l: a pointer to the output value.
 */
void pca_tree_max_length (struct pca_tree *T, double *l) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return;

  /* compare the distance to the left child node length. */
  if (T->dist_a > *l) *l = T->dist_a;

  /* compare the distance to the right child node length. */
  if (T->dist_b > *l) *l = T->dist_b;

  /* traverse the child nodes. */
  pca_tree_max_length (T->child_a, l);
  pca_tree_max_length (T->child_b, l);
}

/* pca_tree_max_x: returns the largest x coordinate in the tree.
 * @T: the tree to operate on.
 * @x: a pointer to the output value.
 */
void pca_tree_max_x (struct pca_tree *T, double *x) {
  /* stop running when we've hit a null. */
  if (!T) return;

  /* compare the current node x-value. */
  if (T->x > *x) *x = T->x;

  /* traverse the child nodes. */
  pca_tree_max_x (T->child_a, x);
  pca_tree_max_x (T->child_b, x);
}

/* pca_tree_max_name: returns the longest leaf node name string.
 * @T: the tree to operate on.
 * @n: a pointer to the output value.
 */
void pca_tree_max_name (struct pca_tree *T, unsigned long *n) {
  /* stop running when we've hit a null. */
  if (!T) return;

  /* is this a leaf node? */
  if (pca_tree_is_leaf (T) && strlen (T->name) > *n) *n = strlen (T->name);

  /* traverse the child nodes. */
  pca_tree_max_name (T->child_a, n);
  pca_tree_max_name (T->child_b, n);
}

/* pca_tree_rescale: rescales the internode distances of a tree.
 * @T: the tree to modify.
 * @imin: the input minimum length.
 * @imax: the input maximum length.
 * @omin: the output minimum length.
 * @omax: the output maximum length.
 */
void pca_tree_rescale (struct pca_tree *T,
                       double imin, double imax,
                       double omin, double omax) {
  /* stop running when we've hit a leaf. */
  if (!T || pca_tree_is_leaf (T)) return;

  /* rescale the current node lengths. */
  T->dist_a = omin + ((omax - omin) / (imax - imin)) * (T->dist_a - imin);
  T->dist_b = omin + ((omax - omin) / (imax - imin)) * (T->dist_b - imin);

  /* traverse the child nodes. */
  pca_tree_rescale (T->child_a, imin, imax, omin, omax);
  pca_tree_rescale (T->child_b, imin, imax, omin, omax);
}

/* pca_tree_register_y: assigns y-axis coordinates to nodes.
 * @T: the tree to operate on.
 * @dy: the height of a node in the tree.
 */
void pca_tree_register_y (struct pca_tree *T, double dy) {
  /* stop running when we've hit a null. */
  if (!T) return;

  /* traverse the tree first. */
  pca_tree_register_y (T->child_a, dy);
  pca_tree_register_y (T->child_b, dy);

  /* is this node a leaf? */
  if (pca_tree_is_leaf (T)) {
    /* the node y-coordinate is based on its index. */
    T->y = ((double) T->idx) * 2.0 * dy + dy;

    /* the node height is zero. */
    T->dy = 0.0;
  }
  else {
    /* the node top y-coordinate is based on its children. */
    T->y = 0.5 * (T->child_a->y + T->child_b->y);

    /* the node height is also based on the kiddies. */
    T->dy = fabs (T->child_a->y - T->child_b->y);
  }
}

/* pca_tree_register_x: assigns x-axis coordinates to nodes.
 * @T: the tree to operate on.
 * @dx: the minimum x-axis position of a node.
 */
void pca_tree_register_x (struct pca_tree *T, double dx) {
  /* stop running when we've hit a null. */
  if (!T) return;

  /* are we at the top? */
  if (T->parent) {
    /* no. which child are we? */
    if (strcmp (T->name, T->parent->child_a->name) == 0) {
      /* first child. */
      T->x = T->parent->x + T->parent->dist_a;
    }
    else {
      /* second child. */
      T->x = T->parent->x + T->parent->dist_b;
    }
  }
  else {
    /* yes. align this node on the left. */
    T->x = dx;
  }

  /* traverse the tree. */
  pca_tree_register_x (T->child_a, dx);
  pca_tree_register_x (T->child_b, dx);
}

/* pca_tree_draw_r: traverses a tree and draws it to a postscript file.
 * @P: the pca structure.
 * @T: the tree to traverse.
 * @fh: the output file handle.
 */
void pca_tree_draw_r (struct pca *P, struct pca_tree *T, FILE *fh) {
  /* declare required variables. */
  char nstr[16];

  /* stop running when we've hit a null. */
  if (!T) return;

  /* are we on a leaf? */
  if (pca_tree_is_leaf (T)) {
    /* should we draw symbols? */
    if (opts_geti (OPTS_S_KEY)) {
      /* draw the symbol. */
      draw_symbol (fh, pca_group_index (P, T->name), P->n_groups, T->x, T->y);
    }

    /* move to the beginning of the node text. */
    draw_moveto (fh, T->x + 4.0, T->y - 4.0);

    /* print the node name. */
    draw_text (fh, 16,  T->name);
  }
  else {
    /* draw the vertical line of the node. */
    draw_line (fh,
      T->x, T->y - T->dy / 2.0,
      T->x, T->y + T->dy / 2.0);

    /* draw the horizontal line to the first child. */
    draw_line (fh,
      T->x, T->y - T->dy / 2.0,
      T->x + T->dist_a, T->y - T->dy / 2.0);

    /* draw the horizontal line to the second child. */
    draw_line (fh,
      T->x, T->y + T->dy / 2.0,
      T->x + T->dist_b, T->y + T->dy / 2.0);

    /* build the p-value/bootstrap number string. */
    if (floor (T->p) == ceil (T->p)) {
      /* bootstrap number => integer. */
      snprintf (nstr, 16, "%.0lf", T->p);
    }
    else {
      /* p-value => float. */
      if (T->p < 1.0e-2) {
        /* scientific notation. */
        snprintf (nstr, 16, "%.1le", T->p);
      }
      else {
        /* floating notation. */
        snprintf (nstr, 16, "%.2lf", T->p);
      }
    }

    /* print the number string. */
    draw_moveto (fh, T->x + 2.0, T->y - 2.0);
    draw_text (fh, 12, nstr);
  }

  /* traverse the tree. */
  pca_tree_draw_r (P, T->child_a, fh);
  pca_tree_draw_r (P, T->child_b, fh);
}

/* pca_tree_print_r: traverses a tree and prints it to a char buffer.
 * @T: the tree to traverse.
 * @buf: the character buffer.
 */
void pca_tree_print_r (struct pca_tree *T, char **buf) {
  /* declare required variables. */
  unsigned long i, i1, i2, j, di, dj;
  char nstr[16];

  /* stop running when we've hit a null. */
  if (!T) return;

  /* calculate the indices. */
  i = (unsigned long) T->y;
  j = (unsigned long) floor (T->x);

  /* are we on a leaf? */
  if (pca_tree_is_leaf (T)) {
    /* loop through the characters of the node name. */
    for (dj = 0; dj < strlen (T->name); dj++) {
      /* store the character in the buffer. */
      buf[i][j + dj] = T->name[dj];
    }
  }
  else {
    /* find the height extents. */
    i1 = i - (unsigned long) floor (T->dy / 2.0);
    i2 = i + (unsigned long) ceil (T->dy / 2.0);

    /* draw the vertical line of the node. */
    for (di = i1; di <= i2; di++) buf[di][j] = '|';

    /* use special characters for node corners. */
    buf[i1][j] = '+';
    buf[i2][j] = '+';

    /* draw the horizontal line to the first child. */
    for (dj = 1; dj < (unsigned long) ceil (T->dist_a); dj++)
      buf[i1][j + dj] = '-';

    /* draw the horizontal line to the second child. */
    for (dj = 1; dj < (unsigned long) ceil (T->dist_b); dj++)
      buf[i2][j + dj] = '-';

    /* build the p-value/bootstrap number string. */
    if (floor (T->p) == ceil (T->p)) {
      /* bootstrap number => integer. */
      snprintf (nstr, 16, "%.0lf", T->p);
    }
    else {
      /* p-value => float. */
      if (T->p < 1.0e-2) {
        /* scientific notation. */
        snprintf (nstr, 16, "%.1le", T->p);
      }
      else {
        /* floating notation. */
        snprintf (nstr, 16, "%.2lf", T->p);
      }
    }

    /* print the number string. */
    for (dj = 0; dj < strlen (nstr); dj++)
      buf[i][j + dj + 1] = nstr[dj];
  }

  /* traverse the tree. */
  pca_tree_print_r (T->child_a, buf);
  pca_tree_print_r (T->child_b, buf);
}

/* pca_tree_draw: draws a tree in a postscript document.
 * @P: the pca structure.
 * @T: the tree to draw.
 * @fname: the output filename.
 */
int pca_tree_draw (struct pca *P, struct pca_tree *T, const char *fname) {
  /* declare required variables. */
  double lmin = DBL_MAX, lmax = 0.0, xmax = 0.0;
  unsigned long n, cmax = 0, idx = 0;
  FILE *fh;

  /* get the number of leaf nodes. */
  n = pca_tree_num_leaves (T);

  /* calculate the node length bounds. */
  pca_tree_min_length (T, &lmin);
  pca_tree_max_length (T, &lmax);

  /* rescale the tree node lengths. */
  pca_tree_rescale (T, lmin, lmax, 20.0, 200.0);

  /* number the tree leaf nodes. */
  pca_tree_enum_leaves (T, &idx);

  /* build up the coordinates of the nodes. */
  pca_tree_register_y (T, 20.0);
  pca_tree_register_x (T, 20.0);

  /* get the extents of the tree nodes. */
  pca_tree_max_x (T, &xmax);
  pca_tree_max_name (T, &cmax);

  /* open the postscript document. */
  fh = draw_document_open (fname,
    xmax + 24.0 * (double) cmax,
    20.0 + 40.0 * (double) n);

  /* ensure we opened the document. */
  if (!fh) {
    /* output an error message and return failure. */
    ferror ("failed to open postscript document for tree");
    return 0;
  }

  /* draw the tree to the postscript document. */
  pca_tree_draw_r (P, T, fh);

  /* close the postscript document. */
  draw_document_close (fh);

  /* return success. */
  return 1;
}

/* pca_tree_print: prints a tree to standard output.
 * @T: the tree to print.
 */
int pca_tree_print (struct pca_tree *T) {
  /* declare required variables. */
  double lmin = DBL_MAX, lmax = 0.0, xmax = 0.0;
  unsigned long n, w, h, i, j, cmax = 0, idx = 0;
  char **buf;

  /* get the number of leaf nodes. */
  n = pca_tree_num_leaves (T);

  /* calculate the node length bounds. */
  pca_tree_min_length (T, &lmin);
  pca_tree_max_length (T, &lmax);

  /* rescale the tree node lengths. */
  pca_tree_rescale (T, lmin, lmax, 6.0, 60.0);

  /* number the tree leaf nodes. */
  pca_tree_enum_leaves (T, &idx);

  /* build up the coordinates of the nodes. */
  pca_tree_register_y (T, 1.0);
  pca_tree_register_x (T, 1.0);

  /* find the tree extents. */
  pca_tree_max_x (T, &xmax);
  pca_tree_max_name (T, &cmax);

  /* calculate the size of the buffer string array. */
  w = ((unsigned long) ceil (xmax)) + cmax + 2;
  h = (2 * n) + 1;

  /* allocate memory for the output buffer string array. */
  buf = (char **) malloc (sizeof (char *) * h);

  /* make sure we allocated memory. */
  if (!buf) {
    /* output an error message and return failure. */
    ferror ("failed to allocate memory for character buffer");
    return 0;
  }

  /* loop through the lines of the buffer. */
  for (i = 0; i < h; i++) {
    /* allocate memory for the current line. */
    buf[i] = (char *) malloc (sizeof (char) * w);

    /* make sure the line was allocated. */
    if (!buf[i]) {
      /* output an error message and return failure. */
      ferror ("failed to allocate memory for character buffer line %lu", i);
      return 0;
    }

    /* loop through the characters. */
    for (j = 0; j < w; j++) {
      /* initialize the character. */
      buf[i][j] = ' ';
    }
  }

  /* place the characters into the buffer. */
  pca_tree_print_r (T, buf);

  /* loop through the lines, printing to standard output. */
  for (i = 0; i < h; i++) {
    /* terminate the current line. */
    buf[i][w - 1] = '\0';

    /* print the current line. */
    fprintf (stdout, "%s\n", buf[i]);
    fflush (stdout);

    /* free the current line. */
    free (buf[i]);
    buf[i] = NULL;
  }

  /* free the string array. */
  free (buf);
  buf = NULL;

  /* return success. */
  return 1;
}

