
/*
 * pca-utils-list.c: source code for linked lists.
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

/* list_new: allocate memory for a new linked list.
 * @cmpfn: function pointer to a comparison operator, or NULL
 *         to force pointer-pointer comparisons.
 */
struct list *list_new (int (*cmpfn) (void *d1, void *d2)) {
  /* declare required variables. */
  struct list *l;

  /* allocate a pointer for the list. */
  l = (struct list *) malloc (sizeof (struct list));

  /* make sure we allocated a pointer. */
  if (!l) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate new list pointer");
    return NULL;
  }

  /* set up the list internals. */
  l->cmpfn = cmpfn;
  l->prev = NULL;
  l->next = NULL;
  l->last = l;
  l->d = NULL;
  l->n = 0;

  /* return the list. */
  return l;
}

/* list_copy: duplicates a linked list; doesn't duplicate its data.
 * @la: the list to duplicate.
 */
struct list *list_copy (struct list *la) {
  /* declare required variables. */
  struct list *lb, *iter;

  /* allocate a new list. */
  lb = list_new (la->cmpfn);

  /* make sure we allocated a new list. */
  if (!lb) {
    /* output an error message and return nothing. */
    ferror ("failed to allocate list copy");
    return NULL;
  }

  /* loop through the links of the input list. */
  for (iter = la->next; iter; iter = iter->next) {
    /* append the input link data to the new list. */
    if (!list_append (lb, iter->d)) {
      /* output an error message and return nothing. */
      ferror ("failed to append element to list copy");
      return NULL;
    }
  }

  /* return the new list. */
  return lb;
}

/* list_free: frees up a linked list; doesn't free its data.
 * @l: the list to free.
 */
void list_free (struct list *l) {
  /* don't bother with null lists. */
  if (!l) return;

  /* recurse through the chain. */
  if (l->next)
    list_free (l->next);

  /* free the pointer memory. */
  free (l);
  l = NULL;
}

/* list_find: searches for data in a linked list. (returns -1 for no match)
 * @l: the list the search through.
 * @d: pointer to data to search for.
 */
signed long list_find (struct list *l, void *d) {
  /* declare required variables. */
  signed long count;
  struct list *iter;

  /* loop through the links of the list, keeping a counter. */
  for (count = 0, iter = l->next; iter; iter = iter->next, count++) {
    /* compare the current link data to the search data. */
    if ((l->cmpfn && l->cmpfn (iter->d, d) == 0) || iter->d == d) {
      /* the data match. return the counter. */
      return count;
    }
  }

  /* no match was found. return an indication of failure. */
  return -1;
}

/* list_append: appends data to a linked list.
 * @l: the list to append to.
 * @d: pointer to data to append.
 */
int list_append (struct list *l, void *d) {
  /* declare required variables. */
  struct list *ladd;

  /* make sure the input structures are defined. */
  if (!l || !d) {
    /* output an error message and return failure. */
    ferror ("input data structures to list append are invalid");
    return 0;
  }

  /* allocate a pointer to the new list link. */
  ladd = (struct list *) malloc (sizeof (struct list));

  /* make sure we allocated a new list link. */
  if (!ladd) {
    /* output an error message and return failure. */
    ferror ("failed to allocate new list link for append");
    return 0;
  }

  /* set up the new link pointers. */
  ladd->prev = l->last;
  ladd->next = NULL;
  ladd->last = NULL;
  ladd->d = d;

  /* add the new link into the chain. */
  l->last->next = ladd;
  l->last = ladd;
  l->n++;

  /* return success. */
  return 1;
}

/* list_append_safe: only appends unique data to a linked list.
 * @l: the list to append to.
 * @d: pointer to data to append.
 */
int list_append_safe (struct list *l, void *d) {
  /* make sure the data doesn't exist in the list. */
  if (list_find (l, d) >= 0) return 1;

  /* it doesn't. append the data. */
  return list_append (l, d);
}

/* list_prepend: prepends data to a linked list.
 * @l: the list to prepend to.
 * @d: pointer to data to prepend.
 */
int list_prepend (struct list *l, void *d) {
  /* declare required variables. */
  struct list *ladd;

  /* make sure the input structures are defined. */
  if (!l || !d) {
    /* output an error message and return failure. */
    ferror ("input data structures to list prepend are invalid");
    return 0;
  }

  /* prepending is the same as appending for empty lists. :) */
  if (!l->n) {
    /* append the data to the list. */
    return list_append (l, d);
  }

  /* allocate a pointer to the new list link. */
  ladd = (struct list *) malloc (sizeof (struct list));

  /* make sure we allocated a new list link. */
  if (!ladd) {
    /* output an error message and return failure. */
    ferror ("failed to allocate new list link for prepend");
    return 0;
  }

  /* set up the new link pointers. */
  ladd->prev = l;
  ladd->next = l->next;
  ladd->last = NULL;
  ladd->d = d;

  /* add the new link into the chain. */
  l->next->prev = ladd;
  l->next = ladd;
  l->n++;

  /* return success. */
  return 1;
}

/* list_prepend_safe: only prepends unique data to a linked list.
 * @l: the list to prepend to.
 * @d: pointer to data to prepend.
 */
int list_prepend_safe (struct list *l, void *d) {
  /* make sure the data doesn't exist in the list. */
  if (list_find (l, d) >= 0) return 1;

  /* it doesn't. prepend the data. */
  return list_prepend (l, d);
}

/* list_insert: inserts data into a linked list.
 * @l: the list to insert into.
 * @li: the list link to insert after.
 * @d: pointer to data to insert.
 */
int list_insert (struct list *l, struct list *li, void *d) {
  /* declare required variables. */
  struct list *lins, *lj;

  /* make sure the input data is defined. */
  if (!l || !li || !d) {
    /* output an error message and return failure. */
    ferror ("input data structures to list insert are invalid");
    return 0;
  }

  /* get the next link after the insert link. */
  lj = li->next;

  /* handle special linkage cases. */
  if (li == l) {
    /* the insert link is the list link => prepend. */
    return list_prepend (l, d);
  }
  else if (li == l->last) {
    /* the insert link is the last link => append. */
    return list_append (l, d);
  }

  /* allocate a pointer to the new list link. */
  lins = (struct list *) malloc (sizeof (struct list));

  /* make sure we allocated a new list link. */
  if (!lins) {
    /* output an error message and return failure. */
    ferror ("failed to allocate new list link for insert");
    return 0;
  }

  /* set up the new link pointers. */
  lins->prev = li;
  lins->next = lj;
  lins->last = NULL;
  lins->d = d;

  /* add the new link into the chain. */
  li->next = lins;
  lj->prev = lins;
  l->n++;

  /* return success. */
  return 1;
}

/* list_delete: deletes data from a linked list.
 * @l: the list to delete from.
 * @d: pointer to data to delete.
 */
int list_delete (struct list *l, void *d) {
  /* declare required variables. */
  struct list *iter;

  /* make sure the input data are defined. */
  if (!l || !d) {
    /* output an error message and return nothing. */
    ferror ("input data structures to list delete are invalid");
    return 0;
  }

  /* loop through the links of the list. */
  for (iter = l->next; iter; iter = iter->next) {
    /* compare the current link data to the delete data. */
    if ((l->cmpfn && l->cmpfn (iter->d, d) == 0) || iter->d == d) {
      /* the data match. pass over the link backwards. */
      if (iter->prev)
        iter->prev->next = iter->next;

      /* and pass over the link forwards. */
      if (iter->next)
        iter->next->prev = iter->prev;

      /* fix the last link pointer if this is the last link. */
      if (iter == l->last)
        l->last = iter->prev;

      /* free the current link. */
      free (iter);
      iter = NULL;
      l->n--;

      /* break from the loop. */
      break;
    }
  }

  /* return success. */
  return 1;
}

/* list_get: gets data from a linked list by its list index.
 * @l: the list to get from.
 * @i: the index of the data to get.
 */
void *list_get (struct list *l, unsigned long i) {
  /* declare required variables. */
  unsigned long count;
  struct list *iter;

  /* avoid out-of-bounds errors. */
  if (i >= l->n) {
    /* output an error message and return nothing. */
    ferror ("list index out of bounds (%lu >= %lu)", i, l->n);
    return NULL;
  }

  /* loop through the links of the list, keeping a counter. */
  for (count = 0, iter = l->next; iter; iter = iter->next, count++) {
    /* if the counter matches, return the current link data. */
    if (count == i) return iter->d;
  }

  /* no match found. */
  ferror ("unexpected error: this should not happen");
  return NULL;
}

/* list_foreach: perform the same action on every data element of a list.
 * @l: the list to perform actions on.
 * @fn: function pointer to the action to perform.
 */
void list_foreach (struct list *l, void (*fn) (void *d)) {
  /* declare required variables. */
  struct list *iter;

  /* loop through the links of the list. */
  for (iter = l->next; iter; iter = iter->next) {
    /* perform the action given by the function pointer. */
    if (iter->d) fn (iter->d);
  }
}

