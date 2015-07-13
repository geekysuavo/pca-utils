
/*
 * pca-utils-list.h: header code for linked lists.
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

#ifndef __PCA_UTILS_LIST__
#define __PCA_UTILS_LIST__

/* include the pca utils header. */
#include "pca-utils.h"

/* begin function definitions below. */

struct list *list_new (int (*cmpfn) (void *d1, void *d2));

struct list *list_copy (struct list *la);

void list_free (struct list *l);

signed long list_find (struct list *l, void *d);

int list_append (struct list *l, void *d);

int list_append_safe (struct list *l, void *d);

int list_prepend (struct list *l, void *d);

int list_prepend_safe (struct list *l, void *d);

int list_insert (struct list *l, struct list *li, void *d);

int list_delete (struct list *l, void *d);

void *list_get (struct list *l, unsigned long i);

void list_foreach (struct list *l, void (*fn) (void *d));

#endif

