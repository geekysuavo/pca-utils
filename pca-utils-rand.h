
/*
 * pca-utils-rand.h: header code for random number generation.
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

#ifndef __PCA_UTILS_RAND__
#define __PCA_UTILS_RAND__

/* begin function definitions below. */

void rand_init (unsigned long long seed);

unsigned long long randull (void);

double randf (void);

double randnormal (void);

#endif

