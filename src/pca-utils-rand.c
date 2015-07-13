
/*
 * pca-utils-rand.c: source code for random number generation.
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

/* global random variable storage values, used on each iteration of 
 * randf(), randull(), etc.
 */
unsigned long long rand_u, rand_v, rand_w;

/* rand_init: initializes the random number generator.
 * @seed: the seed value, usuall the time.
 */
void rand_init (unsigned long long seed) {
  /* initialize based on the seed. */
  rand_v = 4101842887655102017LL;
  rand_u = seed ^ rand_v; randull ();
  rand_v = rand_u; randull ();
  rand_w = rand_v; randull ();
}

/* randull: return a uniformly distributed random deviate inside an
 * unsigned long long datatype.
 */
unsigned long long randull (void) {
  rand_u = rand_u * 2862933555777941757LL + 7046029254386353087LL;
  rand_v ^= rand_v >> 17;
  rand_v ^= rand_v << 31;
  rand_v ^= rand_v >> 8;
  rand_w = 4294957665U * (rand_w & 0xffffffff) + (rand_w >> 32);
  unsigned long long x = rand_u ^ (rand_u << 21);
  x ^= x >> 35;
  x ^= x << 4;
  return (x + rand_v) ^ rand_w;
}

/* randf: return a uniformly distributed random deviate inside a
 * double datatype. range is [0,1].
 */
double randf (void) {
  return 5.42101086242752217E-20 * randull ();
}

/* randnormal: returns a unit-variance, zero-mean normal deviate.
 */
double randnormal (void) {
  double x1, x2, w;
  do {
    x1 = 2.0 * randf () - 1.0;
    x2 = 2.0 * randf () - 1.0;
    w = x1 * x1 + x2 * x2;
  } while (w >= 1.0);
  w = sqrt (-2.0 * log (w) / w);
  return x1 * w / sqrt (2.0 * M_PI);
}

