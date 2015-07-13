
/*
 * pca-utils-getopt.h: header code for argument parsing.
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

#ifndef __PCA_UTILS_GETOPT__
#define __PCA_UTILS_GETOPT__

/* define the command-line option strings and characters. */
#define OPTS_L_INPUT  "input"
#define OPTS_S_INPUT  'i'
#define OPTS_L_OUTPUT "output"
#define OPTS_S_OUTPUT 'o'
#define OPTS_L_METRIC "metric"
#define OPTS_S_METRIC 'm'
#define OPTS_L_COUNT  "count"
#define OPTS_S_COUNT  'n'
#define OPTS_L_DIVS   "divisions"
#define OPTS_S_DIVS   'D'
#define OPTS_L_PC1    "component1"
#define OPTS_S_PC1    '1'
#define OPTS_L_PC2    "component2"
#define OPTS_S_PC2    '2'
#define OPTS_L_PC3    "component3"
#define OPTS_S_PC3    '3'
#define OPTS_L_MEAN   "mean"
#define OPTS_S_MEAN   'u'
#define OPTS_L_VAR    "var"
#define OPTS_S_VAR    'v'
#define OPTS_L_ROT    "rotation"
#define OPTS_S_ROT    'r'
#define OPTS_L_LABEL  "label"
#define OPTS_S_LABEL  'L'
#define OPTS_L_HEADER "header"
#define OPTS_S_HEADER 'H'
#define OPTS_L_KEY    "key"
#define OPTS_S_KEY    'k'
#define OPTS_L_HELP   "help"
#define OPTS_S_HELP   'h'

/* begin function definitions below. */

int opts_init (int argc, char **argv);

char *opts_gets (const char opt);

long opts_geti (const char opt);

double opts_getf (const char opt);

#endif

