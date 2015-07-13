
/*
 * pca-utils-getopt.c: source code for argument parsing.
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

/* define the option string used to parse command-line arguments. */
#define OPTS_OPSTR "i:o:m:n:D:1:2:3:u:v:r:L:Hkh"

/* define the array of option structures for parsing. */
static struct option long_options[] = {
  { OPTS_L_INPUT,  required_argument, 0, OPTS_S_INPUT },
  { OPTS_L_OUTPUT, required_argument, 0, OPTS_S_OUTPUT },
  { OPTS_L_METRIC, required_argument, 0, OPTS_S_METRIC },
  { OPTS_L_COUNT,  required_argument, 0, OPTS_S_COUNT },
  { OPTS_L_DIVS,   required_argument, 0, OPTS_S_DIVS },
  { OPTS_L_PC1,    required_argument, 0, OPTS_S_PC1 },
  { OPTS_L_PC2,    required_argument, 0, OPTS_S_PC2 },
  { OPTS_L_PC3,    required_argument, 0, OPTS_S_PC3 },
  { OPTS_L_MEAN,   required_argument, 0, OPTS_S_MEAN },
  { OPTS_L_VAR,    required_argument, 0, OPTS_S_VAR },
  { OPTS_L_ROT,    required_argument, 0, OPTS_S_ROT },
  { OPTS_L_LABEL,  required_argument, 0, OPTS_S_LABEL },
  { OPTS_L_HEADER, no_argument,       0, OPTS_S_HEADER },
  { OPTS_L_KEY,    no_argument,       0, OPTS_S_KEY },
  { OPTS_L_HELP,   no_argument,       0, OPTS_S_HELP },
  { 0, 0, 0, 0 }
};

/* declare the locally stored command-line arguments. */
char **opts_argv;
int opts_argc;

/* opts_set: stores an option argument at its index in a local array.
 * @opt: the option char index.
 * @arg: the option argument.
 */
int opts_set (const char opt, const char *arg) {
  /* declare required variables. */
  unsigned int opti;

  /* get the integer cast of the option char. */
  opti = (unsigned int) opt;

  /* ensure we successfully cast to an integer. */
  if (!opti) {
    /* print a failure message and exit the function. */
    ferror ("failed to cast option char index '%c' to integer", opt);
    return 0;
  }

  /* check if the arguments array needs to be resized. */
  if (opti >= opts_argc) {
    /* make just enough room. i know, it's wasteful. */
    opts_argc = opti + 1;

    /* reallocate memory for the option arguments array. */
    opts_argv = (char **) realloc (opts_argv, sizeof (char *) * opts_argc);

    /* ensure we successfully allocated new memory. */
    if (!opts_argv) {
      /* print a failure message and exit the function. */
      ferror ("failed to reallocate option argument string memory");
      return 0;
    }
  }

  /* allocate memory for the argument string. */
  opts_argv[opti] = (char *) malloc (sizeof (char) * (strlen (arg) + 2));

  /* ensure we successfully allocated string memory. */
  if (!opts_argv[opti]) {
    /* print a failure message and exit the function. */
    ferror ("failed to allocate option argument string memory");
    return 0;
  }

  /* copy the argument string into the array. */
  strcpy (opts_argv[opti], arg);

  /* return success. */
  return 1;
}

/* opts_init: initializes the local argument array for later use.
 * @argc: the number of arguments passed to main().
 * @argv: the string array containing arguments from main().
 */
int opts_init (int argc, char **argv) {
  /* declare required variables. */
  int idx, optidx;

  /* allocate an initial small pointer for the arguments string array. */
  opts_argv = (char **) malloc (sizeof (char *));
  opts_argv[0] = NULL;
  opts_argc = 1;

  /* loop until we don't have any more argument to parse. */
  while (1) {
    /* initialize the option index to zero. */
    optidx = 0;

    /* parse a new option from the argv array. */
    idx = getopt_long (argc, argv, OPTS_OPSTR, long_options, &optidx);

    /* break when done or on any errors. */
    if (idx == -1) break;

    /* determine which option was passed to the program. */
    switch (idx) {
      /* required-argument options: store the argument. */
      case OPTS_S_INPUT:
      case OPTS_S_OUTPUT:
      case OPTS_S_METRIC:
      case OPTS_S_COUNT:
      case OPTS_S_DIVS:
      case OPTS_S_PC1:
      case OPTS_S_PC2:
      case OPTS_S_PC3:
      case OPTS_S_MEAN:
      case OPTS_S_VAR:
      case OPTS_S_ROT:
      case OPTS_S_LABEL:
        /* attempt to store the argument in the array. */
        if (!opts_set (idx, optarg)) {
          /* print an error message. */
          ferror ("failed to store passed argument '%s' to '%c'",
                  optarg, (char) idx);

          /* return failure. */
          return 0;
        }

        /* don't fall through. */
        break;

      /* non-argument options: store a default argument. */
      case OPTS_S_HEADER:
      case OPTS_S_KEY:
      case OPTS_S_HELP:
        /* set a default numerically parseable TRUE argument. */
        if (!opts_set (idx, "1")) {
          /* print an error message and return failure. */
          ferror ("failed to store default argument to '%c'", (char) idx);
          return 0;
        }

        /* don't fall through. */
        break;

      /* unacceptable options: fail out. */
      default:
        /* print a failure message and return failure. */
        ferror ("unable to parse command-line arguments");
        return 0;
    }
  }

  /* return success. */
  return 1;
}

/* opts_gets: retrieve a stored argument string from the local array.
 * @opt: the option char index.
 */
char *opts_gets (const char opt) {
  /* declare required variables. */
  unsigned int opti = (unsigned int) opt;

  /* return the argument string, or NULL if its undefined. */
  return (opti && opti < opts_argc ? opts_argv[opti] : NULL);
}

/* opts_geti: retrieve a stored argument integer from the local array.
 * @opt: the option char index.
 */
long opts_geti (const char opt) {
  /* declare required variables. */
  char *str = opts_gets (opt);

  /* return the translated integer, or zero if its undefined. */
  return (str ? atol (str) : 0);
}

/* opts_getf: retrieve a stored argument float from the local array.
 * @opt: the option char index.
 */
double opts_getf (const char opt) {
  /* declare required variables. */
  char *str = opts_gets (opt);

  /* return the translated double, or zero if its undefined. */
  return (str ? atof (str) : 0.0);
}

