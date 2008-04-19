/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*
*  Distributed under the terms of the GNU General Public License (GPL)
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "misc.h"

/***************************************************/
void die(char *errormessage) {
  /*This function prints the error message and raises SIGABRT.*/

  fprintf(stderr, "\a%s\n", errormessage);
  abort();
}

/***************************************************/

/* MEMLEAK, use free */
void *safeCalloc( int count, int size ) {
  /* this function calls calloc with the given inputs, 
     but dies with an error message if a NULL is returned */

  void *newthing=calloc( count, size );
  if (newthing==NULL) {
    die("calloc returned NULL");
    return NULL; /* unreachable. */
  }
  else return newthing;
}

/***************************************************/

/* MEMLEAK, use free */
void *safeMalloc( int count, int size ) {
  /* this function calls malloc with the inputs, which are
     to be provided in calloc notation. If the result is
     NULL, the program dies with an error message.*/

  void *newthing=malloc( count*size );
  if (newthing==NULL) {
    die("malloc returned NULL");
    return NULL; /* unreachable */
  }
  else return newthing;
}

/***************************************************/
/***************************************************/

BIT coinFlip() {
  if (rand() < RAND_MAX/2) {
    return 0;
  }  else {
    return 1;
  }
}
