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
#ifdef HAVE_SSE2
#include <mm_malloc.h>
#endif

void m4ri_die(char *errormessage) {
  /*This function prints the error message and raises SIGABRT.*/

  fprintf(stderr, "\a%s\n", errormessage);
  abort();
}

void m4ri_print_bit_string(int number, int length){
  int i;
  for(i=length-1 ; i>=0; i--) {
    ((1<<i) & number) ? printf("1") : printf("0");
  }
  printf("\n");
}

/* Warning: I assume *destination has RADIX+1 bytes available */
void m4ri_word_to_str( char *destination, word data, int colon) {
  int i;
  int j = 0;

  if (colon == 0) {

    for (i=0; i<RADIX; i++) {
      if (GET_BIT(data,i))
	destination[i]='1';
      else 
	destination[i]='0';
    }
    destination[RADIX]='\0';

  } else {

    for (i=0; i<RADIX; i++) {
      if (GET_BIT(data,i))
	destination[j]='1';
      else 
	destination[j]='0';
      j++;
      if (((i % 4)==3) && (i!=RADIX-1)) {
	destination[j]=':';
	j++;
      }
    }

    destination[j]='\0';
  }
}

void *m4ri_mm_calloc( int count, int size ) {
  /* this function calls calloc with the given inputs, 
     but dies with an error message if a NULL is returned */

#ifdef HAVE_SSE2
  void *newthing = _mm_malloc(count*size, 16);
#else
  void *newthing = calloc(count, size);
#endif
  if (newthing==NULL) {
    m4ri_die("calloc returned NULL");
    return NULL; /* unreachable. */
  }
#ifdef HAVE_SSE2
  char *b = (char*)newthing;
  int i;
  for(i=0; i< count*size; i++) {
    b[i] = 0;
  }
#endif
  return newthing;
}

void *m4ri_mm_malloc( int size ) {
#ifdef HAVE_SSE2
  void *newthing = _mm_malloc(size, 16);
#else
  void *newthing=malloc( size );
#endif  
  if (newthing==NULL) {
    m4ri_die("malloc returned NULL");
    return NULL; /* unreachable */
  }
  else return newthing;
}

BIT m4ri_coin_flip() {
  if (rand() < RAND_MAX/2) {
    return 0;
  }  else {
    return 1;
  }
}
