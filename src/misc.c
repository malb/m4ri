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

word packingmask[RADIX];
word bytemask[RADIX/8];
word sixteenmask[RADIX/16];

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

void m4ri_setup_packing_masks() {
  int i, j;
  word x=1;

  for (i=RADIX-1; i>=0; i--) {
    packingmask[i]=x;
    x<<=1;
  }

  for (i=0; i<RADIX/8; i++) {
    x=0;
    for (j=0; j<8; j++) {
      x|=packingmask[j+i*8];
    }
    bytemask[i]=x; 
  }

  for (i=0; i<RADIX/16; i++) {
    x=0;
    for (j=0; j<16; j++) {
      x|=packingmask[j+i*16];
    }
    sixteenmask[i]=x; 
  }
} 
