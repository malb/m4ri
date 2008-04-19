#ifndef GRAYFLEX_H
#define GRAYFLEX_H

/**
 * gray code generation used by the M4RI algorithm.
 * 
 * @author Martin Albrecht 
 */

/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de> 
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
 
#include "misc.h"

struct codestruct {
  int *ord;
  int *inc;
};

typedef struct codestruct /*renamed as*/ code;

extern code **codebook;

#define MAXKAY 16
#define TWOPOW(i) (1<<(i))

void m2t_print_bit_string(int number, int length);

/**
 * swaps length bits in v naively.
 *
 * WARNING: Uppper bits of return value may contain garbage after
 * operation.
 */
int m2t_swap_bits(int v,int length);

/**
 * Returns the 'number'-th gray code entry for a gray code of length
 * $2^{length}$.
 * 
 * INPUT:
 *     number -- index in the gray code table
 *     length -- length of the gray code
 *
 * OUTPUT:
 *      number-th gray code entry
 *
 * AUTHOR: malb
 *
 * THANKS: Soroosh Yazdani explained the repeated sum idea to me.
 * 
 */

int m2t_gray_code(int number, int length);

/**
 * Fills in 'ord' and 'inc' with gray code data for a gray code of
 * length $2^{length}$.
 *
 * 
 * @param ord will hold gray code data, must be preallocated with correct size
 * @param inc will hold some increment data, must be preallocated with correct size
 *
 * @author Martin Albrecht
 *
 * Robert Miller had the idea for a non-recursive implementation.
 *
 */

void m2t_build_code(int *ord, int *inc, int length);

void m2t_build_all_codes();

void m2t_destroy_all_codes();

/**
 * Return the optimal k for the given parameters. If c != 0 then the
 * optimal k for multiplication is returned, else the optimal k for
 * inversion is returned. The optiomal $k$ here means $0.75 log_2(n)$
 * where $n$ is $min(a,b)$ for inversiona and $b$ for multiplication.
 */

int m2t_opt_k(int a,int b,int c);

#endif //GRAYFLEX_H
