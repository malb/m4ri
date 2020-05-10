/******************************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
 *    Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  version 2 or higher.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "graycode.h"
#include "misc.h"
#include <stdio.h>

code **m4ri_codebook = NULL;

int m4ri_gray_code(int number, int length) {
  int lastbit = 0;
  int res     = 0;
  for (int i = length - 1; i >= 0; --i) {
    int bit = number & (1 << i);
    res |= (lastbit >> 1) ^ bit;
    lastbit = bit;
  }
  return res;
}

void m4ri_build_code(int *ord, int *inc, int l) {
  for (int i = 0; i < (int)__M4RI_TWOPOW(l); ++i) { ord[i] = m4ri_gray_code(i, l); }

  for (int i = l; i > 0; --i) {
    for (int j = 1; j < (int)__M4RI_TWOPOW(i) + 1; ++j) {
      inc[j * __M4RI_TWOPOW(l - i) - 1] = l - i;
    }
  }
}

void m4ri_build_all_codes() {
  if (m4ri_codebook) { return; }
  m4ri_codebook = (code **)m4ri_mm_calloc(__M4RI_MAXKAY + 1, sizeof(code *));

  for (int k = 1; k < __M4RI_MAXKAY + 1; ++k) {
    m4ri_codebook[k]      = (code *)m4ri_mm_calloc(1, sizeof(code));
    m4ri_codebook[k]->ord = (int *)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(int));
    m4ri_codebook[k]->inc = (int *)m4ri_mm_calloc(__M4RI_TWOPOW(k), sizeof(int));
    m4ri_build_code(m4ri_codebook[k]->ord, m4ri_codebook[k]->inc, k);
  }
}

void m4ri_destroy_all_codes() {
  if (!m4ri_codebook) { return; }
  for (int i = 1; i < __M4RI_MAXKAY + 1; ++i) {
    m4ri_mm_free(m4ri_codebook[i]->inc);
    m4ri_mm_free(m4ri_codebook[i]->ord);
    m4ri_mm_free(m4ri_codebook[i]);
  }
  m4ri_mm_free(m4ri_codebook);
  m4ri_codebook = NULL;
}

int m4ri_opt_k(int a, int b, int c) {
  int n   = MIN(a, b);
  int res = MIN(__M4RI_MAXKAY, MAX(1, (int)(0.75 * (1 + log2_floor(n)))));
  return res;
}
