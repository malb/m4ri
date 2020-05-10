/*
 * \author Martin Kunev <martinkunev@gmail.com> original implementation of C99 heap
 * \author Martin Albrecht <martinralbrecht@googlemail.com> adapted to M4RI + DJB map
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <m4ri/djb.h>
#include <m4ri/io.h>
#include <m4ri/misc.h>
#include <m4ri/mzd.h>
#include <m4ri/xor.h>

#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <stdlib.h>

static inline int mzd_compare_rows_revlex(const mzd_t *A, rci_t a, rci_t b) {
  for (wi_t j = A->width - 1; j >= 0; j--) {
    if (A->rows[a][j] < A->rows[b][j]) return 0;
    if (A->rows[a][j] > A->rows[b][j]) return 1;
  }
  return 1;
}

/**
 * \brief A Heap
 */

typedef struct heap {
  unsigned int size;  /*!< Size of the allocated memory (in number of items) */
  unsigned int count; /*!<  Count of the elements in the heap */
  rci_t *data;        /*!< Array with the elements */
} heap_t;

// Returns the biggest element in the heap
#define heap_front(h) (*(h)->data)

void heap_free(heap_t *h) {
  free(h->data);
  free(h);
}

static const unsigned int heap_base_size = 4;

// Prepares the heap for use
heap_t *heap_init() {
  heap_t *h = (heap_t *)malloc(sizeof(heap_t));
  if (h == NULL) m4ri_die("malloc failed.\n");
  *h = (heap_t){.size = heap_base_size, .count = 0, .data = malloc(sizeof(rci_t) * heap_base_size)};
  if (h->data == NULL) m4ri_die("malloc failed.\n");
  return h;
}

// Inserts element to the heap
void heap_push(struct heap *restrict h, rci_t value, const mzd_t *A) {
  unsigned int index, parent;

  // Resize the heap if it is too small to hold all the data
  if (h->count == h->size) {
    h->size <<= 1;
    h->data = realloc(h->data, sizeof(rci_t) * h->size);
    if (h->data == NULL) m4ri_die("realloc failed.\n");
  }

  // Find out where to put the element and put it
  for (index = h->count++; index; index = parent) {
    parent = (index - 1) >> 1;
    if (mzd_compare_rows_revlex(A, h->data[parent], value)) break;
    h->data[index] = h->data[parent];
  }
  h->data[index] = value;
}

// Removes the biggest element from the heap
void heap_pop(struct heap *restrict h, const mzd_t *A) {
  unsigned int index, swap, other;

  // Remove the biggest element
  rci_t temp = h->data[--h->count];

  // Resize the heap if it's consuming too much memory
  if ((h->count <= (h->size >> 2)) && (h->size > heap_base_size)) {
    h->size >>= 1;
    h->data = realloc(h->data, sizeof(rci_t) * h->size);
    if (h->data == NULL) m4ri_die("realloc failed.\n");
  }

  // Reorder the elements
  for (index = 0; 1; index = swap) {
    // Find the child to swap with
    swap = (index << 1) + 1;
    if (swap >= h->count) break;  // If there are no children, the heap is reordered
    other = swap + 1;
    if ((other < h->count) && mzd_compare_rows_revlex(A, h->data[other], h->data[swap]))
      swap = other;
    if (mzd_compare_rows_revlex(A, temp, h->data[swap]))
      break;  // If the bigger child is less than or equal to its parent, the heap is reordered

    h->data[index] = h->data[swap];
  }
  h->data[index] = temp;
}

djb_t *djb_compile(mzd_t *A) {
  heap_t *h = heap_init();
  rci_t m   = A->nrows;
  rci_t n   = A->ncols;

  djb_t *z = djb_init(m, n);

  for (rci_t i = 0; i < m; i++) heap_push(h, i, A);  // sort by mzd_compare_rows_revlex

  while (n > 0) {
    if (mzd_read_bit(A, heap_front(h), n - 1) == 0) {
      --n;
      continue;
    }

    rci_t temp = heap_front(h);
    heap_pop(h, A);

    if (m >= 2 && mzd_read_bit(A, heap_front(h), n - 1)) {
      mzd_row_add(A, heap_front(h), temp);
      djb_push_back(z, temp, heap_front(h), source_target);
    } else {
      mzd_write_bit(A, temp, n - 1, 0);
      djb_push_back(z, temp, n - 1, source_source);
    }
    heap_push(h, temp, A);
  }
  heap_free(h);

  return z;
}

void djb_apply_mzd(djb_t *m, mzd_t *W, const mzd_t *V) {
  assert(W->width == V->width);
  rci_t i = m->length;
  while (i > 0) {
    --i;
    if (m->srctyp[i] == source_source) {
      _mzd_combine(mzd_row(W, m->target[i]), mzd_row(V, m->source[i]), W->width);
    } else {
      _mzd_combine(mzd_row(W, m->target[i]), mzd_row(W, m->source[i]), W->width);
    }
  }
}

void djb_print(djb_t *m) {
  int *iszero = m4ri_mm_malloc(sizeof(int) * m->nrows);
  for (int i = 0; i < m->nrows; ++i) iszero[i] = 1;

  rci_t i = m->length;
  while (i > 0) {
    --i;
    if (iszero[m->target[i]]) {
      if (m->srctyp[i] == source_source) {
        printf("cpy src[%d] to dst[%d]\n", m->source[i], m->target[i]);
      } else {
        printf("cpy dst[%d] to dst[%d]\n", m->source[i], m->target[i]);
      }
      iszero[m->target[i]] = 0;
    } else {
      if (m->srctyp[i] == source_source) {
        printf("add src[%d] to dst[%d]\n", m->source[i], m->target[i]);
      } else {
        printf("add dst[%d] to dst[%d]\n", m->source[i], m->target[i]);
      }
    }
  }
  m4ri_mm_free(iszero);
}
