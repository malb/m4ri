/**
 * Dan Bernstein's "Optimizing linear maps mod 2"
 *
 * This code is a port of Dan's sort1.cpp which is available at
 * http://binary.cr.yp.to/linearmod2.html
 *
 * It makes use of a binary heap written by Martin Kunev which is available at
 * https://gist.github.com/martinkunev/1365481
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RI_DJB_H
#define M4RI_DJB_H

#include <m4ri/mzd.h>

typedef enum {
  source_target,
  source_source
} srctyp_t;
 

typedef struct {
  rci_t nrows;
  rci_t ncols;
  rci_t *target;
  rci_t *source;
  srctyp_t *srctyp;
  rci_t length;
  wi_t allocated;
} djb_t;

#define M4RI_DJB_BASE_SIZE 64

static inline djb_t *djb_init(rci_t nrows, rci_t ncols) {
  /* we want to use realloc, so we call unaligned malloc */
  djb_t *m = malloc(sizeof(djb_t));
  if (m == NULL)
    m4ri_die("malloc failed.\n");

  *m = (djb_t){
    .nrows = nrows,
    .ncols = ncols,
    .target = malloc(sizeof(rci_t) * M4RI_DJB_BASE_SIZE),
    .source = malloc(sizeof(rci_t) * M4RI_DJB_BASE_SIZE),
    .srctyp = malloc(sizeof(srctyp_t) *   M4RI_DJB_BASE_SIZE),
    .length = 0,
    .allocated = M4RI_DJB_BASE_SIZE,
  };
  if (m->target == NULL || m->source == NULL || m->srctyp == NULL)
    m4ri_die("malloc failed.\n");
  return m;
}

static inline void djb_free(djb_t *m) {
  free(m->target);
  free(m->source);
  free(m->srctyp);
  free(m);
}

static inline void djb_push_back(djb_t *h, rci_t target, rci_t source, int srctyp) {
  assert((target <= h->nrows) && 
         ((source < h->ncols) | (srctyp != source_source)) &&
         ((source < h->nrows) | (srctyp != source_target)));
  if (h->length >= h->allocated) {
    h->allocated += M4RI_DJB_BASE_SIZE;
    h->target = realloc(h->target, h->allocated*sizeof(rci_t));
    h->source = realloc(h->source, h->allocated*sizeof(rci_t));
    h->srctyp = realloc(h->srctyp, h->allocated*sizeof(srctyp_t));
  }
  h->target[h->length] = target;
  h->source[h->length] = source;
  h->srctyp[h->length] = srctyp;
  h->length++;
}

/**
 * Compile a new DJB linear map data structure from A.
 *
 * \param A
 */

djb_t *djb_compile(mzd_t *A);

/**
 * W = m*V
 *
 * Apply the linear map m to V and write the result in W.
 *
 * \param m Linear map.
 */

void djb_apply_mzd(djb_t *m, mzd_t *W, const mzd_t *V);


/**
 *
 */

static inline void djb_info(djb_t *h) {
  double save = (double)h->length / (double)(h->nrows * h->ncols);
  printf("%d x %d linear map in %d xors (cost: %.5f)\n",h->nrows, h->ncols, h->length, save);
}


#endif //M4RI_DJB_H
