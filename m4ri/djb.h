/**
 * \file djb.h
 *
 * \brief Dan Bernstein's "Optimizing linear maps mod 2"
 *
 * This code is a port of sort1.cpp available at http://binary.cr.yp.to/linearmod2.html
 *
 * Given a matrix A djb_compile(A) will compute a djb_t data structure which realises A with
 * (heuristically) (m * n)/(log m - loglog m) XORs.
 *
 * It makes use of a binary heap written by Martin Kunev which is available at
 * https://gist.github.com/martinkunev/1365481
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RI_DJB_H
#define M4RI_DJB_H

#include <m4ri/mzd.h>

/**
 * \brief Specify source type of addition
 */

typedef enum {
  source_target,  //< add from target matrix
  source_source   //< add from source matrix
} srctyp_t;

/**
 * \brief DJB's optimized linear maps mod 2
 */

typedef struct {
  rci_t nrows;      /*!< Number of rows of map */
  rci_t ncols;      /*!< Number of columns of map */
  rci_t *target;    /*!< target row at index i */
  rci_t *source;    /*!< source row at index i */
  srctyp_t *srctyp; /*!< source type at index i */
  rci_t length;     /*!< length of target, source and srctype */
  wi_t allocated;   /*!< how much did we allocate already */
} djb_t;

/**
 * Standard allocation chunk
 */

#define M4RI_DJB_BASE_SIZE 64

/**
 * Allocate a new DJB linear map
 *
 * \param nrows Number of rows
 * \param ncols Number of columns
 */

static inline djb_t *djb_init(rci_t nrows, rci_t ncols) {
  /* we want to use realloc, so we call unaligned malloc */
  djb_t *m = (djb_t *)malloc(sizeof(djb_t));
  if (m == NULL) m4ri_die("malloc failed.\n");

  m->nrows     = nrows;
  m->ncols     = ncols;
  m->target    = (rci_t *)malloc(sizeof(rci_t) * M4RI_DJB_BASE_SIZE);
  m->source    = (rci_t *)malloc(sizeof(rci_t) * M4RI_DJB_BASE_SIZE);
  m->srctyp    = (srctyp_t *)malloc(sizeof(srctyp_t) * M4RI_DJB_BASE_SIZE);
  m->length    = 0;
  m->allocated = M4RI_DJB_BASE_SIZE;

  if (m->target == NULL || m->source == NULL || m->srctyp == NULL) m4ri_die("malloc failed.\n");
  return m;
}

/**
 * Free a DJB linear maps
 *
 * \param m Map
 */

static inline void djb_free(djb_t *m) {
  free(m->target);
  free(m->source);
  free(m->srctyp);
  free(m);
}

/**
 * Add a new operation out[target] ^= srctype[source] to queue.
 *
 * \param z DJB linear map.
 * \param target Output index
 * \param source Input index
 * \param srctyp Type of input (source_source or source_target)
 */

static inline void djb_push_back(djb_t *z, rci_t target, rci_t source, srctyp_t srctyp) {
  assert((target < z->nrows) && ((source < z->ncols) | (srctyp != source_source)) &&
         ((source < z->nrows) | (srctyp != source_target)));
  if (z->length >= z->allocated) {
    z->allocated += M4RI_DJB_BASE_SIZE;
    z->target = (rci_t *)realloc(z->target, z->allocated * sizeof(rci_t));
    z->source = (rci_t *)realloc(z->source, z->allocated * sizeof(rci_t));
    z->srctyp = (srctyp_t *)realloc(z->srctyp, z->allocated * sizeof(srctyp_t));
  }
  z->target[z->length] = target;
  z->source[z->length] = source;
  z->srctyp[z->length] = srctyp;
  z->length++;
}

/**
 * Compile a new DJB linear map from A.
 *
 * \param A
 */

djb_t *djb_compile(mzd_t *A);

/**
 * \brief W = m*V
 *
 * Apply the linear map m to V and write the result in W.
 *
 * \param z  DJB linear map.
 * \param W  Output matrix
 * \param V  Input matrix
 */

void djb_apply_mzd(djb_t *z, mzd_t *W, const mzd_t *V);

/**
 * Print infomrmation on linear map mA
 */

static inline void djb_info(const djb_t *z) {
  double save = (double)z->length / (double)(z->nrows * z->ncols);
  printf("%d x %d linear map in %d xors (cost: %.5f)\n", z->nrows, z->ncols, z->length, save);
}

#endif  // M4RI_DJB_H
