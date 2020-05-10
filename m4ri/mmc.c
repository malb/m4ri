/*******************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2007, 2008 Gregory Bard <bard@fordham.edu>
 *    Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
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
 *
 ********************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mmc.h"

#if __M4RI_ENABLE_MMC
/**
 * The actual memory block cache.
 */

mmb_t m4ri_mmc_cache[__M4RI_MMC_NBLOCKS];
#endif  // __M4RI_ENABLE_MMC

/**
 * \brief Allocate size bytes.
 *
 * \param size Number of bytes.
 *
 * \return pointer to allocated memory block.
 */

void *m4ri_mmc_malloc(size_t size) {

#if __M4RI_ENABLE_MMC
  void *ret = NULL;

#if __M4RI_HAVE_OPENMP
#pragma omp critical(mmc)
  {
#endif
    mmb_t *mm = m4ri_mmc_cache;
    if (size <= __M4RI_MMC_THRESHOLD) {
      for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
        if (mm[i].size == size) {
          ret        = mm[i].data;
          mm[i].data = NULL;
          mm[i].size = 0;
          break;
        }
      }
    }
#if __M4RI_HAVE_OPENMP
  }
#endif
  if (ret)
    return ret;
  else
    return m4ri_mm_malloc(size);

#else  // __M4RI_ENABLE_MMC

  return m4ri_mm_malloc(size);

#endif  // __M4RI_ENABLE_MMC
}

/**
 * \brief Free the data pointed to by condemned of the given size.
 *
 * \param condemned Pointer to memory.
 * \param size Number of bytes.
 */
void m4ri_mmc_free(void *condemned, size_t size) {
#if __M4RI_ENABLE_MMC

#if __M4RI_HAVE_OPENMP
#pragma omp critical(mmc)
  {
#endif
    static int j = 0;
    mmb_t *mm    = m4ri_mmc_cache;
    if (size < __M4RI_MMC_THRESHOLD) {
      for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
        if (mm[i].size == 0) {
          mm[i].size = size;
          mm[i].data = condemned;
          goto done;
        }
      }
      m4ri_mm_free(mm[j].data);
      mm[j].size = size;
      mm[j].data = condemned;
      j          = (j + 1) % __M4RI_MMC_NBLOCKS;
    } else {
      m4ri_mm_free(condemned);
    }
  done:;
#if __M4RI_HAVE_OPENMP
  }
#endif  // __M4RI_HAVE_OPENMP
#else   // __M4RI_ENABLE_MMC
  m4ri_mm_free(condemned);
#endif  // __M4RI_ENABLE_MMC
}

/**
 * \brief Cleans up memory block cache.

 *
 * This function is called automatically when the shared library is unloaded.
 *
 * \warning Not thread safe.
 */
void m4ri_mmc_cleanup(void) {
#if __M4RI_ENABLE_MMC

#if __M4RI_HAVE_OPENMP
#pragma omp critical(mmc)
  {
#endif
    mmb_t *mm = m4ri_mmc_cache;
    for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
      if (mm[i].size) m4ri_mm_free(mm[i].data);
      mm[i].size = 0;
    }
#if __M4RI_HAVE_OPENMP
  }
#endif  // __M4RI_HAVE_OPENMP
#endif  // __M4RI_ENABLE_MMC
}
