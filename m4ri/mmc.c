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

#if __M4RI_ENABLE_MMC && !__M4RI_HAVE_OPENMP
#define M4RI_MMC_CACHE_ACTIVE 1
#else
#define M4RI_MMC_CACHE_ACTIVE 0
#endif

/**
 * \brief Allocate size bytes.
 *
 * \param size Number of bytes.
 *
 * \return pointer to allocated memory block.
 */

void *m4ri_mmc_malloc(size_t size) {

#if M4RI_MMC_CACHE_ACTIVE
  if (size == 0) return m4ri_mm_malloc(size);

  mmb_t *mm           = m4ri_mmc_cache;
  volatile mmb_t *vmm = (volatile mmb_t *)m4ri_mmc_cache;
  if (size <= __M4RI_MMC_THRESHOLD) {
    for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
      if (mm[i].size == size && mm[i].data != NULL) {
        void *ret   = mm[i].data;
        vmm[i].size = 0;
        vmm[i].data = NULL;
        return ret;
      }
    }
  }
  return m4ri_mm_malloc(size);

#else  // M4RI_MMC_CACHE_ACTIVE

  return m4ri_mm_malloc(size);

#endif  // M4RI_MMC_CACHE_ACTIVE
}

/**
 * \brief Free the data pointed to by condemned of the given size.
 *
 * \param condemned Pointer to memory.
 * \param size Number of bytes.
 */
void m4ri_mmc_free(void *condemned, size_t size) {
#if M4RI_MMC_CACHE_ACTIVE
  if (condemned == NULL || size == 0 || size >= __M4RI_MMC_THRESHOLD) {
    m4ri_mm_free(condemned);
    return;
  }

  mmb_t *mm           = m4ri_mmc_cache;
  /* Keep each slot valid at every interruptible store: size is published last. */
  volatile mmb_t *vmm = (volatile mmb_t *)m4ri_mmc_cache;
  int empty            = -1;
  for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (mm[i].data == condemned) {
      /* A caller may retry after interruption during or after publication. */
      if (mm[i].size == 0) vmm[i].size = size;
      return;
    }
    if (empty < 0 && mm[i].data == NULL) empty = i;
  }
  if (empty >= 0) {
    vmm[empty].size = 0;
    vmm[empty].data = condemned;
    vmm[empty].size = size;
    return;
  }
  m4ri_mm_free(condemned);
#else   // M4RI_MMC_CACHE_ACTIVE
  m4ri_mm_free(condemned);
#endif  // M4RI_MMC_CACHE_ACTIVE
}

/**
 * \brief Cleans up memory block cache.

 *
 * This function is called automatically when the shared library is unloaded.
 *
 * \warning Not thread safe.
 */
void m4ri_mmc_cleanup(void) {
#if M4RI_MMC_CACHE_ACTIVE

  mmb_t *mm           = m4ri_mmc_cache;
  volatile mmb_t *vmm = (volatile mmb_t *)m4ri_mmc_cache;
  for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (mm[i].data) {
      void *condemned = mm[i].data;
      vmm[i].size    = 0;
      vmm[i].data    = NULL;
      m4ri_mm_free(condemned);
    } else {
      vmm[i].size = 0;
    }
  }
#endif  // M4RI_MMC_CACHE_ACTIVE
}
