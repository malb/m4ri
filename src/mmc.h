/**
 * \file mmc.h
 * \brief The mmc memory management functions check a cache for re-usable unused memory before asking the system for it.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */

#ifndef M4RI_MMC_H
#define M4RI_MMC_H

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

#include <m4ri/misc.h>

void *m4ri_mmc_malloc(size_t size);
void m4ri_mmc_free(void *condemned, size_t size);
void m4ri_mmc_cleanup(void);

/**
 * \brief Enable memory block cache (default: enabled).
 */
#define __M4RI_ENABLE_MMC

/**
 * \brief Number of blocks that are cached.
 */
#define __M4RI_MMC_NBLOCKS 16

/**
 * \brief Maximal size of blocks stored in cache.
 */
#define __M4RI_MMC_THRESHOLD __M4RI_CPU_L3_CACHE

/**
 * \brief Tuple of pointer to allocated memory block and it's size.
 */
typedef struct _mm_block {
  /**
   * Size in bytes of the data.
   */
  size_t size;

  /**
   * Pointer to buffer of data.
   */
  void *data;

} mmb_t;

/**
 * \brief Allocate an array of count times size zeroed bytes.
 *
 * \param count Number of elements.
 * \param size Number of bytes per element.
 *
 * \return Pointer to allocated memory block.
 */
static inline void *m4ri_mmc_calloc(size_t count, size_t size) {
  size_t total_size = count * size;
  void *ret = m4ri_mmc_malloc(total_size);
  memset((char*)ret, 0, total_size);
  return ret;
}

#endif // M4RI_MMC_H
