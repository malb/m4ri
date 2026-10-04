/*
 * test_mmc.c
 *
 * Testing the memory management cache.
 *
 * Copyright (C) 2026
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */

#include "testing.h"
#include <m4ri/mmc.h>
#include <stdio.h>
#include <string.h>

#if __M4RI_ENABLE_MMC && !__M4RI_HAVE_OPENMP
extern mmb_t m4ri_mmc_cache[__M4RI_MMC_NBLOCKS];

static int cache_has_inconsistent_empty_slot(void) {
  for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (m4ri_mmc_cache[i].data == NULL && m4ri_mmc_cache[i].size != 0) return 1;
  }
  return 0;
}

static int cache_is_empty(void) {
  for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (m4ri_mmc_cache[i].data != NULL || m4ri_mmc_cache[i].size != 0) return 0;
  }
  return 1;
}

static int cache_pointer_count(void *p) {
  int count = 0;
  for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (m4ri_mmc_cache[i].data == p) ++count;
  }
  return count;
}

static void cache_remove_duplicate_pointers(void *p) {
  int found = 0;
  for (int i = 0; i < __M4RI_MMC_NBLOCKS; ++i) {
    if (m4ri_mmc_cache[i].data != p) continue;
    if (!found) {
      found = 1;
    } else {
      m4ri_mmc_cache[i].data = NULL;
      m4ri_mmc_cache[i].size = 0;
    }
  }
}
#endif

static int test_mmc_direct_alloc(void) {
  const size_t size = 257;
  void *p           = m4ri_mmc_malloc(size);
  memset(p, 0xa5, size);
  m4ri_mmc_free(p, size);
  return 0;
}

#if __M4RI_ENABLE_MMC && !__M4RI_HAVE_OPENMP
static int test_mmc_cache_reuse(void) {
  const size_t size = 257;
  m4ri_mmc_cleanup();

  void *p = m4ri_mmc_malloc(size);
  memset(p, 0xa5, size);
  m4ri_mmc_free(p, size);
  if (cache_has_inconsistent_empty_slot()) return 1;

  void *q = m4ri_mmc_malloc(size);
  if (q != p) return 1;
  if (cache_has_inconsistent_empty_slot()) return 1;

  m4ri_mmc_free(q, size);
  m4ri_mmc_cleanup();
  return !cache_is_empty();
}

static int test_mmc_repairs_stale_empty_slot(void) {
  const size_t size = 257;
  m4ri_mmc_cleanup();

  m4ri_mmc_cache[0].data = NULL;
  m4ri_mmc_cache[0].size = 2 * size;

  void *p = m4ri_mmc_malloc(size);
  m4ri_mmc_free(p, size);
  if (m4ri_mmc_cache[0].data != p) return 1;
  if (m4ri_mmc_cache[0].size != size) return 1;

  p = m4ri_mmc_malloc(size);
  m4ri_mmc_free(p, size);
  m4ri_mmc_cleanup();
  return cache_has_inconsistent_empty_slot();
}

static int test_mmc_cleanup_reclaims_unpublished_slot(void) {
  const size_t size = 257;
  m4ri_mmc_cleanup();

  m4ri_mmc_cache[0].data = m4ri_mm_malloc(size);
  m4ri_mmc_cache[0].size = 0;
  m4ri_mmc_cleanup();

  return !cache_is_empty();
}

static int test_mmc_free_retry_is_idempotent(void) {
  const size_t size = 257;
  int status        = 0;
  m4ri_mmc_cleanup();

  void *p                     = m4ri_mm_malloc(size);
  m4ri_mmc_cache[0].size      = 0;
  m4ri_mmc_cache[0].data      = p;
  m4ri_mmc_free(p, size);
  status += cache_pointer_count(p) != 1;
  status += m4ri_mmc_cache[0].size != size;

  cache_remove_duplicate_pointers(p);
  m4ri_mmc_cache[1]      = m4ri_mmc_cache[0];
  m4ri_mmc_cache[0].data = NULL;
  m4ri_mmc_cache[0].size = 0;
  m4ri_mmc_free(p, size);
  status += cache_pointer_count(p) != 1;

  cache_remove_duplicate_pointers(p);
  m4ri_mmc_cleanup();
  return status || !cache_is_empty();
}
#endif

#if __M4RI_HAVE_OPENMP
static int test_mmc_openmp_disables_cache(void) {
  return __M4RI_ENABLE_MMC;
}
#endif

int main(int argc, char *argv[]) {
  int status = 0;

  status += test_mmc_direct_alloc();

#if __M4RI_ENABLE_MMC && !__M4RI_HAVE_OPENMP
  status += test_mmc_cache_reuse();
  status += test_mmc_repairs_stale_empty_slot();
  status += test_mmc_cleanup_reclaims_unpublished_slot();
  status += test_mmc_free_retry_is_idempotent();
#endif

#if __M4RI_HAVE_OPENMP
  status += test_mmc_openmp_disables_cache();
#endif

  if (!status) {
    printf("All tests passed.\n");
  } else {
    printf("TEST FAILED!\n");
    return 1;
  }
  return 0;
}
