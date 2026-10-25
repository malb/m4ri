/**
 * \file interrupt.h
 *
 * \brief Interrupting long computations.
 *
 * Callers may install functions with m4ri_set_interrupt_check() which
 * M4RI calls at points where it is safe to leave a computation, e.g.
 * by longjmp(), and which let it cancel OpenMP parallel loops early.
 */

#ifndef M4RI_INTERRUPT_H
#define M4RI_INTERRUPT_H

#include <m4ri/m4ri_config.h>
#include <stddef.h>

/**
 * \brief Defined if m4ri_set_interrupt_check() is available.
 */

#define __M4RI_HAVE_INTERRUPT_CHECK 1

/**
 * \brief Type of a user supplied interrupt check.
 */

typedef void (*m4ri_interrupt_check_t)(void);

/**
 * \brief Type of a user supplied query whether an interrupt is pending.
 */

typedef int (*m4ri_interrupt_pending_t)(void);

/**
 * \brief Install functions to interrupt long computations.
 *
 * check is called periodically. It is only ever called from the
 * thread which called into M4RI, outside of any OpenMP parallel region
 * and while M4RI holds no locks. It is therefore safe for it to not
 * return, e.g. to longjmp() out of the computation in order to
 * interrupt it. This is not true for asynchronous signal handlers,
 * which may fire while OpenMP worker threads are still running or
 * while a lock in malloc() is held.
 *
 * pending is optional. It is polled by all threads inside OpenMP
 * parallel loops, so it must be cheap and thread-safe. Once it
 * returns nonzero, it must keep doing so until check is called. These
 * loops are then cancelled (if cancellation is enabled, see
 * omp_get_cancellation() and OMP_CANCELLATION) and their results are
 * incomplete, so the next call of check must not return. M4RI calls
 * m4ri_die() if it does.
 *
 * If the computation is left, memory allocated by it is leaked and the
 * matrices it operated on are left in an unspecified state.
 *
 * \param check The function to call or NULL to remove both functions.
 * \param pending The query or NULL.
 */

void m4ri_set_interrupt_check(m4ri_interrupt_check_t check, m4ri_interrupt_pending_t pending);

/**
 * \brief Call the check function installed by
 * m4ri_set_interrupt_check(), if any, unless we are inside an OpenMP
 * parallel region.
 *
 * Only call this where it is safe to never return, see
 * m4ri_set_interrupt_check().
 */

void m4ri_check_interrupt(void);

/**
 * \brief The query installed by m4ri_set_interrupt_check(); do not use
 * directly.
 */

extern m4ri_interrupt_pending_t m4ri_interrupt_pending_query;

/**
 * \brief Return whether an interrupt is pending, see
 * m4ri_set_interrupt_check(). May be called from any thread.
 */

static inline int m4ri_interrupt_pending(void) {
  m4ri_interrupt_pending_t const f = m4ri_interrupt_pending_query;
  return f != NULL && f();
}

/**
 * \brief Cancel the enclosing OpenMP parallel for loop if an interrupt
 * is pending, see m4ri_set_interrupt_check().
 *
 * Use at the start of the loop body and call m4ri_check_interrupt()
 * after the loop.
 */

#if __M4RI_HAVE_OPENMP
#define __M4RI_CANCEL_FOR_IF_INTERRUPTED()                                                         \
  do {                                                                                             \
    if (m4ri_interrupt_pending()) { _Pragma("omp cancel for") }                                    \
    _Pragma("omp cancellation point for")                                                          \
  } while (0)
#else
#define __M4RI_CANCEL_FOR_IF_INTERRUPTED() \
  do {                                     \
  } while (0)
#endif

#endif  // M4RI_INTERRUPT_H
