#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "testing.h"
#include <m4ri/m4ri_config.h>
#include <m4ri/m4ri.h>
#include <setjmp.h>
#include <stdlib.h>
#include <unistd.h>

#if __M4RI_HAVE_OPENMP
#include <omp.h>
#endif

static unsigned long ncalls;
static unsigned long jump_after;
static int bad_context;
static jmp_buf env;
static unsigned long npolls;
static unsigned long pending_after;

static void count_calls(void) {
#if __M4RI_HAVE_OPENMP
  if (omp_get_level() != 0) bad_context = 1;
#endif
  ++ncalls;
}

static void jump_out(void) {
  count_calls();
  if (ncalls >= jump_after) longjmp(env, 1);
}

/* pending after a number of polls; stays pending */
static int pending(void) {
  return __atomic_add_fetch(&npolls, 1, __ATOMIC_RELAXED) > pending_after;
}

static void jump_if_pending(void) {
  count_calls();
  if (__atomic_load_n(&npolls, __ATOMIC_RELAXED) > pending_after) longjmp(env, 1);
}

/* reference results are computed without any interrupt check installed */
static mzd_t *solve(mzd_t const *A, mzd_t const *B) {
  mzd_t *Acopy = mzd_copy(NULL, A);
  mzd_t *X     = mzd_copy(NULL, B);
  mzd_solve_left(Acopy, X, 0, 0);
  mzd_free(Acopy);
  return X;
}

static int check_called(char const *what, unsigned long before) {
  printf("interrupt check %-12s calls: %6lu ", what, ncalls - before);
  if (ncalls == before || bad_context) {
    printf("FAILED\n");
    return 1;
  }
  printf("passed\n");
  return 0;
}

int test_interrupt_check_called(rci_t n) {
  int status = 0;
  mzd_t *A   = mzd_init(n, n);
  mzd_t *B   = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);
  for (rci_t i = 0; i < n; ++i) mzd_write_bit(A, i, i, 1);

  mzd_t *C0 = mzd_mul_m4rm(NULL, A, B, 0);
  mzd_t *X0 = solve(A, B);

  ncalls      = 0;
  bad_context = 0;
  m4ri_set_interrupt_check(count_calls, NULL);

  unsigned long before = ncalls;
  mzd_t *E             = mzd_copy(NULL, A);
  mzd_echelonize_m4ri(E, 1, 0);
  status += check_called("echelonize", before);

  before    = ncalls;
  mzd_t *C1 = mzd_mul_m4rm(NULL, A, B, 0);
  status += check_called("mul_m4rm", before);

  before    = ncalls;
  mzd_t *C2 = mzd_mul(NULL, A, B, 512);
  status += check_called("mul", before);

  before    = ncalls;
  mzp_t *P  = mzp_init(n);
  mzp_t *Q  = mzp_init(n);
  mzd_t *LU = mzd_copy(NULL, A);
  mzd_pluq(LU, P, Q, 0);
  status += check_called("pluq", before);

  before    = ncalls;
  mzd_t *X1 = solve(A, B);
  status += check_called("solve_left", before);

  m4ri_set_interrupt_check(NULL, NULL);

  /* a check which returns must not change the results */
  printf("interrupt check results unchanged ");
  if (!mzd_equal(C0, C1) || !mzd_equal(C0, C2) || !mzd_equal(X0, X1)) {
    printf("FAILED\n");
    status += 1;
  } else {
    printf("passed\n");
  }

  mzd_free(E);
  mzd_free(C0);
  mzd_free(C1);
  mzd_free(C2);
  mzd_free(LU);
  mzd_free(X0);
  mzd_free(X1);
  mzp_free(P);
  mzp_free(Q);
  mzd_free(A);
  mzd_free(B);
  return status;
}

int test_interrupt_longjmp(rci_t n, unsigned long after) {
  mzd_t *A = mzd_init(n, n);
  mzd_t *B = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);
  for (rci_t i = 0; i < n; ++i) mzd_write_bit(A, i, i, 1);
  mzd_t *X0 = solve(A, B);

  ncalls      = 0;
  jump_after  = after;
  bad_context = 0;
  int jumped  = 0;
  /* leaked when we jump out, as documented */
  mzd_t *volatile Acopy = mzd_copy(NULL, A);
  mzd_t *volatile X     = mzd_copy(NULL, B);

  m4ri_set_interrupt_check(jump_out, NULL);
  if (setjmp(env) == 0) {
    mzd_solve_left(Acopy, X, 0, 0);
  } else {
    jumped = 1;
  }
  m4ri_set_interrupt_check(NULL, NULL);

  /* computations after the jump, including OpenMP parallel ones, must
     work as before */
  mzd_t *X1 = solve(A, B);

  printf("interrupt longjmp n: %5d after: %3lu ", n, after);
  int status = 0;
  if (!jumped || bad_context || !mzd_equal(X0, X1)) {
    printf("FAILED\n");
    status = 1;
  } else {
    printf("passed\n");
  }
  mzd_free(X0);
  mzd_free(X1);
  mzd_free(A);
  mzd_free(B);
  return status;
}

int test_interrupt_set(void) {
  int status = 0;
  printf("interrupt set/remove ");

  /* no check installed: nothing is called */
  m4ri_set_interrupt_check(count_calls, pending);
  m4ri_set_interrupt_check(NULL, pending);
  ncalls   = 0;
  npolls   = 0;
  mzd_t *A = mzd_init(512, 512);
  mzd_randomize(A);
  mzd_echelonize_m4ri(A, 1, 0);
  mzd_free(A);
  if (ncalls != 0 || npolls != 0 || m4ri_interrupt_pending()) status = 1;

  printf(status ? "FAILED\n" : "passed\n");
  return status;
}

#if __M4RI_HAVE_OPENMP
/* an interrupt becoming pending inside a parallel loop cancels it */
int test_interrupt_cancel(rci_t n, unsigned long after) {
  mzd_t *A = mzd_init(n, n);
  mzd_t *B = mzd_init(n, n);
  mzd_randomize(A);
  mzd_randomize(B);
  mzd_t *C0 = mzd_mul_m4rm(NULL, A, B, 0);

  ncalls        = 0;
  npolls        = 0;
  pending_after = after;
  bad_context   = 0;
  int jumped    = 0;

  m4ri_set_interrupt_check(jump_if_pending, pending);
  if (setjmp(env) == 0) {
    mzd_free(mzd_mul_m4rm(NULL, A, B, 0));
  } else {
    jumped = 1;
  }
  unsigned long const polls = npolls;
  m4ri_set_interrupt_check(NULL, NULL);

  /* computations after the jump must work as before */
  mzd_t *C1 = mzd_mul_m4rm(NULL, A, B, 0);

  /* every row of every column block is polled once without
     cancellation; with it, each thread stops after at most one more
     poll per row once the interrupt is pending */
  int cancel_ok = 1;
  if (omp_get_cancellation())
    cancel_ok = polls <= after + 2 * (unsigned long)omp_get_max_threads();

  printf("interrupt cancel n: %5d after: %6lu polls: %8lu ", n, after, polls);
  int status = 0;
  if (!jumped || bad_context || !cancel_ok || !mzd_equal(C0, C1)) {
    printf("FAILED\n");
    status = 1;
  } else {
    printf("passed\n");
  }
  mzd_free(C0);
  mzd_free(C1);
  mzd_free(A);
  mzd_free(B);
  return status;
}
#endif

int main(int argc, char **argv) {
  int status = 0;

#if __M4RI_HAVE_OPENMP
  /* OMP_CANCELLATION is only read at startup */
  if (!omp_get_cancellation() && getenv("M4RI_TEST_INTERRUPT_REEXEC") == NULL) {
    setenv("OMP_CANCELLATION", "true", 1);
    setenv("M4RI_TEST_INTERRUPT_REEXEC", "1", 1);
    execv(argv[0], argv);
  }
  printf("OpenMP cancellation: %s\n", omp_get_cancellation() ? "enabled" : "disabled");
#endif

  srandom(17);

  status += test_interrupt_set();
  status += test_interrupt_check_called(1024);
  status += test_interrupt_check_called(2500);
  status += test_interrupt_longjmp(2500, 1);
  status += test_interrupt_longjmp(2500, 7);
  status += test_interrupt_longjmp(2500, 40);
#if __M4RI_HAVE_OPENMP
  status += test_interrupt_cancel(4000, 1);
  status += test_interrupt_cancel(4000, 3000);
#endif

  if (!status) {
    printf("All tests passed.\n");
  } else {
    return 1;
  }

  return 0;
}
