/******************************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
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

#include "interrupt.h"
#include "misc.h"

#if __M4RI_HAVE_OPENMP
#include <omp.h>
#endif

static m4ri_interrupt_check_t m4ri_interrupt_check = NULL;
m4ri_interrupt_pending_t m4ri_interrupt_pending_query = NULL;

void m4ri_set_interrupt_check(m4ri_interrupt_check_t check, m4ri_interrupt_pending_t pending) {
  m4ri_interrupt_check         = check;
  m4ri_interrupt_pending_query = check ? pending : NULL;
}

void m4ri_check_interrupt(void) {
  m4ri_interrupt_check_t f = m4ri_interrupt_check;
  if (f == NULL) return;
#if __M4RI_HAVE_OPENMP
  /* omp_get_level() also counts inactive (e.g. single-threaded)
     parallel regions, which must not be left via longjmp() either. */
  if (omp_get_level() != 0) return;
#endif
  int const pending = m4ri_interrupt_pending();
  f();
  /* parallel loops may have been cancelled, so we cannot continue */
  if (pending) m4ri_die("m4ri_check_interrupt: check returned while an interrupt is pending\n");
}
