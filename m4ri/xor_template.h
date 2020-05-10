#include <m4ri/m4ri_config.h>
#include <m4ri/misc.h>

/**
 * Compute c[i] += sum(t[j][i], 0 <= j < N) for 0 <= i < wide
 *
 * \todo the non SSE2 version of this code is slow, replace by code from mzd_process_rows8
 *
 * \warn Assumes __M4RI_ALIGNMENT(c, 16) == __M4RI_ALIGNMENT(t[i], 16)
 */

static inline void __M4RI_TEMPLATE_NAME(_mzd_combine)(word *m, word const *t[N], wi_t wide) {
  assert(1 <= N && N <= 8);

#if __M4RI_HAVE_SSE2

  assert((__M4RI_ALIGNMENT(m, 16) == 8) | (__M4RI_ALIGNMENT(m, 16) == 0));

  switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
  case 8: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[7], 16));
  case 7: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[6], 16));
  case 6: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[5], 16));
  case 5: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[4], 16));
  case 4: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[3], 16));
  case 3: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[2], 16));
  case 2: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[1], 16));
  case 1: assert(__M4RI_ALIGNMENT(m, 16) == __M4RI_ALIGNMENT(t[0], 16));
  };

  if (__M4RI_UNLIKELY(__M4RI_ALIGNMENT(m, 16) == 8)) {
    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++ ^ *t[6]++ ^ *t[7]++;
      break;
    case 7: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++ ^ *t[6]++; break;
    case 6: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++; break;
    case 5: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++; break;
    case 4: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++; break;
    case 3: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++; break;
    case 2: *m++ ^= *t[0]++ ^ *t[1]++; break;
    case 1: *m++ ^= *t[0]++; break;
    };
    wide--;
  }

  __m128i *m__ = (__m128i *)m;
  __m128i *t__[N];

  switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
  case 8: t__[N - 8] = (__m128i *)t[N - 8];
  case 7: t__[N - 7] = (__m128i *)t[N - 7];
  case 6: t__[N - 6] = (__m128i *)t[N - 6];
  case 5: t__[N - 5] = (__m128i *)t[N - 5];
  case 4: t__[N - 4] = (__m128i *)t[N - 4];
  case 3: t__[N - 3] = (__m128i *)t[N - 3];
  case 2: t__[N - 2] = (__m128i *)t[N - 2];
  case 1: t__[N - 1] = (__m128i *)t[N - 1];
  };

  __m128i xmm0, xmm1, xmm2, xmm3;

  wi_t i = 0;
  for (; i + 4 <= (wide >> 1); i += 4) {
    xmm0 = m__[0];
    xmm1 = m__[1];
    xmm2 = m__[2];
    xmm3 = m__[3];
    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      xmm0 = _mm_xor_si128(xmm0, t__[7][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[7][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[7][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[7][3]);
      t__[7] += 4;
    case 7:
      xmm0 = _mm_xor_si128(xmm0, t__[6][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[6][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[6][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[6][3]);
      t__[6] += 4;
    case 6:
      xmm0 = _mm_xor_si128(xmm0, t__[5][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[5][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[5][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[5][3]);
      t__[5] += 4;
    case 5:
      xmm0 = _mm_xor_si128(xmm0, t__[4][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[4][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[4][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[4][3]);
      t__[4] += 4;
    case 4:
      xmm0 = _mm_xor_si128(xmm0, t__[3][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[3][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[3][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[3][3]);
      t__[3] += 4;
    case 3:
      xmm0 = _mm_xor_si128(xmm0, t__[2][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[2][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[2][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[2][3]);
      t__[2] += 4;
    case 2:
      xmm0 = _mm_xor_si128(xmm0, t__[1][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[1][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[1][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[1][3]);
      t__[1] += 4;
    case 1:
      xmm0 = _mm_xor_si128(xmm0, t__[0][0]);
      xmm1 = _mm_xor_si128(xmm1, t__[0][1]);
      xmm2 = _mm_xor_si128(xmm2, t__[0][2]);
      xmm3 = _mm_xor_si128(xmm3, t__[0][3]);
      t__[0] += 4;
    }
    m__[0] = xmm0;
    m__[1] = xmm1;
    m__[2] = xmm2;
    m__[3] = xmm3;
    m__ += 4;
  }

  for (; i < (wide >> 1); i++) {
    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm1 = _mm_xor_si128(*t__[2]++, *t__[3]++);
      xmm2 = _mm_xor_si128(*t__[4]++, *t__[5]++);
      xmm3 = _mm_xor_si128(*t__[6]++, *t__[7]++);
      xmm0 = _mm_xor_si128(xmm0, xmm1);
      xmm2 = _mm_xor_si128(xmm2, xmm3);
      xmm0 = _mm_xor_si128(xmm0, xmm2);
      xmm0 = _mm_xor_si128(*m__, xmm0);
      break;
    case 7:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm1 = _mm_xor_si128(*t__[2]++, *t__[3]++);
      xmm0 = _mm_xor_si128(xmm0, *t__[4]++);
      xmm1 = _mm_xor_si128(xmm1, *t__[5]++);
      xmm0 = _mm_xor_si128(xmm0, *t__[6]++);
      xmm0 = _mm_xor_si128(xmm0, xmm1);
      xmm0 = _mm_xor_si128(*m__, xmm0);
      break;
    case 6:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm1 = _mm_xor_si128(*t__[2]++, *t__[3]++);
      xmm0 = _mm_xor_si128(xmm0, *t__[4]++);
      xmm1 = _mm_xor_si128(xmm1, *t__[5]++);
      xmm0 = _mm_xor_si128(xmm0, xmm1);
      xmm0 = _mm_xor_si128(*m__, xmm0);
      break;
    case 5:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm1 = _mm_xor_si128(*t__[2]++, *t__[3]++);
      xmm0 = _mm_xor_si128(xmm0, *t__[4]++);
      xmm0 = _mm_xor_si128(xmm0, xmm1);
      xmm0 = _mm_xor_si128(*m__, xmm0);
      break;
    case 4:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm1 = _mm_xor_si128(*t__[2]++, *t__[3]++);
      xmm0 = _mm_xor_si128(xmm0, xmm1);
      xmm0 = _mm_xor_si128(*m__, xmm0);
      break;
    case 3:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm1 = _mm_xor_si128(*m__, *t__[2]++);
      xmm0 = _mm_xor_si128(xmm0, xmm1);
      break;
    case 2:
      xmm0 = _mm_xor_si128(*t__[0]++, *t__[1]++);
      xmm0 = _mm_xor_si128(*m__, xmm0);
      break;
    case 1: xmm0 = _mm_xor_si128(*m__, *t__[0]++); break;
    };
    *m__++ = xmm0;
  }

  if (wide & 0x1) {
    m = (word *)m__;
    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8: t[N - 8] = (word *)t__[N - 8];
    case 7: t[N - 7] = (word *)t__[N - 7];
    case 6: t[N - 6] = (word *)t__[N - 6];
    case 5: t[N - 5] = (word *)t__[N - 5];
    case 4: t[N - 4] = (word *)t__[N - 4];
    case 3: t[N - 3] = (word *)t__[N - 3];
    case 2: t[N - 2] = (word *)t__[N - 2];
    case 1: t[N - 1] = (word *)t__[N - 1];
    }

    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++ ^ *t[6]++ ^ *t[7]++;
      break;
    case 7: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++ ^ *t[6]++; break;
    case 6: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++; break;
    case 5: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++; break;
    case 4: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++; break;
    case 3: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++; break;
    case 2: *m++ ^= *t[0]++ ^ *t[1]++; break;
    case 1: *m++ ^= *t[0]++; break;
    }
  }
  return;
#else

  for (wi_t i = 0; i < wide; i++) {
    switch (N) { /* we rely on the compiler to optimise this switch away, it reads nicer than #if */
    case 8:
      *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++ ^ *t[6]++ ^ *t[7]++;
      break;
    case 7: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++ ^ *t[6]++; break;
    case 6: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++ ^ *t[5]++; break;
    case 5: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++ ^ *t[4]++; break;
    case 4: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++ ^ *t[3]++; break;
    case 3: *m++ ^= *t[0]++ ^ *t[1]++ ^ *t[2]++; break;
    case 2: *m++ ^= *t[0]++ ^ *t[1]++; break;
    case 1: *m++ ^= *t[0]++; break;
    }
  }

  return;
#endif  // __M4RI_HAVE_SSE2
}
