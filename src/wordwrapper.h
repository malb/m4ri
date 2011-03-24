/**
 * \file wordwrapper.h
 *
 * \brief C++ class wrapper for a word.
 * 
 * \author Carlo Wood <carlo@alinoe.com>
 *
 * To use the wrapper class, configure with CC (not CXX) set to a C++ compiler.
 * For example:
 * 
 * CFLAGS="-O2" CC="g++" ./configure --enable-debug
 */
/******************************************************************************
*
*                 M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2010 Carlo Wood <carlo@alinoe.com>
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

#include <limits.h>

struct wi_t {
    int M_val;

    wi_t(void) { }
    wi_t(unsigned int val) : M_val(val) { }
    wi_t(wi_t const& wi) : M_val(wi.M_val) { }
    wi_t& operator=(unsigned int val) { M_val = val; return *this; }
    wi_t& operator=(wi_t const& wi) { M_val = wi.M_val; return *this; }

    wi_t& operator+=(wi_t const& wi) { M_val += wi.M_val; return *this; }
    wi_t& operator-=(wi_t const& wi) { M_val -= wi.M_val; return *this; }

    int operator%(unsigned int m) const { return M_val % m; }
    wi_t operator/(unsigned int d) const { wi_t result; result.M_val = M_val / d; return result; }

    friend bool operator==(wi_t const& wi1, wi_t const& wi2) { return wi1.M_val == wi2.M_val; }
    friend bool operator==(unsigned int wi1, wi_t const& wi2) { return (int)wi1 == wi2.M_val; }
    friend bool operator==(wi_t const& wi1, unsigned int wi2) { return wi1.M_val == (int)wi2; }
    friend bool operator!=(wi_t const& wi1, wi_t const& wi2) { return wi1.M_val != wi2.M_val; }
    friend bool operator!=(unsigned int wi1, wi_t const& wi2) { return (int)wi1 != wi2.M_val; }
    friend bool operator!=(wi_t const& wi1, unsigned int wi2) { return wi1.M_val != (int)wi2; }
    friend bool operator<(wi_t const& wi1, wi_t const& wi2) { return wi1.M_val < wi2.M_val; }
    friend bool operator<(unsigned int wi1, wi_t const& wi2) { return (int)wi1 < wi2.M_val; }
    friend bool operator<(wi_t const& wi1, unsigned int wi2) { return wi1.M_val < (int)wi2; }
    friend bool operator>(wi_t const& wi1, wi_t const& wi2) { return wi1.M_val > wi2.M_val; }
    friend bool operator>(unsigned int wi1, wi_t const& wi2) { return (int)wi1 > wi2.M_val; }
    friend bool operator>(wi_t const& wi1, unsigned int wi2) { return wi1.M_val > (int)wi2; }
    friend bool operator<=(wi_t const& wi1, wi_t const& wi2) { return wi1.M_val <= wi2.M_val; }
    friend bool operator<=(unsigned int wi1, wi_t const& wi2) { return (int)wi1 <= wi2.M_val; }
    friend bool operator<=(wi_t const& wi1, unsigned int wi2) { return wi1.M_val <= (int)wi2; }
    friend bool operator>=(wi_t const& wi1, wi_t const& wi2) { return wi1.M_val >= wi2.M_val; }
    friend bool operator>=(unsigned int wi1, wi_t const& wi2) { return (int)wi1 >= wi2.M_val; }
    friend bool operator>=(wi_t const& wi1, unsigned int wi2) { return wi1.M_val >= (int)wi2; }

    friend wi_t operator>>(wi_t const& wi, int shift) { assert(shift >= 0 && shift < 32); wi_t result; result.M_val = wi.M_val >> shift; return result; }

    wi_t& operator++() { ++M_val; return *this; }
    wi_t operator++(int) { wi_t cur(*this); ++M_val; return cur; }
    wi_t& operator--() { --M_val; return *this; }
    wi_t operator--(int) { wi_t cur(*this); --M_val; return cur; }

    friend wi_t operator-(wi_t const& wi1, wi_t const& wi2) { wi_t result; result.M_val = wi1.M_val - wi2.M_val; return result; }
    friend wi_t operator+(wi_t const& wi1, wi_t const& wi2) { wi_t result; result.M_val = wi1.M_val + wi2.M_val; return result; }

    int val(void) const { return M_val; }

  private:
    wi_t(int val) : M_val(val) { }
    wi_t& operator=(int val) { M_val = val; return *this; }
};

class rci_t;

struct Radix_t {
    int M_val;

    Radix_t(unsigned int val) : M_val(val) { assert(M_val % 64 == 0); assert(val != 0); }
    Radix_t(Radix_t const& r1) : M_val(r1.M_val) { }
    Radix_t& operator=(Radix_t const& r1) { M_val = r1.M_val; return *this; }
    Radix_t& operator=(int val) { assert(val % 64 == 0); assert(val != 0); M_val = val; return *this; }

    friend bool operator==(Radix_t const& r1, Radix_t const& r2) { return r1.M_val == r2.M_val; }
    friend bool operator==(unsigned int r1, Radix_t const& r2) { return (int)r1 == r2.M_val; }
    friend bool operator==(Radix_t const& r1, unsigned int r2) { return r1.M_val == (int)r2; }
    friend bool operator!=(Radix_t const& r1, Radix_t const& r2) { return r1.M_val != r2.M_val; }
    friend bool operator<(Radix_t const& r1, Radix_t const& r2) { return r1.M_val < r2.M_val; }
    friend bool operator<(unsigned int r1, Radix_t const& r2) { return (int)r1 < r2.M_val; }
    friend bool operator<(Radix_t const& r1, unsigned int r2) { return r1.M_val < (int)r2; }
    friend bool operator>(Radix_t const& r1, Radix_t const& r2) { return r1.M_val > r2.M_val; }
    friend bool operator>(unsigned int r1, Radix_t const& r2) { return (int)r1 > r2.M_val; }
    friend bool operator>(Radix_t const& r1, unsigned int r2) { return r1.M_val > (int)r2; }
    friend bool operator<=(Radix_t const& r1, Radix_t const& r2) { return r1.M_val <= r2.M_val; }
    friend bool operator<=(unsigned int r1, Radix_t const& r2) { return (int)r1 <= r2.M_val; }
    friend bool operator<=(Radix_t const& r1, unsigned int r2) { return r1.M_val <= (int)r2; }
    friend bool operator>=(Radix_t const& r1, Radix_t const& r2) { return r1.M_val >= r2.M_val; }
    friend bool operator>=(unsigned int r1, Radix_t const& r2) { return (int)r1 >= r2.M_val; }
    friend bool operator>=(Radix_t const& r1, unsigned int r2) { return r1.M_val >= (int)r2; }

    Radix_t& operator+=(int offset) { M_val += offset; assert(offset % 64 == 0); return *this; }
    Radix_t& operator-=(int offset) { M_val -= offset; assert(offset % 64 == 0); return *this; }
    Radix_t& operator*=(int factor) { M_val *= factor; return *this; }
    friend int operator+(Radix_t const& r1, int offset) { return r1.M_val + offset; }
    friend int operator+(int offset, Radix_t const& r1) { return offset + r1.M_val; }
    friend int operator-(Radix_t const& r1, int offset) { return r1.M_val - offset; }
    friend int operator-(int offset, Radix_t const& r1) { return offset - r1.M_val; }
    friend Radix_t operator*(Radix_t const& r1, unsigned int factor) { Radix_t result(r1); result *= factor; return result; }
    friend Radix_t operator*(unsigned int factor, Radix_t const& r2) { Radix_t result(r2); result *= factor; return result; }
    friend int operator%(int offset, Radix_t const& r1) { return offset % r1.M_val; }

    int val(void) const { return M_val; }
};

extern int m4ri_opt_k(int, int,int);

class rci_t {
  private:
    int M_val;
  public:
    rci_t(void) { }
    rci_t(rci_t const& rci) : M_val(rci.M_val) { }
    ~rci_t() { }
    rci_t& operator=(rci_t const& rci) { M_val = rci.M_val; return *this; }
    rci_t& operator=(int zero) { assert(zero == 0); M_val = 0; return *this; }
    rci_t& operator=(unsigned int val) { M_val = val; return *this; }
    rci_t& operator=(Radix_t const& radix) { M_val = radix.M_val; return *this; }
    rci_t(int zero) : M_val(zero) { assert(zero == 0); }
    rci_t(uint64_t val) : M_val((int)val) { assert(val <= (uint64_t)INT_MAX); }
    rci_t(unsigned int val) : M_val(val) { }
    rci_t(Radix_t const& radix) : M_val(radix.M_val) { }
    friend bool operator==(rci_t const& rci1, rci_t const& rci2) { return rci1.M_val == rci2.M_val; }
    friend bool operator==(rci_t const& rci1, int zero) { assert(zero == 0); return rci1.M_val == 0; }
    friend bool operator==(rci_t const& rci1, unsigned int literal) { return rci1.M_val == literal; }
    friend bool operator==(Radix_t const& r1, rci_t const& rci) { return r1.M_val == rci.M_val; }
    friend bool operator==(rci_t const& rci, Radix_t const& r2) { return r2.M_val == rci.M_val; }
    friend bool operator!=(rci_t const& rci1, rci_t const& rci2) { return rci1.M_val != rci2.M_val; }
    friend bool operator!=(Radix_t const& r1, rci_t const& rci) { return r1.M_val != rci.M_val; }
    friend bool operator!=(rci_t const& rci, Radix_t const& r2) { return r2.M_val != rci.M_val; }
    friend bool operator<(rci_t const& rci1, rci_t const& rci2) { return rci1.M_val < rci2.M_val; }
    friend bool operator<(Radix_t const& r1, rci_t const& rci2) { return r1.M_val < rci2.M_val; }
    friend bool operator<(rci_t const& rci1, Radix_t const& r2) { return rci1.M_val < r2.M_val; }
    friend bool operator>(rci_t const& rci1, rci_t const& rci2) { return rci1.M_val > rci2.M_val; }
    friend bool operator>(Radix_t const& r1, rci_t const& rci2) { return r1.M_val > rci2.M_val; }
    friend bool operator>(rci_t const& rci1, Radix_t const& r2) { return rci1.M_val > r2.M_val; }
    friend bool operator<=(rci_t const& rci1, rci_t const& rci2) { return rci1.M_val <= rci2.M_val; }
    friend bool operator<=(Radix_t const& r1, rci_t const& rci2) { return r1.M_val <= rci2.M_val; }
    friend bool operator<=(rci_t const& rci1, Radix_t const& r2) { return rci1.M_val <= r2.M_val; }
    friend bool operator>=(rci_t const& rci1, rci_t const& rci2) { return rci1.M_val >= rci2.M_val; }
    friend bool operator>=(Radix_t const& r1, rci_t const& rci2) { return r1.M_val >= rci2.M_val; }
    friend bool operator>=(rci_t const& rci1, Radix_t const& r2) { return rci1.M_val >= r2.M_val; }
    friend bool operator>(rci_t const& rci1, int zero) { assert(zero == 0); return rci1.M_val > 0; }
    friend bool operator>=(rci_t const& rci1, int zero) { assert(zero == 0); return rci1.M_val >= 0; }
    friend bool operator<(rci_t const& rci1, int val) { return rci1.M_val < val; }
    friend bool operator<=(rci_t const& rci1, int val) { return rci1.M_val <= val; }
    int operator%(Radix_t radix) const { return M_val % radix.M_val; }
    rci_t operator%(rci_t const& rci) const { rci_t result; result.M_val = M_val % rci.M_val; return result; }
    int operator%(int m) const { return M_val % m; }
    int operator%(unsigned int m) const { return M_val % m; }
    friend rci_t operator|(rci_t const& rci1, rci_t const& rci2) { rci_t result; result.M_val = rci1.M_val | rci2.M_val; return result; }
    friend rci_t operator+(rci_t const& rci1, rci_t const& rci2) { rci_t result; result.M_val = rci1.M_val + rci2.M_val; return result; }
    friend rci_t operator+(Radix_t const& r1, rci_t const& rci2) { rci_t result; result.M_val = r1.M_val + rci2.M_val; return result; }
    friend rci_t operator+(rci_t const& rci1, Radix_t const& r2) { rci_t result; result.M_val = rci1.M_val + r2.M_val; return result; }
    rci_t& operator+=(rci_t const& rci) { M_val += rci.M_val; return *this; }
    rci_t& operator+=(Radix_t const& r1) { M_val += r1.M_val; return *this; }
    rci_t& operator+=(int val) { M_val += val; return *this; }
    rci_t& operator-=(rci_t const& rci) { M_val -= rci.M_val; return *this; }
    rci_t& operator-=(Radix_t const& r1) { M_val -= r1.M_val; return *this; }
    rci_t& operator-=(int val) { M_val -= val; return *this; }
    friend rci_t operator-(rci_t const& rci1, rci_t const& rci2) { rci_t result; result.M_val = rci1.M_val - rci2.M_val; return result; }
    friend rci_t operator-(Radix_t const& r1, rci_t const& rci2) { rci_t result; result.M_val = r1.M_val - rci2.M_val; return result; }
    friend rci_t operator-(rci_t const& rci1, Radix_t const& r2) { rci_t result; result.M_val = rci1.M_val - r2.M_val; return result; }
    friend rci_t operator-(rci_t const& rci, int offset) { rci_t result; result.M_val = rci.M_val - offset; return result; }
    friend rci_t operator+(rci_t const& rci, int offset) { rci_t result; result.M_val = rci.M_val + offset; return result; }
    rci_t& operator++() { ++M_val; return *this; }
    rci_t operator++(int) { rci_t cur(*this); ++M_val; return cur; }
    rci_t& operator--() { --M_val; return *this; }
    rci_t operator--(int) { rci_t cur(*this); --M_val; return cur; }
    int val(void) const { return M_val; }
    friend int m4ri_opt_k(rci_t a, rci_t b, int c) { return ::m4ri_opt_k(a.M_val, b.M_val, c); }
    friend int m4ri_opt_k(int a, rci_t b, rci_t c) { return ::m4ri_opt_k(a, b.M_val, c.M_val); }

    friend rci_t operator*(int factor, rci_t const& rci) { assert(factor >= 1 && factor <= 64); rci_t result(rci); result.M_val *= factor; return result; }
    friend rci_t operator*(unsigned int factor, rci_t const& rci) { rci_t result(rci); result.M_val *= factor; return result; }
    friend rci_t operator*(uint64_t factor, rci_t const& rci) { assert(factor >= 1 && factor <= 64); rci_t result(rci); result.M_val *= factor; return result; }
    friend size_t operator*(wi_t width, rci_t const& rci) { return (size_t)rci.M_val * width.val(); }
    friend double operator*(double factor, rci_t const& rci) { return factor * rci.M_val; }
    friend double operator/(rci_t const& rci, double denom) { return rci.M_val / denom; }
    friend rci_t operator/(rci_t const& rci, int div) { rci_t result; result.M_val = rci.M_val / div; return result; }
    friend rci_t operator/(rci_t const& rci, unsigned int div) { rci_t result; result.M_val = rci.M_val / div; return result; }
    rci_t& operator/=(int two) { assert(two == 2); M_val /= 2; return *this; }
    rci_t& operator*=(int two) { assert(two == 2); M_val *= 2; return *this; }

    // NOT operator. Returns true if the index (or difference) is zero.
    bool operator!(void) const { return !M_val; }
    // Automatic conversion to boolean.
    operator bool(void) const { return M_val != 0; }

  private:
    // Disallow conversion to int (this protects us from accidental conversion to bool (see above) and from there to int without us noticing that).
    operator int(void) const { assert(false); return 0; }
};

inline rci_t operator*(wi_t wi, Radix_t const& r1)
{
  rci_t result((unsigned int)(wi.M_val * r1.M_val));
  return result;
}

inline rci_t operator*(Radix_t const& r1, wi_t wi)
{
  rci_t result((unsigned int)(wi.M_val * r1.M_val));
  return result;
}

inline wi_t operator/(rci_t const& rci, Radix_t const& r1)
{ 
  wi_t result;
  result.M_val = rci.val() / r1.M_val;
  return result;
}

class wordPtr;

class word
{
  private:
    bool M_initialized;
    uint64_t M_word;

  public:
    // Default constructor. Construct uninitialized word.
    word(void) : M_initialized(false), M_word(0xdead12344321deadUL) { }
    // Construct a zeroed word from the int 0.
    word(int value) : M_initialized(true), M_word(0) { assert(value == 0); }
    // Construct a word from a given uint64_t integer value.
    explicit word(uint64_t value) : M_initialized(true), M_word(value) { }

    // Copy constructor.
    word(word const& w) : M_initialized(w.M_initialized), M_word(w.M_word) { assert(M_initialized); }
    // Destructor.
    ~word() { M_initialized = false; M_word = 0xdeaddeaddeaddeadUL; }

    // Assignment operators.
    word& operator=(word const& w) { assert(w.M_initialized); M_initialized = w.M_initialized; M_word = w.M_word;  return *this; }
    // Assign 0 to a word.
    word& operator=(int value)
    {
      assert(value == 0);			// Only 0 may be assigned.
      M_initialized = true;
      M_word = 0;
      return *this;
    }

    // Compare two words.
    friend bool operator==(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return w1.M_word == w2.M_word; }
    friend bool operator!=(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return w1.M_word != w2.M_word; }

    // Invert all bits in a word.
    word operator~(void) const { return word(~M_word); }

    // Convert word as boolean to a mask with all zeroes (false) or all ones (true), by negating it.
    word operator-(void) const
    {
      assert((M_word & ~1UL) == 0);
      return word(-M_word);
    }

    // Bit-wise binary operators.
    friend word operator^(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return word(w1.M_word ^ w2.M_word); }
    friend word operator&(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return word(w1.M_word & w2.M_word); }
    friend word operator|(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return word(w1.M_word | w2.M_word); }
    word& operator^=(word const& w) { assert(M_initialized && w.M_initialized); M_word ^= w.M_word; return *this; }
    word& operator&=(word const& w) { assert(M_initialized && w.M_initialized); M_word &= w.M_word; return *this; }
    word& operator|=(word const& w) { assert(M_initialized && w.M_initialized); M_word |= w.M_word; return *this; }

    // Shift operators.
    friend word operator<<(word const& w, size_t shift) { assert(w.M_initialized); assert(shift < 64); return word(w.M_word << shift); }
    friend word operator<<(word const& w, int shift) { assert(w.M_initialized); assert(shift >= 0 && shift < 64); return word(w.M_word << shift); }
    friend word operator<<(word const& w, rci_t rci) { assert(w.M_initialized); int shift = rci.val(); assert(shift >= 0 && shift < 64); return word(w.M_word << shift); }
    friend word operator>>(word const& w, size_t shift) { assert(w.M_initialized); assert(shift < 64); return word(w.M_word >> shift); }
    friend word operator>>(word const& w, int shift) { assert(w.M_initialized); assert(shift >= 0 && shift < 64); return word(w.M_word >> shift); }
    friend word operator>>(word const& w, rci_t rci) { assert(w.M_initialized); int shift = rci.val(); assert(shift >= 0 && shift < 64); return word(w.M_word >> shift); }
    word& operator<<=(int shift) { assert(M_initialized); assert(shift >= 0 && shift < 64); M_word <<= shift; return *this; }
    word& operator<<=(rci_t rci) { assert(M_initialized); int shift = rci.val(); assert(shift >= 0 && shift < 64); M_word <<= shift; return *this; }
    word& operator>>=(int shift) { assert(M_initialized); assert(shift >= 0 && shift < 64); M_word >>= shift; return *this; }
    word& operator>>=(rci_t rci) { assert(M_initialized); int shift = rci.val(); assert(shift >= 0 && shift < 64); M_word >>= shift; return *this; }

    // Initialize an array of words with zero.
    static void init_array(wordPtr const& a, wi_t size);

    // Perform explicit conversions.
    BIT convert_to_BIT(void) const
    {
      assert(M_initialized);
      assert((M_word & ~1UL) == 0);			// May only be 0 or 1.
      return M_word;
    }
    int convert_to_int(void) const
    {
      assert(M_initialized);
      assert(M_word <= 0x7fffffffU);			// Make sure the value doesn't exceed the maximum value of an int.
      return M_word;
    }
    uint64_t convert_to_uint64_t(void) const { assert(M_initialized); return M_word; }

    // NOT operator. Returns true if all bits are zero.
    bool operator!(void) const { return !M_word; }
    // Automatic conversion to boolean.
    operator bool(void) const { assert(M_initialized); return M_word != 0; }

  private:
    // Disallow conversion to int (this protects us from accidental conversion to bool (see above) and from there to int without us noticing that).
    operator int(void) const { assert(false); return 0; }
};

#define CONVERT_TO_BIT(w) ((w).convert_to_BIT())
#define CONVERT_TO_INT(w) ((w).convert_to_int())
#define CONVERT_TO_UINT64_T(w) ((w).convert_to_uint64_t())
#define CONVERT_TO_WORD(i) word((uint64_t)(i))

#include "parity.h"

class wordPtr {
  private:
    word* M_val;

  public:
    wordPtr(void) : M_val(NULL) { }
    wordPtr(wordPtr const& wp) : M_val(wp.M_val) { }
    wordPtr(void* wp) : M_val((word*)wp) { }
    ~wordPtr() { }
    wordPtr& operator=(wordPtr const& wp) { M_val = wp.M_val; return *this; }
    friend bool operator==(wordPtr const& wp1, wordPtr const& wp2) { return wp1.M_val == wp2.M_val; }
    friend bool operator!=(wordPtr const& wp1, wordPtr const& wp2) { return wp1.M_val == wp2.M_val; }
    friend bool operator<(wordPtr const& wp1, wordPtr const& wp2) { return wp1.M_val < wp2.M_val; }
    friend bool operator<=(wordPtr const& wp1, wordPtr const& wp2) { return wp1.M_val <= wp2.M_val; }
    friend bool operator>(wordPtr const& wp1, wordPtr const& wp2) { return wp1.M_val > wp2.M_val; }
    friend bool operator>=(wordPtr const& wp1, wordPtr const& wp2) { return wp1.M_val >= wp2.M_val; }
    friend wordPtr operator+(wordPtr const& wp, wi_t wi) { wordPtr result(wp); result.M_val += wi.val(); return result; }
    friend wordPtr operator-(wordPtr const& wp, wi_t wi) { wordPtr result(wp); result.M_val -= wi.val(); return result; }
    //friend wordPtr operator+(wordPtr const& wp, int wi) { wordPtr result(wp); result.M_val += wi; return result; }
    //friend wordPtr operator-(wordPtr const& wp, int wi) { wordPtr result(wp); result.M_val -= wi; return result; }
    word parity64(wordPtr buf) { return ::parity64(buf.M_val); }
    word& operator*() const { return *M_val; }
    word& operator[](wi_t const& wi) const { return M_val[wi.val()]; }
    wordPtr& operator++() { ++M_val; return *this; }
    wordPtr operator++(int) { wordPtr cur(*this); ++M_val; return cur; }
    wordPtr& operator--() { --M_val; return *this; }
    wordPtr operator--(int) { wordPtr cur(*this); --M_val; return cur; }
    wordPtr& operator+=(wi_t const& wi) { M_val += wi.val(); return *this; }
    wordPtr& operator-=(wi_t const& wi) { M_val -= wi.val(); return *this; }

    word* val(void) const { return M_val; }
};

inline void word::init_array(wordPtr const& a, wi_t size)
{
  for (wi_t i = 0U; i < size; ++i)
    a[i] = 0;
}

typedef wordPtr wordConstPtr;

class wordPtrPtr {
  private:
    wordPtr* M_val;

  public:
    wordPtrPtr(void) : M_val(NULL) { }
    wordPtrPtr(wordPtrPtr const& wpp) : M_val(wpp.M_val) { }
    ~wordPtrPtr() { }
    wordPtrPtr& operator=(wordPtrPtr& wpp) { M_val = wpp.M_val; return *this; }
    wordPtrPtr& operator=(wordPtr* wpp) { M_val = wpp; return *this; }
    friend bool operator==(wordPtrPtr const& wpp1, wordPtrPtr const& wpp2) { return wpp1.M_val == wpp2.M_val; }
    friend bool operator!=(wordPtrPtr const& wpp1, wordPtrPtr const& wpp2) { return wpp1.M_val == wpp2.M_val; }
    friend bool operator<(wordPtrPtr const& wpp1, wordPtrPtr const& wpp2) { return wpp1.M_val < wpp2.M_val; }
    friend bool operator<=(wordPtrPtr const& wpp1, wordPtrPtr const& wpp2) { return wpp1.M_val <= wpp2.M_val; }
    friend bool operator>(wordPtrPtr const& wpp1, wordPtrPtr const& wpp2) { return wpp1.M_val > wpp2.M_val; }
    friend bool operator>=(wordPtrPtr const& wpp1, wordPtrPtr const& wpp2) { return wpp1.M_val >= wpp2.M_val; }
    wordPtr& operator*() const { return *M_val; }
    wordPtr& operator[](rci_t const& rci) const { return M_val[rci.val()]; }
    wordPtrPtr& operator++() { ++M_val; return *this; }
    wordPtrPtr operator++(int) { wordPtrPtr cur(*this); ++M_val; return cur; }
    wordPtrPtr& operator--() { --M_val; return *this; }
    wordPtrPtr operator--(int) { wordPtrPtr cur(*this); --M_val; return cur; }

    wordPtr* val(void) const { return M_val; }
};

