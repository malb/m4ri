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

static inline uint64_t reverse(uint64_t v)
{
  v = ((v >>  1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
  v = ((v >>  2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
  v = ((v >>  4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
  v = ((v >>  8) & 0x00FF00FF00FF00FFULL) | ((v & 0x00FF00FF00FF00FFULL) << 8);
  v = ((v >> 16) & 0x0000FFFF0000FFFFULL) | ((v & 0x0000FFFF0000FFFFULL) << 16);
  v =  (v >> 32)                          |  (v                          << 32);
  return v;
}

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

    // Reverse the bits of this word.
    word& reverse(void) { M_word = ::reverse(M_word); return *this; }

    // Compare two words.
    friend bool operator==(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return w1.M_word == w2.M_word; }
    friend bool operator!=(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return w1.M_word != w2.M_word; }

    // Invert all bits in a word.
    word operator~(void) const { return word(~M_word); }

    // Convert word as boolean to a mask with all zeroes (false) or all ones (true), by negating it.
    word operator-(void) const
    {
      uint64_t reversed_word = ::reverse(M_word);
      assert((reversed_word & ~1UL) == 0);
      return word(-reversed_word);
    }

    // Bit-wise binary operators.
    friend word operator^(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return word(w1.M_word ^ w2.M_word); }
    friend word operator&(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return word(w1.M_word & w2.M_word); }
    friend word operator|(word const& w1, word const& w2) { assert(w1.M_initialized && w2.M_initialized); return word(w1.M_word | w2.M_word); }
    word& operator^=(word const& w) { assert(M_initialized && w.M_initialized); M_word ^= w.M_word; return *this; }
    word& operator&=(word const& w) { assert(M_initialized && w.M_initialized); M_word &= w.M_word; return *this; }
    word& operator|=(word const& w) { assert(M_initialized && w.M_initialized); M_word |= w.M_word; return *this; }

    // Allow conversion of a word that has the integer value 2^n into a mask with n bits set, by substracting 1.
    friend word operator-(word const& w, int value)
    {
      uint64_t reversed_word = ::reverse(w.M_word);
      assert(w.M_initialized);
      assert((reversed_word & (reversed_word - 1)) == 0);	// Must be a power of 2.
      assert(value == 1);					// May only substract 1.
      return word(::reverse(reversed_word - 1));
    }

    // Shift operators.
    friend word operator<<(word const& w, size_t shift) { assert(w.M_initialized); return word(w.M_word >> shift); }
    friend word operator<<(word const& w, int shift) { assert(w.M_initialized); return word(w.M_word >> shift); }
    friend word operator>>(word const& w, size_t shift) { assert(w.M_initialized); return word(w.M_word << shift); }
    friend word operator>>(word const& w, int shift) { assert(w.M_initialized); return word(w.M_word << shift); }
    word& operator<<=(int shift) { assert(M_initialized); M_word >>= shift; return *this; }
    word& operator>>=(int shift) { assert(M_initialized); M_word <<= shift; return *this; }

    // Initialize an array of words with zero.
    static void init_array(word* a, size_t size)
    {
      for (size_t i = 0; i < size; ++i)
	a[i] = 0;
    }

    // Perform explicit conversions.
    BIT convert_to_BIT(void) const
    {
      uint64_t reversed_word = ::reverse(M_word);
      assert(M_initialized);
      assert((reversed_word & ~1UL) == 0);			// May only be 0 or 1.
      return reversed_word;
    }
    int convert_to_int(void) const
    {
      uint64_t reversed_word = ::reverse(M_word);
      assert(M_initialized);
      assert(reversed_word <= 0x7fffffffU);			// Make sure the value doesn't exceed the maximum value of an int.
      return reversed_word;
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

