/*******************************************************************
*
*                M4RI: Linear Algebra over GF(2)
*
*    Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include "m4ri_config.h"

#if __M4RI_HAVE_LIBPNG
#include <png.h>
#endif //__M4RI_HAVE_LIBPNG


#include "io.h"
#include "echelonform.h"

void mzd_info(const mzd_t *A, int do_rank) {
  printf("nrows: %6zu, ncols: %6zu, density: %6.5f, hash: 0x%016zx",(size_t)A->nrows,(size_t)A->ncols,mzd_density(A,1),mzd_hash(A));
  if(do_rank) {
    mzd_t *AA = mzd_copy(NULL, A);
    printf(", rank: %6zu\n",(size_t)mzd_echelonize(AA,0));
    mzd_free(AA);
  } else {
    printf("\n");
  }
}

#define SAFECHAR (m4ri_radix + m4ri_radix / 4 + 1)

void mzd_print(mzd_t const *M) {
  char temp[SAFECHAR];
  for (rci_t i = 0; i < M->nrows; ++i) {
    printf("[");
    word *row = M->rows[i];
    for (wi_t j = 0; j < M->width - 1; ++j) {
      m4ri_word_to_str(temp, row[j], 1);
      printf("%s|", temp);
    }
    row = row + M->width - 1;
    int const wide = (M->ncols % m4ri_radix) ? M->ncols % m4ri_radix : m4ri_radix;
    for (int j = 0; j < wide; ++j) {
      if(j != 0 && (j % 4) == 0)
        printf(":");
      if (__M4RI_GET_BIT(*row, j)) 
        printf("1");
      else
        printf(" ");
    }
    printf("]\n");
  }
}

void mzd_print_row(mzd_t const *M, const rci_t i) {
  char temp[SAFECHAR];
  printf("[");
  word *row = M->rows[i];
  for (wi_t j = 0; j < M->width - 1; ++j) {
    m4ri_word_to_str(temp, row[j], 1);
    printf("%s|", temp);
  }
  row = row + M->width - 1;
  int const wide = (M->ncols % m4ri_radix) ? M->ncols % m4ri_radix : m4ri_radix;
  for (int j = 0; j < wide; ++j) {
    if(j != 0 && (j % 4) == 0)
      printf(":");
    if (__M4RI_GET_BIT(*row, j)) 
      printf("1");
    else
      printf(" ");
  }
  printf("]\n");
}


#if __M4RI_HAVE_LIBPNG
#define PNGSIGSIZE 8

mzd_t * mzd_from_png(const char *fn, int verbose) {
  int retval = 0;
  mzd_t *A = NULL;
  png_byte pngsig[PNGSIGSIZE];

  FILE *fh = fopen(fn,"rb");

  if (!fh) {
    if (verbose)
      printf("Could not open file '%s' for reading\n",fn);
    return NULL;
  };
    
  if (fread((char*)pngsig, PNGSIGSIZE, 1, fh) != 1) {
    if (verbose)
      printf("Could not read file '%s'\n",fn);
    retval = 1;
    goto from_png_close_fh;
  }

  if (png_sig_cmp(pngsig, 0, PNGSIGSIZE) != 0) {
    if (verbose)
      printf("'%s' is not a PNG file.\n",fn);
    retval = 2;
    goto from_png_close_fh;
  }

  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    if (verbose)
      printf("failed to initialise PNG read struct.\n");
    retval = 3;
    goto from_png_close_fh;
  }
  png_set_user_limits(png_ptr, 0x7fffffffL,  0x7fffffffL);

  png_infop info_ptr = png_create_info_struct(png_ptr);

  if (!info_ptr) {
    if (verbose)
      printf("failed to initialise PNG info struct\n");
    retval = 3;
    goto from_png_destroy_read_struct;
  }

  png_init_io(png_ptr, fh);
  png_set_sig_bytes(png_ptr, PNGSIGSIZE);
  png_read_info(png_ptr, info_ptr);
 
  const png_uint_32 m                = png_get_image_height(png_ptr, info_ptr);
  const png_uint_32 n                = png_get_image_width(png_ptr, info_ptr);
  const png_uint_32 bit_depth        = png_get_bit_depth(png_ptr, info_ptr);
  const png_uint_32 channels         = png_get_channels(png_ptr, info_ptr);
  const png_uint_32 color_type       = png_get_color_type(png_ptr, info_ptr);
  const png_uint_32 compression_type = png_get_compression_type(png_ptr, info_ptr);
  const png_uint_32 interlace_type   = png_get_interlace_type(png_ptr, info_ptr);

  if (interlace_type != PNG_INTERLACE_NONE) {
    if (verbose)
      printf("interlaced images not supported\n");
    goto from_png_destroy_read_struct;
  };

  if (verbose)
    printf("reading %u x %u matrix (bit depth: %u, channels: %u, color type: %u, compression type: %u)\n",(unsigned int)m,
           (unsigned int)n, (unsigned int)bit_depth, (unsigned int)channels, (unsigned int)color_type, (unsigned int)compression_type);

  if(color_type != 0 && color_type != 3) {
    if (verbose)
      printf("only graycscale and palette colors are supported.\n");
    goto from_png_destroy_read_struct;
  }
    
  A = mzd_init(m, n);
  const word bitmask_end = A->high_bitmask;
  png_bytep row = m4ri_mm_calloc(sizeof(char),n/8+1);

  word tmp;
  wi_t j;

  png_set_packswap(png_ptr);
  //png_set_invert_mono(png_ptr);

  for(rci_t i=0; i<m; i++) {
    png_read_row(png_ptr, row, NULL);
    for(j=0; j<A->width-1; j++) {
      tmp = ((word)row[8*j+7])<<56 | ((word)row[8*j+6])<<48 \
        |   ((word)row[8*j+5])<<40 | ((word)row[8*j+4])<<32 \
        |   ((word)row[8*j+3])<<24 | ((word)row[8*j+2])<<16 \
        |   ((word)row[8*j+1])<< 8 | ((word)row[8*j+0])<< 0;
      A->rows[i][j] = ~tmp;
    }
    tmp = 0;
    switch((n/8 + ((n%8) ? 1 : 0))%8) {
    case 0: tmp |= ((word)row[8*j+7])<<56;
    case 7: tmp |= ((word)row[8*j+6])<<48;
    case 6: tmp |= ((word)row[8*j+5])<<40;
    case 5: tmp |= ((word)row[8*j+4])<<32;
    case 4: tmp |= ((word)row[8*j+3])<<24;
    case 3: tmp |= ((word)row[8*j+2])<<16;
    case 2: tmp |= ((word)row[8*j+1])<< 8;
    case 1: tmp |= ((word)row[8*j+0])<< 0;
    };
    A->rows[i][j] |= (~tmp & bitmask_end);
  }

  m4ri_mm_free(row);
  png_read_end(png_ptr, NULL);

 from_png_destroy_read_struct: 
  png_destroy_read_struct(&png_ptr, &info_ptr,(png_infopp)0);

 from_png_close_fh: 
  fclose(fh);

  if (retval != 0 && A) {
    mzd_free(A);
    return NULL;
  } else {
    return A;
  }
}

int mzd_to_png(const mzd_t *A, const char *fn, int compression_level, const char *comment, int verbose) {
  FILE *fh = fopen(fn, "wb");

  if (!fh) {
    if(verbose)
      printf("Could not open file '%s' for writing\n",fn);
    return 1;
  }

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    if(verbose)
      printf("failed to initialise PNG write struct.\n");
    fclose(fh);
    return 3;
  }
  png_set_user_limits(png_ptr, 0x7fffffffL,  0x7fffffffL);

  png_infop info_ptr = png_create_info_struct(png_ptr);

  if (!info_ptr) {
    if (verbose)
      printf("failed to initialise PNG info struct\n");
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fh);
    return 3;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    if (verbose)
      printf("error writing PNG file\n");
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fh);
    return 1;
  }

  png_init_io(png_ptr, fh);
  png_set_compression_level(png_ptr, compression_level);

  png_set_IHDR(png_ptr, info_ptr, A->ncols, A->nrows, 1, \
               PNG_COLOR_TYPE_GRAY, \
               PNG_INTERLACE_NONE, \
               PNG_COMPRESSION_TYPE_DEFAULT, \
               PNG_FILTER_TYPE_DEFAULT);

  png_text txt_ptr[3];

  char pdate[21];
  time_t ptime=time(NULL);
  struct tm *ltime=localtime(&ptime);
  sprintf(pdate,"%04d/%02d/%02d %02d:%02d:%02d",ltime->tm_year+1900,ltime->tm_mon+1,ltime->tm_mday,ltime->tm_hour,ltime->tm_min,ltime->tm_sec);

  txt_ptr[0].key="Software";
  txt_ptr[0].text="M4RI";
  txt_ptr[0].compression=PNG_TEXT_COMPRESSION_NONE;
  txt_ptr[1].key="Date";
  txt_ptr[1].text=pdate;
  txt_ptr[1].compression=PNG_TEXT_COMPRESSION_NONE;
  txt_ptr[2].key="Comment";
  txt_ptr[2].text=(char*)comment;
  txt_ptr[2].compression=PNG_TEXT_COMPRESSION_NONE;

  png_set_text(png_ptr, info_ptr, txt_ptr, 3);

  png_write_info(png_ptr, info_ptr);

  png_set_packswap(png_ptr);
  png_set_invert_mono(png_ptr);
  
  png_bytep row = m4ri_mm_calloc(sizeof(char),A->ncols/8+8);

  wi_t j=0;
  word tmp = 0;
  for(rci_t i=0; i<A->nrows; i++) {
    word *rowptr = A->rows[i];
    for(j=0; j<A->width-1; j++) {
      tmp = rowptr[j];
      row[8*j+0] = (png_byte)((tmp>> 0) & 0xff);
      row[8*j+1] = (png_byte)((tmp>> 8) & 0xff);
      row[8*j+2] = (png_byte)((tmp>>16) & 0xff);
      row[8*j+3] = (png_byte)((tmp>>24) & 0xff);
      row[8*j+4] = (png_byte)((tmp>>32) & 0xff);
      row[8*j+5] = (png_byte)((tmp>>40) & 0xff);
      row[8*j+6] = (png_byte)((tmp>>48) & 0xff);
      row[8*j+7] = (png_byte)((tmp>>56) & 0xff);
    }
    tmp = rowptr[j];
    switch( (A->ncols/8 + ((A->ncols%8) ? 1 : 0)) %8 ) {
    case 0: row[8*j+7] = (png_byte)((tmp>>56) & 0xff);
    case 7: row[8*j+6] = (png_byte)((tmp>>48) & 0xff);
    case 6: row[8*j+5] = (png_byte)((tmp>>40) & 0xff); 
    case 5: row[8*j+4] = (png_byte)((tmp>>32) & 0xff); 
    case 4: row[8*j+3] = (png_byte)((tmp>>24) & 0xff);
    case 3: row[8*j+2] = (png_byte)((tmp>>16) & 0xff);
    case 2: row[8*j+1] = (png_byte)((tmp>> 8) & 0xff);
    case 1: row[8*j+0] = (png_byte)((tmp>> 0) & 0xff);
    };
    png_write_row(png_ptr, row);
  }
  m4ri_mm_free(row);

  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fh);
  return 0;
}

#endif //__M4RI_HAVE_LIBPNG

mzd_t *mzd_from_jcf(const char *fn, int verbose) {
  int retval = 0;
  mzd_t *A = NULL;
  FILE *fh = fopen(fn,"r");

  rci_t m,n;
  long p  = 0;
  long nonzero = 0;

  if (!fh) {
    if (verbose)
      printf("Could not open file '%s' for reading\n",fn);
    return NULL;
  }

  if (fscanf(fh, "%d %d %ld\n%ld\n\n",&m,&n,&p,&nonzero) != 4) {
    if (verbose)
      printf("File '%s' does not seem to be in JCF format.",fn);
    retval = 1;
    goto from_jcf_close_fh;
  }

  if(p != 2) {
    if (verbose)
      printf("Expected p==2 but found p==%ld\n",p);
    retval = 1;
    goto from_jcf_close_fh;
  }

  if (verbose)
    printf("reading %lu x %lu matrix with at most %ld non-zero entries (density at most: %6.5f)\n",
           (unsigned long)m, (unsigned long)n, (unsigned long)nonzero, ((double)nonzero)/((double)m*n));

  A = mzd_init(m,n);
  
  long i = -1;
  long j = 0;

  while(fscanf(fh,"%ld\n",&j) == 1) {
    if (j<0) {
      i++, j = -j;
    }
    if (((j-1) >= n) || (i>= m))
      m4ri_die("trying to write to (%ld,%ld) in %ld x %ld matrix\n", i, j-1, m, n);
    mzd_write_bit(A, i, j-1, 1);
  };

 from_jcf_close_fh:
  fclose(fh);

  if(retval != 0 && A) {
    mzd_free(A);
    return NULL;
  } else {
    return A;
  }
}

mzd_t *mzd_from_str(rci_t m, rci_t n, const char *str) {
  int idx = 0;
  mzd_t *A = mzd_init(m, n);
  for(rci_t i=0; i<A->nrows; i++) {
    for(rci_t j=0; j<A->ncols; j++) {
      mzd_write_bit(A, i, j, str[idx++] == '1');
    }
  }
  return A;
}
