/**
 * \file m4ri.h
 * \brief Main include file for the M4RI library.
 * 
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */
/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
*       Copyright (C) 2007,2008 Martin Albrecht <malb@informatik.uni-bremen.de> 
*
*  Distributed under the terms of the GNU General Public License (GPL)
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
#ifndef PERMUTATION_H
#define PERMUTATION_H

typedef struct {
  size_t *values;
  size_t size;
} permutation;

#endif //PERMUTATION_H
