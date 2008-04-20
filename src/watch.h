#ifndef WATCH_H
#define WATCH_H

/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org> 
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

#include <time.h>
#include <math.h>

extern clock_t watchstart;
extern clock_t watchstop;
extern long watchrunning;

extern double logtotal, logtotalsquare, logcount;

void m4ri_watch_seed();

double m4ri_sqrt(double n);

void m4ri_watch_start();

void m4ri_watch_stop();

clock_t m4ri_watch_get();

void m4ri_watch_clear_logs();

void m4ri_watch_store(clock_t watch);

double m4ri_watch_get_average_time();

double m4ri_watch_get_sigma();

long m4ri_watch_get_count();

#endif //WATCH_H



