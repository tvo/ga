/* ---------------------------------------------------------------------------

   ga.h global header file of simple_ga.cpp, 
   a simple genetic algorithm implementation

   Copyright (C) 2005  Hans-Gerhard Gross, EWI TUDelft, NL

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details: write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

   ------------------------------------------------------------------------------ */

#include <stdlib.h>

// ---------- RANDOM GENERATOR drand48() used under WIN32
#ifdef WIN32 // win32 does not know drand48()
inline double drand48(void) {
    return (double) (rand() % RAND_MAX) / RAND_MAX;
    // define RAND_MAX: 2147483647 in stdlib.h in LINUX
}
#endif

// ---------- return the minimum of two values
template <class T> 
    inline T Min (T a, T b) {
	if (a < b) 
	    return a;
	else
	    return b;
    }

// CONSTANTS

const double       dflt_fitness        = (double) 0.0;   
const float        dflt_crossover_rate = (float) 0.5;  // uniform crossover
const float        dflt_mutation_rate  = (float) 0.1;
const unsigned int dflt_pa_size        = 40;  // even number! population
const unsigned int dflt_ch_size        = 40;  // even number! population         
const unsigned int dflt_tour_size      = 3;   // for tournament selection





