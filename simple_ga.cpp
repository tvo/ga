/* ---------------------------------------------------------------------------

simple_ga.cpp
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <time.h>

#include "ga.h"
#include "simple_ga.h"

// ----------------------------------------------------------------
//
//             SIMPLE_GA INDIVIDUAL --- IMPLEMENTATION
//
// ----------------------------------------------------------------

// ---------- constructor
Simple_Invid::Simple_Invid(void) {
	chr = NULL;
	csz = 0;
	gaf = dflt_fitness;
}

// ---------- initializer
void Simple_Invid::Init(unsigned int size) {
	csz = size; // set size
	chr = new unsigned char[csz]; // create a chrom for invid
}

// ---------- destructor
Simple_Invid::~Simple_Invid(void) {
	delete[] chr; // delete the chrom of invid
}

// ---------- copy operator
Simple_Invid& Simple_Invid::operator=(Simple_Invid& I) {
	if (this != &I) {
		delete[] chr; // delete TO chrom
		csz = I.csz; // copy size from FROM
		gaf = I.gaf; // copy fitness from FROM
		chr = new unsigned char[csz]; // create new chrom
		for (unsigned int i = 0; i < csz; i++)
			*(chr + i) = *(I.chr + i); // copy chrom values from FROM into TO
	}
	return *this;
}

// ---------- initiate chrom
void Simple_Invid::Randomise(double(*FF)(unsigned char*, unsigned int)) {
	for (unsigned int i = 0; i < csz; i++)
		chr[i] = (char) rand(); // initialize chrom randomly
	if (FF != NULL) {
		gaf = (*FF)(chr, csz); // invoke fitness function
	}
}

// ---------- initiate with data from file
void Simple_Invid::Randomise(double(*FF)(unsigned char*, unsigned int), char* fname) {
	unsigned int counter = 0;
	ifstream ins;
	ins.open(fname); // read chrom from (binary) file
	while ((ins.read((char*) chr + counter, 1)) && (counter < csz))
		counter++;
	if (FF != NULL)
		gaf = (*FF)(chr, csz); // invoke fitness function
	ins.close();
}

// ---------- mutate chrom
void Simple_Invid::Mutation(float pm) {
	static char mask;
	static unsigned int i, j;
	for (i = 0; i < csz; i++) { // for each byte in the variable
		mask = 0; // set mask to zero
		for (j = 0; j < 8; j++) { // put a 1 into each randomly
			mask <<= 1; // selected bit to determine
			if (drand48() < pm) // which bit is swapped
				mask += 1;
		}
		chr[i] ^= mask; // apply mask to every chrom byte
	}
}

// ---------- recombination of two invids
void Simple_Invid::Crossover(Simple_Invid& I, float pc) {
	if (this != &I) {
		char mask_xor, mask_and;
		unsigned int i, j;
		for (i = 0; i < csz; i++) { // for each byte in the chrom
			mask_and = 0; // set and_mask to zero
			for (j = 0; j < 8; j++) {
				mask_and <<= 1; // put a 1 into each randomly
				if (drand48() < pc) // determine crossover locations
					mask_and += 1; // selected positions where
			} // bits are swapped
			mask_xor = chr[i] ^ I.chr[i];
			mask_xor &= ~mask_and; // swap the bits
			chr[i] ^= mask_xor;
			I.chr[i] ^= mask_xor;
		}
	}
}

// ---------- call the fitness function
double Simple_Invid::Fitness(double(*FF)(unsigned char*, unsigned int)) {
	gaf = (*FF)(chr, csz);
	return gaf;
}

// ---------- save invid to ofstream
int Simple_Invid::Save(ofstream& fdout) {
	fdout.write((char*) &gaf, sizeof(gaf));
	fdout.write((char*) &csz, sizeof(csz));
	fdout.write((char*) &chr, csz);
	return 0;
}

// ---------- load invid from ifstream        
int Simple_Invid::Load(ifstream& fdin) {
	fdin.read((char*) &gaf, sizeof(gaf));
	fdin.read((char*) &csz, sizeof(csz));
	fdin.read((char*) &chr, csz);
	return 0;
}

// -----------------------------------------------------------------
//
//             SIMPLE_GA POPULATION --- IMPLEMENTATION
//
// -----------------------------------------------------------------

// ---------- constructor
Simple_Popul::Simple_Popul(void) {
	SS = bubble_sort_max; // default is maximisation
	pc = dflt_crossover_rate;
	pm = dflt_mutation_rate;
	tour = dflt_tour_size;
	pasz = dflt_pa_size; // number of parents
	chsz = dflt_ch_size; // number of children
	psz = pasz + chsz; // total number of invids
	pop = new Simple_Invid[psz]; // create new population
	ppop = new Simple_Invid*[psz]; // create pointerlist
	// connect pointerlist to the population
	for (unsigned int i = 0; i < psz; i++) {
		ppop[i] = &(pop[i]);
	}
}

// ---------- initialiser
void Simple_Popul::Init(double(*f)(unsigned char*, unsigned int), unsigned int size) {
	FF = f;
	csz = size;
	for (unsigned int i = 0; i < psz; i++) {
		ppop[i]->Init(csz); // init individuals using pointers
	}
	cpop = ppop + pasz; // init pointerlist for children
}

// ---------- destructor
Simple_Popul::~Simple_Popul(void) {
	delete[] pop;
	delete[] ppop;
}

// ---------- initiate individuals in the population
void Simple_Popul::Randomise(void) {
	static bool seeded = false;
	if (!seeded) {
		seeded = true;
		srand(time(NULL)); // init random generator (Windows)
#ifndef WIN32
		srand48(time(NULL)); // Linux uses different rand init
#endif
	}
	for (unsigned int i = 0; i < pasz + chsz; i++)
		ppop[i]->Randomise(FF);
	SS(ppop, psz); // sort the population min or max
}

// ---------- Overwrite Invid <number> (parents only) with data from file 
void Simple_Popul::Randomise(unsigned int number, char* fname) {
	if ((number >= 0) && (number <= pasz)) {
		ppop[number]->Randomise(FF, fname);
	}
}

// ---------- perform <number> generations
void Simple_Popul::Generation(unsigned int number) {
	// for number generations
	for (unsigned int loop = 0; loop < number; loop++) {
		// for each child
		for (unsigned int child = 0; child < chsz - 1; child++) {

			// clone two parents
			int p1 = Tournament_Selection(ppop, pasz, tour);
			int p2 = Tournament_Selection(ppop, pasz, tour);
			*(cpop[child]) = *(ppop[p1]);
			*(cpop[child + 1]) = *(ppop[p2]);

			// crossover the clones
			cpop[child]->Crossover(*(cpop[child + 1]), pc);

			// mutate the two new children
			cpop[child]->Mutation(pm);
			cpop[child + 1]->Mutation(pm);

			// call fitness function for both new invids
			cpop[child]->Fitness(FF);
			cpop[child + 1]->Fitness(FF);
		}
		SS(ppop, psz); // sort population min or max
	}
}

// ---------- save popul to filestream
//            note, that the fitness is not saved (must be re-run)
int Simple_Popul::Save(ofstream& fdout) {
	fdout.write((char*) &pasz, sizeof(pasz));
	fdout.write((char*) &chsz, sizeof(chsz));
	fdout.write((char*) &psz, sizeof(psz));
	fdout.write((char*) &csz, sizeof(csz));
	fdout.write((char*) &pc, sizeof(pc));
	fdout.write((char*) &pm, sizeof(pm));
	fdout.write((char*) &tour, sizeof(tour));
	// save the invids
	for (unsigned int i = 0; i < psz; i++) {
		ppop[i]->Save(fdout);
	}
	return 0;
}

// ---------- load popul from filestream
//            note, that the fitness is not loaded (must be re-run)
int Simple_Popul::Load(ifstream& fdin) {
	fdin.read((char*) &pasz, sizeof(pasz));
	fdin.read((char*) &chsz, sizeof(chsz));
	fdin.read((char*) &psz, sizeof(psz));
	fdin.read((char*) &csz, sizeof(csz));
	fdin.read((char*) &pc, sizeof(pc));
	fdin.read((char*) &pm, sizeof(pm));
	fdin.read((char*) &tour, sizeof(tour));
	// delete actual population
	delete[] pop;
	delete[] ppop;
	pop = new Simple_Invid[psz]; // create new population
	ppop = new Simple_Invid*[psz]; // create pointerlist
	// connect pointerlist to the population
	for (unsigned int i = 0; i < psz; i++) {
		ppop[i] = &(pop[i]);
	}
	// load the invids
	for (unsigned int j = 0; j < psz; j++) {
		ppop[j]->Load(fdin);
	}
	return 0;
}

// ---------- write the popul's chrom's to a text file
//            (for cluster analysis for example)
int Simple_Popul::Write(ofstream& fdout) {
	for (unsigned int i = 0; i < pasz; i++) {
		for (unsigned int j = 0; j < csz; j++) {
			fdout << (float) (ppop[i]->chr[j]) << " ";
		}
		fdout << endl;
	}
	return 0;
}

// ---------- change the number of parents
void Simple_Popul::Pa(unsigned int size) {
	// set new population size
	psz = size + chsz;
	// create two new lists (objects and pointers)
	Simple_Invid* tmp = new Simple_Invid[psz];
	Simple_Invid** ptmp = new Simple_Invid*[psz];
	// connect pointerlist to the population
	for (unsigned int i = 0; i < psz; i++) {
		ptmp[i] = &(tmp[i]);
	}
	// initiate the new individuals by using pointerlist
	for (unsigned int j = 0; j < psz; j++) {
		ptmp[j]->Init(csz);
	}
	// copy the old individuals into the new ones
	unsigned int copy_size = Min(size, pasz) + chsz;
	for (unsigned int k = 0; k < copy_size; k++) {
		*(ptmp[k]) = *(ppop[k]); // copy the individuals into new list
	}
	delete[] pop;
	delete[] ppop;
	// set up pointers properly
	pasz = size;
	pop = tmp;
	ppop = ptmp;
	cpop = ptmp + pasz;
}

// ---------- change the number of children
void Simple_Popul::Ch(unsigned int size) {
	// set new population size
	psz = size + pasz;
	// create two new lists (objects and pointers)
	Simple_Invid* tmp = new Simple_Invid[psz];
	Simple_Invid** ptmp = new Simple_Invid*[psz];
	// connect pointerlist to the population
	for (unsigned int i = 0; i < psz; i++) {
		ptmp[i] = &(tmp[i]);
	}
	// initiate the new individuals by using pointerlist
	for (unsigned int j = 0; j < psz; j++) {
		ptmp[j]->Init(csz);
	}
	// copy the old individuals into the new ones
	unsigned int copy_size = Min(size, chsz) + pasz;
	for (unsigned int k = 0; k < copy_size; k++) {
		*(ptmp[k]) = *(ppop[k]); // copy the individuals into new list
	}
	delete[] pop;
	delete[] ppop;
	// set up pointers properly
	chsz = size;
	pop = tmp;
	ppop = ptmp;
	cpop = ptmp + pasz;
}

// ---------- Simple_Population printout
ostream& operator<<(ostream& out, Simple_Popul& p) {
	//for (unsigned int i=0; i<p.psz; i++)
	out << *(p.ppop[0]) << " ";
	out << endl;
	return out;
}

static bool compare_max(Simple_Invid* a, Simple_Invid* b) {
	return *a > *b;
}

static bool compare_min(Simple_Invid* a, Simple_Invid* b) {
	return *a < *b;
}

// ---------- FRIEND bubble sort for the population pointer list
void bubble_sort_max(Simple_Invid** List, unsigned int size) {
	std::sort(List, List + size, compare_max);
}

// ---------- FRIEND bubble sort for the population pointer list
void bubble_sort_min(Simple_Invid** List, unsigned int size) {
	std::sort(List, List + size, compare_min);
}

// ---------- FRIEND TOURNAMENT Selection Algorithm
unsigned int Tournament_Selection(Simple_Invid** List, unsigned int size, unsigned int tour) {
	unsigned int n_pos;
	unsigned int pos = (unsigned int) rand() % size;
	double fit = List[pos]->Gaf();
	for (unsigned int i = 0; i < tour; i++) {
		n_pos = (unsigned int) rand() % size;
		if (List[n_pos]->Gaf() > fit) {
			fit = List[n_pos]->Gaf();
			pos = n_pos;
		}
	}
	return pos; // return position in the population
}

