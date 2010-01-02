/* ---------------------------------------------------------------------------

simple_ga.h header file of simple_ga.cpp, 
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


USAGE: You will most likely only want to use the class Simple_Popul. 
This uses Simple_Invid. 

                             Create a new Population 
                             with standard parameters defined in ga.h
Simple_Popul P;
P.Init(&<your_fitness_function>,<your_chromosome_size>);
                             The fitness function always takes the form
			     double <your_fitness_function> (char* chrom, unsigned int sz);
			     chrom is the chromosome string
			     sz is the size of that string
			     remember it is a function pointer

                             Change any ga parameters:
P.Pm(float <your_pm>);       Set mutation probability to your preferred value
P.Pc(float <your_pc>);       Set crossover probability to your preferred value
P.Pa(unsigned int <your_Pa); Change the number of parents to your preferred number
                             (only even numbers are permitted)
P.Ch(unsigned int <your_Ch); Change the number of children to your preferred number
                             (only even numbers are permitted)

                             Initialize population with random values
P.Randomise();
                             Initialize one individual with values from file
P.Randomise(unsigned int <invid>, char* <filename>);

                             Perform evolution for a number of generations
P.Generation(unsigned int <your_#_generations>);

                             Save and load entire population (binary)
P.Save(ofstream* <filedescriptor>);
P.Load(ifstream* <filedescriptor>);
                             Write entire population into text file
P.Write(ofstream* <filedescriptor>);

P.Maximise();                Do maximisation (long execution times are favored)
P.Minimise();                Do minimisation (short execution times are favored)

In the fitness function, the string of chromosome char* must be transferred 
into the data types that are used in the tested program. chrom is a pointer 
to a character array, so you have to cast the char-pointer into a 
type-x-pointer and then access the value in there, for example:

double Fitness_Function (char* chrom, unsigned int size) {
  int x;
  int y;
  x = *((int*)chrom+0); // cast char-pointer to int-pointer and read value 
  y = *((int*)chrom+sizeof(int));
  // timing start
  tested_function(int x, int y);
  // timing stop
}
------------------------------------------------------------------------------ */

#include <iostream>
#include <fstream>

#ifndef __SIMPLE_GA_H__
#define __SIMPLE_GA_H__

#undef SS

using namespace std;

class Simple_Invid;
class Simple_Popul;

void bubble_sort_max(Simple_Invid** List, unsigned int size);
void bubble_sort_min(Simple_Invid** List, unsigned int size);
void survivalA(Simple_Invid** list, unsigned int psize, unsigned int csize);
unsigned int Tournament_Selection(Simple_Invid** List, unsigned int size,
		unsigned int tour);

class Simple_Invid { // DECLARATION of the Individual

protected:
	unsigned char* chr; // pointer to chromosome
	unsigned int csz; // size of the chromosome
	double gaf; // fitness for the ga process
	int approxLvl; // GACov

public:
	friend class Simple_Popul;
	Simple_Invid (void); // constructor
	~Simple_Invid (void); // destructor
	void Init (unsigned int); // initialiser

	Simple_Invid& operator= (Simple_Invid&);

	// Small rewrite to accomodate for approximation level
	int operator< (Simple_Invid& I) {
		if(approxLvl > I.approxLvl)
		return 1;
		else if(approxLvl == I.approxLvl)
		return gaf < I.gaf;
		else
		return 0;
	};

	// Small rewrite to accomodate for approximation level
	int operator<= (Simple_Invid& I) {//return gaf <= I.gaf;};
		return approxLvl < I.approxLvl || (approxLvl == I.approxLvl && gaf < I.gaf);
	};

	int operator> (Simple_Invid& I) {return gaf > I.gaf;};
	int operator>= (Simple_Invid& I) {return gaf >= I.gaf;};
	int operator== (Simple_Invid& I) {return approxLvl == I.approxLvl && gaf == I.gaf;};
	int operator!= (Simple_Invid& I) {return gaf != I.gaf;};

	void Randomise (double (*FF) (unsigned char*, unsigned int));
	void Randomise (double (*FF) (unsigned char*, unsigned int), char*);
	void Crossover (Simple_Invid&, float);
	void Mutation (float);
	double Fitness (double (*FF) (unsigned char*, unsigned int));

	double Gaf (void) {return gaf;}; // return ga fitness

	int Save (ofstream&); // save invid to ofstream
	int Load (ifstream&); // load invid from ifstream

	// GACov
	unsigned char* getChromData();
	void setApproxLvl(int al);
	int getApproxLvl();

private:
	// may be template and usable for inherited classes
	friend ostream& operator<< (ostream& out, Simple_Invid& I) {
		out << "(" << I.approxLvl << ")" << I.gaf << " ";
		//out << (int)(I.chr[0]) << "-" << (int)(I.chr[1]) << "-" << (int)(I.chr[2]) << " (" << I.approxLvl << ") " << I.gaf << endl;
		return out;
	};
};

class Simple_Popul { // DECLARATION of the Population

protected:

	Simple_Invid* pop; // list of individuals
	Simple_Invid** ppop; // pointerlist of invids
	Simple_Invid** cpop; // pointerlist of children

	unsigned int pasz; // number of parents
	unsigned int chsz; // number of children
	unsigned int psz; // total number of invids
	unsigned int csz; // size of the invids' chromosome

	float pc; // crossover probability
	float pm; // mutation probability

	unsigned int tour; // for tournament selection

	double (*FF) (unsigned char*, unsigned int); // pointer to fitness function

	void (*SS) (Simple_Invid**, unsigned int);
	// pointer to sorting algorithm

	// GACov
	int* checkedBranches;
	int branches;
	int looking;
	int branch;
	int nodeExecuted;
	double fitness;
	int globalAL;
	Simple_Invid* branchData;
	Simple_Invid currentInvid;

public: // basic functionality

	Simple_Popul (void); // constructor
	~Simple_Popul (void); // destructor
	void Init (double (*f) (unsigned char*, unsigned int), unsigned int size);
	// Initializer for the population: assign fitness function,
	// and assign the size of the chromosome (unsigned int size)

	void Randomise (void); // Initialize chroms randomly
	void Randomise (unsigned int, char*); // Initialize chroms
	void Generation (unsigned int number); // Perform number of generations

	int Save (ofstream&); // save popul to ofstream binary
	int Load (ifstream&); // load popul from ifstream binary
	int Write(ofstream&); // write chrom's to ofstream text

public: // change population parameter

	float Pm (void) {return pm;}; // return probabilities
	float Pc (void) {return pc;};
	void Pm (float p) {pm = p;}; // set probabilities
	void Pc (float p) {pc = p;};

	int Pa (void) {return pasz;}; // return size of pa and ch
	int Ch (void) {return chsz;};

	void Pa (unsigned int); // set size of pa and ch
	void Ch (unsigned int);

	// Optimisation strategy
	void Maximise (void) {survive = 0; SS = bubble_sort_max;};// Maximisation problem
	void Minimise (void) {survive = 0; SS = bubble_sort_min;};// Minimisation problem


	// GACov
	void setBranches(int b);
	void checkBranch(int b);
	void checkBranchLoop(int b, int its, int reqIts);
	double lookingBranch(int b, int pred1, int pred2);
	double lookingBranchAL(int b, int al, int pred1, int pred2);
	double lookingBranchDouble(int b, double pred1, double pred2);
	double lookingBranchDoubleAL(int b, int al, double pred1, double pred2);
	double lookingBranchHammDist(int b, int pred1, int pred2);
	double lookingBranchHammDistAL(int b, int al, int pred1, int pred2);
	double lookingBranchLoop(int b, int its, int reqIts);
	double lookingBranchLoopAL(int b, int al, int its, int reqIts);
	int* getCheckedBranches();
	Simple_Invid* getBranchData();
	double getFitness();

	void (*survive) (Simple_Invid**, unsigned int, unsigned int);

	// only use when survivalA is to be used!
	// USE WITH CAUTION
	void setSurvivalA() {SS = 0; survive = survivalA;};

private:
	// may be template and usable for inherited classes
	friend ostream& operator<< (ostream& out, Simple_Popul& I);
	friend void survivalA(Simple_Invid**, unsigned int);
	friend void bubble_sort_max (Simple_Invid**, unsigned int);
	friend void bubble_sort_min (Simple_Invid**, unsigned int);
	friend unsigned int Tournament_Selection
	(Simple_Invid**, unsigned int, unsigned int);
};

#endif
