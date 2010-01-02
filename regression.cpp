/*
 * regression.cpp
 *
 *  Created on: Jan 1, 2010
 *      Author: Tobi Vollebregt
 */

#include "simple_ga.h"

using std::cout;
using std::endl;

double onemax_fitness(unsigned char* chrom, unsigned int sz) {
	int c = 0;
	for (unsigned int i = 0; i < sz; ++i) {
		int g = chrom[i];
		while (g) {
			c += (g & 1);
			g >>= 1;
		}
	}
	return double(c) / double(8 * sz);
}

int run() {
	Simple_Popul P;

	// parameter setup
	int n = 32;
	P.Init(onemax_fitness, n / 8);
	P.Pm(1.0 / n); // 1 bit/chrom mutates (expected value)
	P.Pc(1.0 / n); // 1 crossover position (half of chrom from one parent, other half from other parent)
	//P.Pa(50);
	//P.Ch(50);
	P.Maximise();

	// initialize population
	P.Randomise();

	// run for a while
	double f;
	int c = 0;
	do {
		P.Generation(1);
		f = P.getFitness();
		++c;
	} while (f < 1.0);

	return c;
}

int main() {
	for (int i = 0; i < 1000; ++i) {
		cout << run() << endl;
	}
}
