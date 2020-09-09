#include "project1.hpp"
#include <iostream>
#include "time.h"
#include <iomanip>
#include <fstream>

using namespace std;

int main() {
	clock_t start, finish;

	if (false) {

		// Compares computation time of the three methods

		Poisson my_solver;
		my_solver.Initialize(3);

		double time_proper = 0.0, time_efficient = 0.0, time_LU = 0.0;
		double n = 100.0;

		start = clock();
		for (int i = 0; i < n*100; i++)
			my_solver.Integrate_proper();
		finish = clock();
		time_proper += (double)finish - start;

		start = clock();
		for (int i = 0; i < n*100; i++)
			my_solver.Integrate_efficient();
		finish = clock();
		time_efficient += (double)finish - start;

		start = clock();
		for (int i = 0; i < n; i++)
			my_solver.LU_decomp();
		finish = clock();
		time_LU += (double)finish - start;

		my_solver.Free();

		time_proper = (time_proper/n)/(double)CLOCKS_PER_SEC*10.0;
		time_efficient = (time_efficient/n)/(double)CLOCKS_PER_SEC*10.0;
		time_LU = (time_LU/n)/(double)CLOCKS_PER_SEC*1000.0;

		ofstream outFile;
		outFile.open("Time.txt");

		outFile << setprecision(8) << "time proper:  " << time_proper << " ms." << endl;
		outFile << setprecision(8) << "time efficient:  " << time_efficient << " ms." << endl;
		outFile << setprecision(8) << "time LU:  " << time_LU << " ms.";

		outFile.close();
	}


	if (false) {

		// Finds the max error of different Ns, and writes them to file

		ofstream outFile;
		outFile.open("error.txt");

		double error;
		double exp;

		for (int i = 10; i <= 70; i++) {
			exp = i/10.0;
			Poisson my_solver;
			my_solver.Initialize(exp);
			my_solver.Integrate_efficient();
			error = my_solver.get_maxerror();

			outFile << setw(15) << setprecision(8) << scientific << (int)pow(10,exp) << " " << error << endl;
		}

		outFile.close();
	}


	
	return 0;
}