#include <iostream>
#include "SOLVER.h"
using namespace std;

int main()
{
	cout << "MNA Solver: SparseSolver:\n\n";

	Circuit circuit;
	Solver solver;
	string file_name;

	cout << "Enter netlist file path: \n=> ";
	cin >> file_name;

	parser(file_name, circuit);
	solver.solve(circuit);

	return 0;
}
