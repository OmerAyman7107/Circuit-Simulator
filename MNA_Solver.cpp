#include <iostream>
#include "SOLVER.h"
using namespace std;

// Entry point of the program.
//
// The flow is:
//   1. Create an empty Circuit object.
//   2. Read a netlist file path from the user.
//   3. Call parser() to fill the Circuit with elements and analysis command.
//   4. Call Solver::solve() to run the requested analysis and export results.
//
// The Circuit object owns all the element data; the Solver only reads from
// it and writes to its own result vectors.
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
