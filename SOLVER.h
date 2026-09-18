#include <iostream>
#include <iomanip>
#include <limits>
#include "PARSER.h"
using namespace std;

// The Solver class drives the numerical analysis.
//
// It owns the solution vectors (DC, AC, transient) and dispatches to the
// appropriate routine based on the analysis command parsed from the netlist.
// All element data lives in the Circuit object that is passed in; the
// Solver only reads it and stores the results.
//
// Layout of the solution vector (from top to bottom):
//   [0 .. external_nodes-1]                         node voltages
//   [external_nodes .. external+internal-1]         internal node voltages
//   [... + vs.size()]                               currents through voltage sources
//   [... + vcvs.size()]                             currents through VCVS branches
//   [... + ccvs.size()]                             currents through CCVS branches
//   [... + inductors.size()]                        inductor branch currents (OP only)
class Solver
{
public:
	VectorXd dc_results;
	vector<VectorXd> transient_results;
	vector<VectorXcd> ac_results;

	vector<double> freq_points;
	vector<double> time_points;

	// Top-level dispatcher. Chooses the analysis routine based on
	// circuit.analysis.type.
	void solve(Circuit& circuit)
	{
		if (circuit.analysis.type == AnalysisType::OP)
			solveOP(circuit);
		else if (circuit.analysis.type == AnalysisType::AC_DEC ||
			circuit.analysis.type == AnalysisType::AC_LIN ||
			circuit.analysis.type == AnalysisType::AC_OCT)
			solveAC(circuit);
		else if (circuit.analysis.type == AnalysisType::TRAN)
			solveTRAN(circuit);
	}

	// Solves a (possibly nonlinear) DC operating point using Newton-Raphson.
	//
	// The first iteration linearizes the matrix using the current guess
	// for every diode voltage (initial guess = 0 V). Each subsequent
	// iteration recomputes the diode companion models around the last
	// solution and repeats until every diode's junction voltage stops
	// changing between iterations, or until 100 iterations are exhausted.
	//
	// This routine is called directly by solveOP() and is effectively the
	// per-time-step inner loop for solveTRAN().
	void solveNR(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, VectorXd& dc_results, Circuit& circuit)
	{
		SparseLU<SparseMatrix<double>> solver;

		solver.compute(G_matrix);

		if (solver.info() != Eigen::Success)
		{
			cout << "failed to solve\n";
			return;
		}

		dc_results = solver.solve(RHS_Ivector);

		// If the circuit is purely linear, the first solve is already
		// the answer and no NR iterations are needed.
		if (!circuit.diodes.empty())
		{
			int cnt = 100;
			while (cnt--)
			{
				// Re-linearize every diode around the current solution
				// and check convergence.
				int converged = circuit.updateDC_stamp(G_matrix, RHS_Ivector, dc_results);

				solver.factorize(G_matrix);
				if (solver.info() != Eigen::Success)
				{
					cout << "failed to solve circuit during newton-raphson iterations\n";
					return;
				}
				dc_results = solver.solve(RHS_Ivector);

				if (converged == circuit.diodes.size())
					break;
			}
		}
	}

	// Prints the DC operating point to stdout: node voltages, branch
	// currents of the voltage sources, and derived device currents for
	// resistors and diodes.
	void printDC_results(VectorXd& dc_results, Circuit& circuit)
	{
		// Print node voltages and source/inductor branch currents.
		for (auto& x : circuit.variable_names)
		{
			if (x.index < dc_results.size())
			{
				cout << x.label << " = " << dc_results(x.index);
				if (x.label[0] == 'V')
					cout << "\tvoltage\n";
				else if (x.label[0] == 'I')
					cout << "\tdevice current\n";
			}
		}

		// Derived currents for resistors, computed from the node voltages.
		for (const auto& x : circuit.resistors)
		{
			double V_node_plus = x.node_plus ? dc_results(x.node_plus - 1) : 0;
			double V_node_minus = x.node_minus ? dc_results(x.node_minus - 1) : 0;

			double device_current;
			device_current = (V_node_plus - V_node_minus) / x.value;
			cout << "I(" << x.name << ") = " << device_current << "\tdevice current\n";
		}

		// Diode currents from the Shockley equation. Internal diodes of
		// a BJT (name starts with "__") are skipped.
		for (const auto& x : circuit.diodes)
		{
			if (!x.name.compare(0, 2, "__"))
				continue;
			double V_node_plus = x.node_plus ? dc_results(x.node_plus - 1) : 0;
			double V_node_minus = x.node_minus ? dc_results(x.node_minus - 1) : 0;

			double device_current;
			device_current = x.reverse_saturation_current * (exp(x.current_voltage / (x.n * x.thermal_voltage)) - 1);
			cout << "I(" << x.name << ") = " << device_current << "\tdevice current\n";
		}
	}

	// Asks the user which solution variable to export and writes the
	// AC sweep results to ac_results.txt. Each line contains the
	// frequency, magnitude in dB, and phase in degrees.
	void exportAC_results(const vector<VectorXcd>& ac_results, const vector<double>& freq_points, const Circuit& circuit)
	{
		cout << "please choose the variable you want export: \n";
		int option;
		// Only node voltages and source currents are listed. The loop
		// breaks when it hits an inductor current label, because
		// inductors don't appear in the AC matrix.
		for (int i = 0; i < circuit.variable_names.size(); ++i)
		{
			if (circuit.variable_names[i].label[2] == 'L') break;
			cout << i + 1 << ". " << circuit.variable_names[i].label << endl;
		}
		cout << "=> ";
		cin >> option;
		row_label_pair output_variable = circuit.variable_names[option - 1];

		ofstream fout("ac_results.txt");

		if (fout.fail())
		{
			cout << "failed to export results\n";
			return;
		}

		fout << circuit.analysis.type_line << endl;
		fout << "FREQUENCY \t\t" << output_variable.label << endl;
		for (int i = 0; i < freq_points.size(); ++i)
		{
			VectorXcd temp = ac_results[i];
			complex<double> value = temp(output_variable.index);

			double mag_dB = 20 * log10(abs(value));
			double phase_deg = arg(value) * 180 / PI;

			fout << freq_points[i] << '\t';
			fout << fixed << setprecision(numeric_limits<double>::max_digits10)
				<< mag_dB << ' ' << phase_deg << endl;
		}
		fout.close();
	}

	// Same as exportAC_results, but for the transient sweep. Only the
	// selected variable's time series is exported (not the full vector).
	void exportTRAN_results(const vector<VectorXd>& transient_results, const vector<double>& time_points, const Circuit& circuit)
	{
		cout << "please choose the variable you want export: \n";
		int option;
		for (int i = 0; i < circuit.variable_names.size(); ++i)
		{
			if (circuit.variable_names[i].label[2] == 'L') break;
			cout << i + 1 << ". " << circuit.variable_names[i].label << endl;
		}
		cout << "=> ";
		cin >> option;
		row_label_pair output_variable = circuit.variable_names[option - 1];

		ofstream fout("transient_results.txt");

		if (fout.fail())
		{
			cout << "failed to export results\n";
			return;
		}

		fout << "TIME \t\t" << output_variable.label << endl;
		for (int i = 0; i < time_points.size(); ++i)
		{
			VectorXd temp = transient_results[i];

			fout << time_points[i] << '\t';
			fout << fixed << setprecision(numeric_limits<double>::max_digits10)
				<< temp(output_variable.index) << endl;
		}
		fout.close();
	}

	// Builds the DC MNA system, solves it (with NR if needed), and
	// prints the result if the current analysis is OP. When called from
	// solveAC or solveTRAN, the OP solution is used internally (to
	// initialize caps/inductors or linearize diodes) but not printed.
	void solveOP(Circuit& circuit)
	{
		SparseMatrix<double> G_matrix;
		VectorXd RHS_Ivector;
		circuit.buildDC(G_matrix, RHS_Ivector, dc_results);
		solveNR(G_matrix, RHS_Ivector, dc_results, circuit);
		if (circuit.analysis.type == AnalysisType::OP)
			printDC_results(dc_results, circuit);
	}

	// Advances the current frequency by one sweep step.
	// Returns false when freq_stop has already been reached.
	//
	//   AC_LIN: linearly spaced, total of "points" values from fstart to fstop.
	//   AC_DEC: "points" per decade — step factor is 10^(1/points).
	//   AC_OCT: "points" per octave  — step factor is  2^(1/points).
	bool get_next_freq_step(const AnalysisCommand& analysis, double& frequency)
	{
		if (frequency >= analysis.freq_stop)
			return false;

		if (analysis.type == AnalysisType::AC_DEC)
		{
			frequency *= pow(10.0, 1.0 / (analysis.points));
			if (frequency >= analysis.freq_stop)
				frequency = analysis.freq_stop;
			return true;
		}

		else if (analysis.type == AnalysisType::AC_OCT)
		{
			frequency *= pow(2.0, 1.0 / (analysis.points));
			if (frequency >= analysis.freq_stop)
				frequency = analysis.freq_stop;
			return true;
		}

		else if (analysis.type == AnalysisType::AC_LIN)
		{
			double step_size = (analysis.freq_stop - analysis.freq_start) / (analysis.points - 1);
			frequency += step_size;
			return true;
		}
		return false;
	}

	// Runs the AC sweep.
	//
	// Steps:
	//   1. Solve the DC operating point first — diodes use that solution
	//      to pick a linearization point for the small-signal analysis.
	//   2. Build the complex MNA matrix at the starting frequency.
	//   3. Factorize once, then solve at each frequency, updating only
	//      the frequency-dependent stamps (C and L) between steps.
	//
	// The matrix pattern doesn't change between frequencies, so the LU
	// pattern is computed once and reused.
	void solveAC(Circuit& circuit)
	{
		solveOP(circuit);
		SparseMatrix<complex<double>> G_matrix;
		VectorXcd RHS_Ivector;
		double frequency = circuit.analysis.freq_start;
		circuit.buildAC(G_matrix, RHS_Ivector, frequency);

		SparseLU<SparseMatrix<complex<double>>> solver;
		solver.analyzePattern(G_matrix);

		bool valid_freq = true;
		while (valid_freq)
		{
			solver.factorize(G_matrix);
			if (solver.info() != Eigen::Success) { cout << "failed to solve AC\n"; break; }
			ac_results.push_back(solver.solve(RHS_Ivector));

			double old_freq = frequency;
			freq_points.push_back(old_freq);
			valid_freq = get_next_freq_step(circuit.analysis, frequency);
			if (valid_freq)
				circuit.updateAC_stamp(G_matrix, RHS_Ivector, frequency, old_freq);
		}
		cout << "AC analysis completed successfully\n\n";
		exportAC_results(ac_results, freq_points, circuit);
	}

	// Runs the transient analysis.
	//
	// Steps:
	//   1. Solve the DC operating point; its solution is the initial state
	//      of the capacitors and inductors (assumed to be the equilibrium
	//      state at t=0).
	//   2. For each time step:
	//        - On the first step, build the full MNA matrix and factorize.
	//        - On later steps, only update the RHS (sources, cap/ind
	//          companion models) and re-factorize — the matrix pattern
	//          never changes.
	//        - Run Newton-Raphson because the diodes make the system
	//          nonlinear. Cap and inductor companion models use the
	//          previous step's solution, which is only updated once per
	//          time step, not on every NR iteration.
	//   3. Export the requested variable to a text file.
	void solveTRAN(Circuit& circuit)
	{
		solveOP(circuit);

		transient_results.push_back(dc_results);
		time_points.push_back(0.0);

		SparseMatrix<double> G_matrix;
		VectorXd RHS_Ivector;

		SparseLU<SparseMatrix<double>> solver;
		VectorXd temp;
		bool first_step = true;
		for (double time = circuit.analysis.tstep;
			time <= circuit.analysis.tstop;
			time += circuit.analysis.tstep)
		{
			if (first_step)
			{
				circuit.buildTRAN(G_matrix, RHS_Ivector, dc_results);
				solver.compute(G_matrix);
				first_step = false;
			}
			else
			{
				circuit.updateTRAN_stamp(G_matrix, RHS_Ivector, temp, time);
				solver.factorize(G_matrix);
			}

			if (solver.info() != Eigen::Success)
			{
				cout << "failed to solve\n";
				return;
			}

			temp = solver.solve(RHS_Ivector);

			// Newton-Raphson loop: re-linearize diodes around the
			// current iterate until every diode's junction voltage
			// converges.
			if (!circuit.diodes.empty())
			{
				int cnt = 100;
				while (cnt--)
				{
					int converged = circuit.updateDC_stamp(G_matrix, RHS_Ivector, temp);

					solver.factorize(G_matrix);
					if (solver.info() != Eigen::Success)
					{
						cout << "failed to solve circuit during newton-raphson iterations\n";
						return;
					}
					temp = solver.solve(RHS_Ivector);

					if (converged == circuit.diodes.size())
						break;
				}
			}

			transient_results.push_back(temp);
			time_points.push_back(time);
		}
		cout << "Transient analysis completed successfully\n\n";
		exportTRAN_results(transient_results, time_points, circuit);
	}
};
