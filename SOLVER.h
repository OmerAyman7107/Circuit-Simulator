#include <iostream>
#include <iomanip>
#include <limits>
#include "PARSER.h"
using namespace std;

class Solver
{
public:
	VectorXd dc_results;
	vector<VectorXd> transient_results;
	vector<VectorXcd> ac_results;

	vector<double> freq_points;
	vector<double> time_points;

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

		// if the circuit contains active elements it will begin to
		// perform newton-raphson iterations after getting the first solution
		if (!circuit.diodes.empty())
		{
			int cnt = 100;
			while (cnt--)
			{
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

	void printDC_results(VectorXd& dc_results, Circuit& circuit)
	{
		// printing the results from the vector 

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

		// printing the device currents 
		for (const auto& x : circuit.resistors)
		{
			double V_node_plus = x.node_plus ? dc_results(x.node_plus - 1) : 0;
			double V_node_minus = x.node_minus ? dc_results(x.node_minus - 1) : 0;

			double device_current;
			device_current = (V_node_plus - V_node_minus) / x.value;
			cout << "I(" << x.name << ") = " << device_current << "\tdevice current\n";
		}

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

	void exportAC_results(const vector<VectorXcd>& ac_results, const vector<double>& freq_points, const Circuit& circuit)
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

		//fout << circuit.analysis.type_line << endl;
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

	void solveOP(Circuit& circuit)
	{
		SparseMatrix<double> G_matrix;
		VectorXd RHS_Ivector;
		circuit.buildDC(G_matrix, RHS_Ivector, dc_results);
		solveNR(G_matrix, RHS_Ivector, dc_results, circuit);
		if (circuit.analysis.type == AnalysisType::OP)
			printDC_results(dc_results, circuit);
	}

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

			// if the circuit contains active elements it will begin to
			// perform newton-raphson iterations after getting the first solution
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

