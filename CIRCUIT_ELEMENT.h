#include <iostream>
#include <string>
#include <complex>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double PI = 3.14159265358979323846;

//////////////////////////////////////////////////
//////////////*** ELEMENT MODELS ***//////////////
//////////////////////////////////////////////////


//////////////*** SOURCES ***//////////////


// voltage source
class Voltage
{
public:
	string name;
	int node_plus;
	int node_minus;
	double dc_magnitude;
	double ac_magnitude;
	double phase;
	bool is_ac = false;
	bool is_tran = false;
	int source_idx; // has to be initialized in the parser

	// the following values are only useful if the source is transient
	double sine_offset = 0.0;
	double sine_amplitude = 0.0;
	double sine_frequency = 0.0;
	double sine_delay = 0.0;
	double sine_theta = 0.0;
	double sine_phase = 0.0;


	Voltage() { is_ac = false; dc_magnitude = ac_magnitude = 0; phase = 0; source_idx = 0; node_plus = node_minus = -1; }

	double source_at(const double& t)
	{
		if (is_tran)
		{
			if (t < sine_delay) return sine_offset + sine_amplitude * sin(sine_phase * PI / 180.0);
			double tau = t - sine_delay;
			double damping = exp(-tau * sine_theta);
			double angle = 2.0 * PI * sine_frequency * tau + sine_phase * PI / 180.0;
			return sine_offset + sine_amplitude * damping * sin(angle);
		}
		return dc_magnitude;  // fallback: constant DC
	}

	void stampDC(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector, const int& external_nodes, const int& internal_nodes)
	{

		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int source_row = external_nodes + internal_nodes + i;

		if (node1 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node1 - 1, 1.0));
			matrix_initializer.push_back(Triplet<double>(node1 - 1, source_row, 1.0));
		}
		if (node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node2 - 1, -1.0));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, source_row, -1.0));
		}
		RHS_Ivector(source_row) += source_at(0);


	}

	void stampAC(vector< Triplet< complex<double> > >& matrix_initializer, VectorXcd& RHS_Ivector, const int& external_nodes, const int& internal_nodes)
	{
		// voltage source value will be a complex number
		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int source_row = external_nodes + internal_nodes + i;

		if (node1 != 0)
		{
			matrix_initializer.push_back(Triplet< complex<double> >(source_row, node1 - 1, 1.0));
			matrix_initializer.push_back(Triplet< complex<double> >(node1 - 1, source_row, 1.0));
		}
		if (node2 != 0)
		{
			matrix_initializer.push_back(Triplet< complex<double> >(source_row, node2 - 1, -1.0));
			matrix_initializer.push_back(Triplet< complex<double> >(node2 - 1, source_row, -1.0));
		}
		double phase_rad = phase * PI / 180.0;
		RHS_Ivector(source_row) += polar(ac_magnitude, phase_rad);
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector, const int& external_nodes, const int& internal_nodes, const double& t)
	{
		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int source_row = external_nodes + internal_nodes + i;

		if (node1 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node1 - 1, 1.0));
			matrix_initializer.push_back(Triplet<double>(node1 - 1, source_row, 1.0));
		}
		if (node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node2 - 1, -1.0));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, source_row, -1.0));
		}
		RHS_Ivector(source_row) += source_at(t);
	}

	void updateTRAN_stamp(VectorXd& RHS_Ivector, const int& external_nodes, const int& internal_nodes, const double& time, const double& dt)
	{
		int& i = source_idx;
		int source_row = external_nodes + internal_nodes + i;
		RHS_Ivector(source_row) += -source_at(time - dt) + source_at(time);
	}
};

// current source
class Current
{
public:
	string name;
	int node_plus;
	int node_minus;
	double dc_magnitude;
	double ac_magnitude;
	double phase;
	bool is_ac = false;
	bool is_tran = false;

	// the following values are only useful if the source is transient
	double sine_offset = 0.0;
	double sine_amplitude = 0.0;
	double sine_frequency = 0.0;
	double sine_delay = 0.0;
	double sine_theta = 0.0;
	double sine_phase = 0.0;


	Current() { is_ac = false; dc_magnitude = ac_magnitude = 0; phase = 0; node_plus = node_minus = -1; }

	double source_at(const double& t)
	{
		if (is_tran)
		{
			if (t < sine_delay) return sine_offset + sine_amplitude * sin(sine_phase * PI / 180.0);
			double tau = t - sine_delay;
			double damping = exp(-tau * sine_theta);
			double angle = 2.0 * PI * sine_frequency * tau + sine_phase * PI / 180.0;
			return sine_offset + sine_amplitude * damping * sin(angle);
		}
		return dc_magnitude;  // fallback: constant DC
	}

	void stampDC(VectorXd& RHS_Ivector)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= (source_at(0));
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += (source_at(0));
	}

	void stampAC(VectorXcd& RHS_Ivector)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		double phase_rad = phase * PI / 180.0;
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= polar(ac_magnitude, phase_rad);
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += polar(ac_magnitude, phase_rad);
	}

	void stampTRAN(VectorXd& RHS_Ivector, const double& t)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		double tran_magnitude = source_at(t);
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= tran_magnitude;
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += tran_magnitude;
	}

	void updateTRAN_stamp(VectorXd& RHS_Ivector, const double& time, const double& dt)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		double tran_magnitude = -source_at(time - dt) + source_at(time);
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= tran_magnitude;
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += tran_magnitude;
	}

};

// voltage controlled voltage source
class VCVS
{
public:
	string name;
	int node_plus;
	int node_minus;
	int control_node_plus;
	int control_node_minus;
	double voltage_gain;
	int source_idx;

	void stampDC(vector<Triplet<double>>& matrix_initializer, const int& vs_size, const int& external_nodes, const int& internal_nodes)
	{
		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int& node_control_1 = control_node_plus;
		int& node_control_2 = control_node_minus;
		int source_row = external_nodes + internal_nodes + vs_size + i;
		if (node1 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node1 - 1, 1));
			matrix_initializer.push_back(Triplet<double>(node1 - 1, source_row, 1));
		}
		if (node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node2 - 1, -1));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, source_row, -1));
		}
		if (node_control_1 != 0)
			matrix_initializer.push_back(Triplet<double>(source_row, node_control_1 - 1, -voltage_gain));
		if (node_control_2 != 0)
			matrix_initializer.push_back(Triplet<double>(source_row, node_control_2 - 1, voltage_gain));

	}

	void stampAC(vector<Triplet< complex<double> >>& matrix_initializer, const int& vs_size, const int& external_nodes, const int& internal_nodes)
	{
		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int& node_control_1 = control_node_plus;
		int& node_control_2 = control_node_minus;
		int source_row = external_nodes + internal_nodes + vs_size + i;
		if (node1 != 0)
		{
			matrix_initializer.push_back(Triplet<complex<double>>(source_row, node1 - 1, 1));
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, source_row, 1));
		}
		if (node2 != 0)
		{
			matrix_initializer.push_back(Triplet<complex<double>>(source_row, node2 - 1, -1));
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, source_row, -1));
		}
		if (node_control_1 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(source_row, node_control_1 - 1, -voltage_gain));
		if (node_control_2 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(source_row, node_control_2 - 1, voltage_gain));

	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, const int& vs_size, const int& external_nodes, const int& internal_nodes)
	{
		stampDC(matrix_initializer, vs_size, external_nodes, internal_nodes);
	}
};

// voltage controlled current source
class VCCS
{
public:
	string name;
	int node_plus;
	int node_minus;
	int control_node_plus;
	int control_node_minus;
	double transconductance;

	void stampDC(vector<Triplet<double>>& matrix_initializer)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		int& node_control_1 = control_node_plus;
		int& node_control_2 = control_node_minus;
		if (node1 != 0)
		{
			if (node_control_1 != 0)
				matrix_initializer.push_back(Triplet<double>(node1 - 1, node_control_1 - 1, transconductance));
			if (node_control_2 != 0)
				matrix_initializer.push_back(Triplet<double>(node1 - 1, node_control_2 - 1, -transconductance));
		}
		if (node2 != 0)
		{
			if (node_control_1 != 0)
				matrix_initializer.push_back(Triplet<double>(node2 - 1, node_control_1 - 1, -transconductance));
			if (node_control_2 != 0)
				matrix_initializer.push_back(Triplet<double>(node2 - 1, node_control_2 - 1, transconductance));
		}
	}

	void stampAC(vector< Triplet<complex<double>> >& matrix_initializer)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		int& node_control_1 = control_node_plus;
		int& node_control_2 = control_node_minus;
		if (node1 != 0)
		{
			if (node_control_1 != 0)
				matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node_control_1 - 1, transconductance));
			if (node_control_2 != 0)
				matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node_control_2 - 1, -transconductance));
		}
		if (node2 != 0)
		{
			if (node_control_1 != 0)
				matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node_control_1 - 1, -transconductance));
			if (node_control_2 != 0)
				matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node_control_2 - 1, transconductance));
		}
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer)
	{
		stampDC(matrix_initializer);
	}

};

// current controlled voltage source
class CCVS
{
public:
	string name;
	int node_plus;
	int node_minus;
	string Vcontrol; // the name of the voltage source that the current passes through
	double transimpedence;
	int source_idx;

	void stampDC(vector<Triplet<double>>& matrix_initializer, const vector<Voltage>& vs
		, const int& vs_size, const int& vcvs_size, const int& external_nodes, const int& internal_nodes)
	{
		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int source_row = external_nodes + internal_nodes + vs_size + vcvs_size + i;
		int control_row = external_nodes + internal_nodes;

		bool valid = false;
		for (; control_row - external_nodes - internal_nodes < vs_size && !valid; ++control_row)
			if (Vcontrol == vs[control_row - external_nodes - internal_nodes].name)
			{
				valid = true;
				break;
			}
		if (valid)
		{
			if (node1 != 0)
			{
				matrix_initializer.push_back(Triplet<double>(source_row, node1 - 1, 1));
				matrix_initializer.push_back(Triplet<double>(node1 - 1, source_row, 1));
			}
			if (node2 != 0)
			{
				matrix_initializer.push_back(Triplet<double>(source_row, node2 - 1, -1));
				matrix_initializer.push_back(Triplet<double>(node2 - 1, source_row, -1));
			}
			matrix_initializer.push_back(Triplet<double>(source_row, control_row, -transimpedence));
		}
		else
		{
			cout << "unable to determine control voltage source\n";
			return;
		}

	}

	void stampAC(vector<Triplet<complex<double>>>& matrix_initializer, const vector<Voltage>& vs
		, const int& vs_size, const int& vcvs_size, const int& external_nodes, const int& internal_nodes)
	{
		int& i = source_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int source_row = external_nodes + internal_nodes + vs_size + vcvs_size + i;
		int control_row = external_nodes + internal_nodes;

		bool valid = false;
		for (; control_row - external_nodes - internal_nodes < vs_size && !valid; ++control_row)
			if (Vcontrol == vs[control_row - external_nodes - internal_nodes].name)
			{
				valid = true;
				break;
			}
		if (valid)
		{
			if (node1 != 0)
			{
				matrix_initializer.push_back(Triplet<complex<double>>(source_row, node1 - 1, 1));
				matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, source_row, 1));
			}
			if (node2 != 0)
			{
				matrix_initializer.push_back(Triplet<complex<double>>(source_row, node2 - 1, -1));
				matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, source_row, -1));
			}
			matrix_initializer.push_back(Triplet<complex<double>>(source_row, control_row, -transimpedence));
		}
		else
		{
			cout << "unable to determine control voltage source\n";
			return;
		}

	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, const vector<Voltage>& vs
		, const int& vs_size, const int& vcvs_size, const int& external_nodes, const int& internal_nodes)
	{
		stampDC(matrix_initializer, vs, vs_size, vcvs_size, external_nodes, internal_nodes);
	}

};

// current controlled current source
class CCCS
{
public:
	string name;
	int node_plus;
	int node_minus;
	string Vcontrol; // the name of the voltage source that the current passes through
	double current_gain;

	void stampDC(vector<Triplet<double>>& matrix_initializer, const vector<Voltage>& vs, const int& external_nodes, const int& internal_nodes)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		int control_row = external_nodes + internal_nodes;

		bool valid = false;
		for (; control_row - external_nodes - internal_nodes < vs.size() && !valid; ++control_row)
			if (Vcontrol == vs[control_row - external_nodes - internal_nodes].name)
			{
				valid = true;
				break;
			}

		if (valid)
		{
			if (node1 != 0)
				matrix_initializer.push_back(Triplet<double>(node1 - 1, control_row, current_gain));
			if (node2 != 0)
				matrix_initializer.push_back(Triplet<double>(node2 - 1, control_row, -current_gain));
		}
		else
		{
			cout << "unable to determine control voltage source\n";
			return;
		}
	}

	void stampAC(vector<Triplet<complex<double>>>& matrix_initializer, const vector<Voltage>& vs, const int& external_nodes, const int& internal_nodes)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		int control_row = external_nodes + internal_nodes;

		bool valid = false;
		for (; control_row - external_nodes - internal_nodes < vs.size() && !valid; ++control_row)
			if (Vcontrol == vs[control_row - external_nodes - internal_nodes].name)
			{
				valid = true;
				break;
			}

		if (valid)
		{
			if (node1 != 0)
				matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, control_row, current_gain));
			if (node2 != 0)
				matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, control_row, -current_gain));
		}
		else
		{
			cout << "unable to determine control voltage source\n";
			return;
		}
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, const vector<Voltage>& vs, const int& external_nodes, const int& internal_nodes)
	{
		stampDC(matrix_initializer, vs, external_nodes, internal_nodes);
	}

};


//////////////*** PASSIVES ***//////////////

// resistors
class Resistor
{
public:
	string name;
	int node_plus;
	int node_minus;
	double value;

	void stampDC(vector<Triplet<double>>& matrix_initializer)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		if (node1 != 0)
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node1 - 1, 1.0 / value));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node2 - 1, 1.0 / value));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node2 - 1, -1.0 / value));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node1 - 1, -1.0 / value));
		}
	}

	void stampAC(vector<Triplet<complex<double>>>& matrix_initializer)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		if (node1 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node1 - 1, 1.0 / value));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node2 - 1, 1.0 / value));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node2 - 1, -1.0 / value));
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node1 - 1, -1.0 / value));
		}
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer)
	{
		stampDC(matrix_initializer);
	}
};

// inductor
class Inductor
{
public:
	string name;
	int node_plus;
	int node_minus;
	double value;
	double initial_condition;
	double previous_voltage;
	double previous_current;
	int element_idx;
	Inductor() { initial_condition = 0; previous_voltage = previous_current = 0; }

	void stampDC(vector<Triplet<double>>& matrix_initializer,
		const int& vs_size, const int& vcvs_size, const int& ccvs_size, const int& external_nodes, const int& internal_nodes)
	{
		int& i = element_idx;
		int& node1 = node_plus;
		int& node2 = node_minus;
		int source_row = external_nodes + internal_nodes + vs_size + vcvs_size + ccvs_size + i;

		if (node1 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node1 - 1, 1.0));
			matrix_initializer.push_back(Triplet<double>(node1 - 1, source_row, 1.0));
		}
		if (node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(source_row, node2 - 1, -1.0));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, source_row, -1.0));
		}
	}

	void stampAC(vector<Triplet<complex<double>>>& matrix_initializer, double& frequency)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double omega = 2 * PI * frequency;
		complex<double> admittance = polar(1 / omega / value, -PI / 2);

		if (node1 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node1 - 1, admittance));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node2 - 1, admittance));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node2 - 1, -admittance));
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node1 - 1, -admittance));
		}
	}

	void updateAC_stamp(SparseMatrix< complex<double> >& G_matrix, VectorXcd& RHS_Ivector, const double& new_freq, const double& old_freq)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double omega = 2 * PI * old_freq;
		complex<double> admittance = polar(1 / omega / value, -PI / 2);

		// destamping the old admittance
		if (node1 != 0)
			G_matrix.coeffRef(node1 - 1, node1 - 1) -= admittance;
		if (node2 != 0)
			G_matrix.coeffRef(node2 - 1, node2 - 1) -= admittance;
		if (node1 != 0 && node2 != 0)
		{
			G_matrix.coeffRef(node1 - 1, node2 - 1) += admittance;
			G_matrix.coeffRef(node2 - 1, node1 - 1) += admittance;
		}

		omega = 2 * PI * new_freq;
		admittance = polar(1 / omega / value, -PI / 2);

		// stamping with the new admittance
		if (node1 != 0)
			G_matrix.coeffRef(node1 - 1, node1 - 1) += admittance;
		if (node2 != 0)
			G_matrix.coeffRef(node2 - 1, node2 - 1) += admittance;
		if (node1 != 0 && node2 != 0)
		{
			G_matrix.coeffRef(node1 - 1, node2 - 1) -= admittance;
			G_matrix.coeffRef(node2 - 1, node1 - 1) -= admittance;
		}
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector, const double& dt)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double G_eq = dt / 2 / value;
		double I_eq = (G_eq * previous_voltage + previous_current);

		// stamping G_eq
		if (node1 != 0)
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node1 - 1, G_eq));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node2 - 1, G_eq));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node2 - 1, -G_eq));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node1 - 1, -G_eq));
		}

		// stamping I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) -= I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) += I_eq;;
		}
	}

	void updateTRAN_stamp(VectorXd& RHS_Ivector, const VectorXd& prev_sol, const double& time, const double& dt)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double G_eq = dt / 2 / value;
		double I_eq = (G_eq * previous_voltage + previous_current);

		double V_new = (node1 ? prev_sol(node1 - 1) : 0) - (node2 ? prev_sol(node2 - 1) : 0);
		double I_new = G_eq * V_new + I_eq;

		previous_voltage = V_new;
		previous_current = I_new;

		// destamping old I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) += I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) -= I_eq;;
		}

		I_eq = (G_eq * V_new + I_new);

		// stamping new I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) -= I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) += I_eq;;
		}

	}

};

// capacitor
class Capacitor
{
public:
	string name;
	int node_plus;
	int node_minus;
	double value;
	double initial_condition;
	double previous_voltage;
	double previous_current;

	Capacitor() { initial_condition = 0; previous_voltage = previous_current = 0; }

	void stampDC()
	{
		// a capacitor is an open circuit in DC, so it has no stamp
	}

	void stampAC(vector<Triplet<complex<double>>>& matrix_initializer, double& frequency)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double omega = 2 * PI * frequency;
		complex<double> admittance = polar(omega * value, PI / 2);

		if (node1 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node1 - 1, admittance));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node2 - 1, admittance));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node2 - 1, -admittance));
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node1 - 1, -admittance));
		}
	}

	void updateAC_stamp(SparseMatrix< complex<double> >& G_matrix, VectorXcd& RHS_Ivector, const double& new_freq, const double& old_freq)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double omega = 2 * PI * old_freq;
		complex<double> admittance = polar(omega * value, PI / 2);

		// destamping the old admittance
		if (node1 != 0)
			G_matrix.coeffRef(node1 - 1, node1 - 1) -= admittance;
		if (node2 != 0)
			G_matrix.coeffRef(node2 - 1, node2 - 1) -= admittance;
		if (node1 != 0 && node2 != 0)
		{
			G_matrix.coeffRef(node1 - 1, node2 - 1) += admittance;
			G_matrix.coeffRef(node2 - 1, node1 - 1) += admittance;
		}

		omega = 2 * PI * new_freq;
		admittance = polar(omega * value, PI / 2);

		// stamping with the new admittance
		if (node1 != 0)
			G_matrix.coeffRef(node1 - 1, node1 - 1) += admittance;
		if (node2 != 0)
			G_matrix.coeffRef(node2 - 1, node2 - 1) += admittance;
		if (node1 != 0 && node2 != 0)
		{
			G_matrix.coeffRef(node1 - 1, node2 - 1) -= admittance;
			G_matrix.coeffRef(node2 - 1, node1 - 1) -= admittance;
		}
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector, const double& dt)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double G_eq = 2 * value / dt;
		double I_eq = -(G_eq * previous_voltage + previous_current);

		// stamping G_eq
		if (node1 != 0)
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node1 - 1, G_eq));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node2 - 1, G_eq));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node2 - 1, -G_eq));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node1 - 1, -G_eq));
		}

		// stamping I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) -= I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) += I_eq;
		}

	}

	void updateTRAN_stamp(VectorXd& RHS_Ivector, const VectorXd& prev_sol, const double& time, const double& dt)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double G_eq = 2 * value / dt;
		double I_eq = -(G_eq * previous_voltage + previous_current);

		double V_new = (node1 ? prev_sol(node1 - 1) : 0) - (node2 ? prev_sol(node2 - 1) : 0);
		double I_new = G_eq * V_new + I_eq;

		previous_voltage = V_new;
		previous_current = I_new;

		// destamping old I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) += I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) -= I_eq;;
		}

		I_eq = -(G_eq * V_new + I_new);

		// stamping new I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) -= I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) += I_eq;;
		}

	}
};

//////////////*** ACTIVES ***//////////////

double pnjlim(double vnew, double vold, double vt, double vcrit)
{
	double vlimited;

	if ((vnew > vcrit) && (abs(vnew - vold) > 2.0 * vt)) {
		if (vold > 0.0) {
			double arg = 1.0 + (vnew - vold) / vt;
			vlimited = (arg > 0.0) ? vold + vt * log(arg) : vcrit;
		}
		else {
			vlimited = vt * log(vnew / vt);
		}
	}
	else {
		vlimited = vnew;
	}

	double VMAX = 40.0 * vt;
	if (vlimited > VMAX) vlimited = VMAX;

	return vlimited;
}

class Diode
{
public:
	string name;
	string model_name;
	int node_plus;
	int node_minus;
	double reverse_saturation_current;
	double thermal_voltage;
	double previous_voltage;
	double current_voltage;
	double n;
	double vt;
	double vcrit;
	Diode()
	{
		reverse_saturation_current = 2.52E-15;
		n = 1;
		thermal_voltage = 25e-3;
		previous_voltage = 0;
		current_voltage = 0;
		vt = n * thermal_voltage;
		vcrit = vt * std::log(vt / (std::sqrt(2.0) * reverse_saturation_current));
	}


	void stampDC(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector)
	{
		// the stamp of the diode is made of two stamps 
		// 1 - current source stamp whose value is determined by the device's parameters
		// 2 - resistor whose value is determined by the device's parameters
		// and during the solving the values of these stamps will be updated after every iteration
		// until the solution converges

		int node1 = node_plus;
		int node2 = node_minus;
		double G_eq = reverse_saturation_current * exp(current_voltage / (n * thermal_voltage)) / (n * thermal_voltage);
		double I_eq = reverse_saturation_current * (exp(current_voltage / (n * thermal_voltage)) - 1) - G_eq * current_voltage;

		// stamping G_eq
		if (node1 != 0)
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node1 - 1, G_eq));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node2 - 1, G_eq));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<double>(node1 - 1, node2 - 1, -G_eq));
			matrix_initializer.push_back(Triplet<double>(node2 - 1, node1 - 1, -G_eq));
		}

		// stamping I_eq
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= I_eq;
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += I_eq;
	}

	void updateDC_stamp(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, VectorXd& results)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		double G_eq = reverse_saturation_current * exp(current_voltage / (n * thermal_voltage)) / (n * thermal_voltage);
		double I_eq = reverse_saturation_current * (exp(current_voltage / (n * thermal_voltage)) - 1) - G_eq * current_voltage;

		// de-stamping old G_eq
		if (node1 != 0)
			G_matrix.coeffRef(node1 - 1, node1 - 1) -= G_eq;
		if (node2 != 0)
			G_matrix.coeffRef(node2 - 1, node2 - 1) -= G_eq;
		if (node1 != 0 && node2 != 0)
		{
			G_matrix.coeffRef(node1 - 1, node2 - 1) += G_eq;
			G_matrix.coeffRef(node2 - 1, node1 - 1) += G_eq;
		}

		// de-stamping old I_eq
		if (node1 != 0)
			RHS_Ivector(node1 - 1) += I_eq;
		if (node2 != 0)
			RHS_Ivector(node2 - 1) -= I_eq;

		// 1. Calculate thermal voltage
		//double vt = n * thermal_voltage;

		// 2. Compute critical voltage (or precompute once on device init: x.vcrit)
		//double vcrit = vt * std::log(vt / (std::sqrt(2.0) * reverse_saturation_current));

		// 3. Get raw voltage from matrix results
		double raw_voltage = (node1 != 0 ? results(node1 - 1) : 0) - (node2 != 0 ? results(node2 - 1) : 0);

		// 4. Update states using pnjlim
		previous_voltage = current_voltage;
		current_voltage = pnjlim(raw_voltage, previous_voltage, vt, vcrit);

		G_eq = reverse_saturation_current * exp(current_voltage / (n * thermal_voltage)) / (n * thermal_voltage);
		I_eq = reverse_saturation_current * (exp(current_voltage / (n * thermal_voltage)) - 1) - G_eq * current_voltage;
		//cout << "G_eq = " << G_eq << "  ,   I_eq = " << I_eq << "  ,  V_d = " << x.current_voltage << endl;
		// re-stamping new G_eq
		if (node1 != 0)
			G_matrix.coeffRef(node1 - 1, node1 - 1) += G_eq;
		if (node2 != 0)
			G_matrix.coeffRef(node2 - 1, node2 - 1) += G_eq;
		if (node1 != 0 && node2 != 0)
		{
			G_matrix.coeffRef(node1 - 1, node2 - 1) -= G_eq;
			G_matrix.coeffRef(node2 - 1, node1 - 1) -= G_eq;
		}

		// re-stamping new I_eq
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= I_eq;
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += I_eq;
	}

	void stampAC(vector<Triplet<complex<double>>>& matrix_initializer)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;

		double G_eq = reverse_saturation_current * exp(current_voltage / (n * thermal_voltage)) / (n * thermal_voltage);

		if (node1 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node1 - 1, G_eq));
		if (node2 != 0)
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node2 - 1, G_eq));
		if (node1 != 0 && node2 != 0)
		{
			matrix_initializer.push_back(Triplet<complex<double>>(node1 - 1, node2 - 1, -G_eq));
			matrix_initializer.push_back(Triplet<complex<double>>(node2 - 1, node1 - 1, -G_eq));
		}
	}

	void stampTRAN(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector)
	{
		stampDC(matrix_initializer, RHS_Ivector);
	}

};

class BJT
{
public:
	string name;
	string model_name;
	int collector_node;
	int base_node;
	int emitter_node;
	double alphaF;
	double alphaR;
	int element_idx;
	BJT()
	{
		alphaF = 0.971;
		alphaR = 0.748;
	}

	void expand_device_model(vector<Voltage>& vs, vector<Diode>& diodes, vector<CCCS>& cccs, const int& external_nodes, const int& internal_nodes)
	{
		Voltage Ic_monitor, Ie_monitor, Ib_monitor;
		Voltage IED_monitor, ICD_monitor;
		CCCS forward_CS, reverse_CS;
		Diode BE_junction, BC_junction;

		int& i = element_idx;
		int B1, B2, B3, C1, E1;
		B1 = external_nodes + 5 * i + 1;
		B2 = external_nodes + 5 * i + 2;
		B3 = external_nodes + 5 * i + 3;

		C1 = external_nodes + 5 * i + 4;
		E1 = external_nodes + 5 * i + 5;

		// external monitoring sources
		Ic_monitor.name = name + "_collector_current";
		Ib_monitor.name = name + "_base_current";
		Ie_monitor.name = name + "_emitter_current";

		Ic_monitor.node_plus = collector_node;
		Ic_monitor.node_minus = C1;

		Ib_monitor.node_plus = base_node;
		Ib_monitor.node_minus = B1;

		Ie_monitor.node_plus = E1;
		Ie_monitor.node_minus = emitter_node;

		Ic_monitor.dc_magnitude = 0;
		Ib_monitor.dc_magnitude = 0;
		Ie_monitor.dc_magnitude = 0;

		vs.push_back(Ic_monitor);
		vs.push_back(Ib_monitor);
		vs.push_back(Ie_monitor);

		// internal monitoring sources
		IED_monitor.name = "__" + name + "_internal_IED_monitor";
		ICD_monitor.name = "__" + name + "_internal_ICD_monitor";

		IED_monitor.node_plus = B1;
		IED_monitor.node_minus = B3;

		ICD_monitor.node_plus = B1;
		ICD_monitor.node_minus = B2;

		IED_monitor.dc_magnitude = 0;
		ICD_monitor.dc_magnitude = 0;

		vs.push_back(IED_monitor);
		vs.push_back(ICD_monitor);

		// CCCSs
		forward_CS.name = "__" + name + "_forward_active_collector_current";
		reverse_CS.name = "__" + name + "_reverse_active_collector_current";

		forward_CS.node_plus = C1;
		reverse_CS.node_plus = E1;

		forward_CS.node_minus = B1;
		reverse_CS.node_minus = B1;

		forward_CS.current_gain = alphaF;
		reverse_CS.current_gain = alphaR;

		forward_CS.Vcontrol = IED_monitor.name;
		reverse_CS.Vcontrol = ICD_monitor.name;

		cccs.push_back(forward_CS);
		cccs.push_back(reverse_CS);


		// Junction diodes
		BE_junction.name = "__" + name + "_BE_junction";
		BC_junction.name = "__" + name + "_BC_junction";

		BE_junction.node_plus = B3;
		BE_junction.node_minus = E1;

		BC_junction.node_plus = B2;
		BC_junction.node_minus = C1;

		diodes.push_back(BE_junction);
		diodes.push_back(BC_junction);
	}

};


enum class AnalysisType { NONE, OP, AC_LIN, AC_DEC, AC_OCT, TRAN };

struct AnalysisCommand
{
	AnalysisType type = AnalysisType::NONE;
	string type_line;

	// only meaningful when type is AC_LIN / AC_DEC / AC_OCT:
	int    points = 0;
	double freq_start = 0.0;
	double freq_stop = 0.0;

	// only meaningful when type is TRAN:
	double tstep = 0.0;
	double tstop = 0.0;
};

struct row_label_pair
{
	int index;
	string label;
	row_label_pair(const int& _index, const string& _label)
	{
		index = _index;
		label = _label;
	}
};

class Circuit
{
public:
	// sources:
	vector<Voltage> vs;//
	vector<Current> cs;
	vector<VCVS> vcvs;//
	vector<VCCS> vccs;
	vector<CCVS> ccvs;//
	vector<CCCS> cccs;

	//passives:
	vector<Resistor> resistors;
	vector<Capacitor> caps;
	vector<Inductor> inductors;

	//actives
	vector<Diode> diodes;
	vector<BJT> bipolar_junction_transistors;

	AnalysisCommand analysis;

	int external_nodes = 0;
	int internal_nodes = 0;

	vector<row_label_pair> variable_names;

	void resolve_multi_element_devices()
	{
		for (int i = 0; i < bipolar_junction_transistors.size(); ++i)
			bipolar_junction_transistors[i].element_idx = i;

		for (auto& x : bipolar_junction_transistors)
			x.expand_device_model(vs, diodes, cccs, external_nodes, internal_nodes);
	}

	void set_indices()
	{
		for (int i = 0; i < vs.size(); ++i)
			vs[i].source_idx = i;

		for (int i = 0; i < vcvs.size(); ++i)
			vcvs[i].source_idx = i;

		for (int i = 0; i < ccvs.size(); ++i)
			ccvs[i].source_idx = i;

		for (int i = 0; i < inductors.size(); ++i)
			inductors[i].element_idx = i;
	}

	void buildDC(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, VectorXd& dc_results)
	{
		int matrix_size = external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size() + (int)inductors.size();

		G_matrix.resize(matrix_size, matrix_size);
		RHS_Ivector = VectorXd::Zero(matrix_size);
		dc_results.resize(matrix_size);

		vector<Triplet<double>> matrix_initializer;
		matrix_initializer.reserve(5 * matrix_size);

		// DC stamping sources
		for (auto& x : vs)
			x.stampDC(matrix_initializer, RHS_Ivector, external_nodes, internal_nodes);
		for (auto& x : cs)
			x.stampDC(RHS_Ivector);
		for (auto& x : vcvs)
			x.stampDC(matrix_initializer, (int)vs.size(), external_nodes, internal_nodes);
		for (auto& x : vccs)
			x.stampDC(matrix_initializer);
		for (auto& x : ccvs)
			x.stampDC(matrix_initializer, vs, (int)vs.size(), (int)vcvs.size(), external_nodes, internal_nodes);
		for (auto& x : cccs)
			x.stampDC(matrix_initializer, vs, external_nodes, internal_nodes);

		// DC stamping passives
		for (auto& x : resistors)
			x.stampDC(matrix_initializer);
		for (auto& x : caps)
			x.stampDC();
		for (auto& x : inductors)
			x.stampDC(matrix_initializer, (int)vs.size(), (int)vcvs.size(), (int)ccvs.size(), external_nodes, internal_nodes);

		// DC stamping actives
		for (auto& x : diodes)
			x.stampDC(matrix_initializer, RHS_Ivector);

		G_matrix.setFromTriplets(matrix_initializer.begin(), matrix_initializer.end());
	}

	void buildAC(SparseMatrix<complex<double>>& G_matrix, VectorXcd& RHS_Ivector, double& fstart)
	{
		int matrix_size = external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size();

		G_matrix.resize(matrix_size, matrix_size);
		RHS_Ivector = VectorXcd::Zero(matrix_size);


		vector<Triplet <complex< double >> > matrix_initializer;
		matrix_initializer.reserve(5 * matrix_size);

		// DC stamping sources
		for (auto& x : vs)
			x.stampAC(matrix_initializer, RHS_Ivector, external_nodes, internal_nodes);
		for (auto& x : cs)
			x.stampAC(RHS_Ivector);
		for (auto& x : vcvs)
			x.stampAC(matrix_initializer, (int)vs.size(), external_nodes, internal_nodes);
		for (auto& x : vccs)
			x.stampAC(matrix_initializer);
		for (auto& x : ccvs)
			x.stampAC(matrix_initializer, vs, (int)vs.size(), (int)vcvs.size(), external_nodes, internal_nodes);
		for (auto& x : cccs)
			x.stampAC(matrix_initializer, vs, external_nodes, internal_nodes);

		// DC stamping passives
		for (auto& x : resistors)
			x.stampAC(matrix_initializer);
		for (auto& x : caps)
			x.stampAC(matrix_initializer, fstart);
		for (auto& x : inductors)
			x.stampAC(matrix_initializer, fstart);

		// DC stamping actives
		for (auto& x : diodes)
			x.stampAC(matrix_initializer);

		G_matrix.setFromTriplets(matrix_initializer.begin(), matrix_initializer.end());
	}

	void buildTRAN(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, VectorXd& dc_results)
	{
		int matrix_size = external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size();

		G_matrix.resize(matrix_size, matrix_size);
		RHS_Ivector = VectorXd::Zero(matrix_size);

		vector<Triplet<double>> matrix_initializer;
		matrix_initializer.reserve(5 * matrix_size);

		for (auto& x : caps)
		{
			double V_plus = x.node_plus ? dc_results(x.node_plus - 1) : 0.0;
			double V_minus = x.node_minus ? dc_results(x.node_minus - 1) : 0.0;
			x.previous_voltage = V_plus - V_minus;
			x.previous_current = 0.0;
		}

		for (auto& x : inductors)
		{
			x.previous_voltage = 0.0;
			x.previous_current = dc_results(external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size() + x.element_idx);
		}

		// DC stamping sources
		for (auto& x : vs)
			x.stampTRAN(matrix_initializer, RHS_Ivector, external_nodes, internal_nodes, analysis.tstep);
		for (auto& x : cs)
			x.stampTRAN(RHS_Ivector, analysis.tstep);
		for (auto& x : vcvs)
			x.stampTRAN(matrix_initializer, (int)vs.size(), external_nodes, internal_nodes);
		for (auto& x : vccs)
			x.stampTRAN(matrix_initializer);
		for (auto& x : ccvs)
			x.stampTRAN(matrix_initializer, vs, (int)vs.size(), (int)vcvs.size(), external_nodes, internal_nodes);
		for (auto& x : cccs)
			x.stampTRAN(matrix_initializer, vs, external_nodes, internal_nodes);

		// DC stamping passives
		for (auto& x : resistors)
			x.stampTRAN(matrix_initializer);
		for (auto& x : caps)
			x.stampTRAN(matrix_initializer, RHS_Ivector, analysis.tstep);
		for (auto& x : inductors)
			x.stampTRAN(matrix_initializer, RHS_Ivector, analysis.tstep);

		// DC stamping actives
		for (auto& x : diodes)
			x.stampTRAN(matrix_initializer, RHS_Ivector);

		G_matrix.setFromTriplets(matrix_initializer.begin(), matrix_initializer.end());
	}

	// DC stamps are only updated during NR iterations
	int updateDC_stamp(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, VectorXd& dc_results)
	{
		int converged = 0;
		for (auto& x : diodes)
		{
			x.updateDC_stamp(G_matrix, RHS_Ivector, dc_results);
			double tol = 1e-6 * max(abs(x.current_voltage), abs(x.previous_voltage)) + 1e-9;
			if (abs(x.current_voltage - x.previous_voltage) <= tol)
				converged++;
		}
		return converged;
	}

	void updateAC_stamp(SparseMatrix< complex<double> >& G_matrix, VectorXcd& RHS_Ivector, const double& new_freq, const double& old_freq)
	{
		for (auto& x : caps)
			x.updateAC_stamp(G_matrix, RHS_Ivector, new_freq, old_freq);
		for (auto& x : inductors)
			x.updateAC_stamp(G_matrix, RHS_Ivector, new_freq, old_freq);
	}

	void updateTRAN_stamp(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, const VectorXd& prev_sol, const double& time)
	{
		// update independent sources' stamps
		for (auto& x : vs)
			x.updateTRAN_stamp(RHS_Ivector, external_nodes, internal_nodes, time, analysis.tstep);
		for (auto& x : cs)
			x.updateTRAN_stamp(RHS_Ivector, time, analysis.tstep);

		// update caps and inductor stamps
		for (auto& x : caps)
			x.updateTRAN_stamp(RHS_Ivector, prev_sol, time, analysis.tstep);
		for (auto& x : inductors)
			x.updateTRAN_stamp(RHS_Ivector, prev_sol, time, analysis.tstep);

	}


	void set_variable_names()
	{
		// printing the results from the vector 
		int i = 0;
		string temp;
		for (; i < external_nodes; ++i)
		{
			temp = "V(" + to_string(1 + i) + ')';
			variable_names.push_back(row_label_pair(i, temp));
		}

		i += internal_nodes;

		for (; i - external_nodes - internal_nodes < vs.size(); ++i)
		{
			if (vs[i - external_nodes - internal_nodes].name.compare(0, 2, "__"))
			{
				temp = "I(" + vs[i - external_nodes - internal_nodes].name + ")";
				variable_names.push_back(row_label_pair(i, temp));
			}
		}

		for (; i - external_nodes - internal_nodes - vs.size() < vcvs.size(); ++i)
		{
			temp = "I(" + vcvs[i - external_nodes - internal_nodes - vs.size()].name + ")";
			variable_names.push_back(row_label_pair(i, temp));
		}
		for (; i - external_nodes - internal_nodes - vs.size() - vcvs.size() < ccvs.size(); ++i)
		{
			temp = "I(" + ccvs[i - external_nodes - internal_nodes - vs.size() - vcvs.size()].name + ")";
			variable_names.push_back(row_label_pair(i, temp));
		}
		for (; i - external_nodes - internal_nodes - vs.size() - vcvs.size() - ccvs.size() < inductors.size(); ++i)
		{
			temp = "I(" + inductors[i - external_nodes - internal_nodes - vs.size() - vcvs.size() - ccvs.size()].name + ")";
			variable_names.push_back(row_label_pair(i, temp));
		}
	}
};
