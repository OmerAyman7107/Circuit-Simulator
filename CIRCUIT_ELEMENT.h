#include <iostream>
#include <string>
#include <complex>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double PI = 3.14159265358979323846;

// This file contains all the element classes and their helper methods
// that are used to build and update the MNA matrix according to the mode
// of analysis.
//
// Layout convention:
//   Each element class provides up to four stamp/update methods:
//     - stampDC()        : contribution to the real DC matrix and RHS
//     - stampAC()        : contribution to the complex AC matrix and RHS
//     - stampTRAN()      : contribution to the real transient matrix and RHS
//     - update*_stamp()  : in-place update of an already-built matrix,
//                          used to avoid rebuilding from scratch
//
//   Node 0 is ground and is never assigned a row or column in the matrix.
//   The expression "node - 1" therefore maps a SPICE node number to its
//   row/column index in the MNA matrix.

//////////////////////////////////////////////////
//////////////*** ELEMENT MODELS ***//////////////
//////////////////////////////////////////////////


//////////////*** SOURCES ***//////////////



// The voltage source class stores the information about itself after
// parsing the netlist: name, node connections, DC and AC magnitudes, and
// optional SINE parameters.
//
// A voltage source requires an extra row/column in the MNA matrix: one
// for the KVL equation that enforces V(n+) - V(n-) = value, and one for
// the branch current that is added as an unknown to the KCL equations at
// the source's nodes.
class Voltage
{
public:
	string name;
	int node_plus;
	int node_minus;

	// DC and AC magnitudes are stored separately because a source can
	// carry both (e.g. "V1 1 0 DC 5 AC 1"). Both default to zero.
	double dc_magnitude;
	double ac_magnitude;

	// Phase of the AC small-signal voltage, in degrees. Defaults to 0.
	double phase;

	// Flags describing which analyses this source participates in.
	// is_ac marks the source as an AC excitation; is_tran marks it as
	// a transient (SIN/SINE) source. Both default to false.
	bool is_ac = false;
	bool is_tran = false;

	// Position of this source's extra matrix row.
	// Every voltage source, VCVS branch, and CCVS branch takes one row
	// after the node rows. Assigned in Circuit::set_indices().
	int source_idx; // has to be initialized in the parser

	// SINE waveform parameters. Only meaningful when is_tran is true.
	double sine_offset = 0.0;
	double sine_amplitude = 0.0;
	double sine_frequency = 0.0;
	double sine_delay = 0.0;
	double sine_theta = 0.0;
	double sine_phase = 0.0;


	Voltage() { is_ac = false; dc_magnitude = ac_magnitude = 0; phase = 0; source_idx = 0; node_plus = node_minus = -1; }

	// Returns the instantaneous value of the source at time t.
	//
	// For a transient source the waveform is
	//     V(t) = offset + amplitude * exp(-theta * tau) * sin(2*pi*f*tau + phase)
	// where tau = t - delay. If t < delay the source is held at its
	// initial value (offset + amplitude * sin(phase)).
	//
	// For a non-transient source the DC magnitude is returned; this makes
	// the same function usable as the OP value.
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

	// Stamps the voltage source into the DC MNA system.
	//
	// Two matrix entries are needed per connected node:
	//   Row source_row: V(node1) - V(node2) = value
	//   Column source_row at each node: current injected by the branch
	// The RHS gets the source value at t = 0.
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

	// Same structure as stampDC, but with complex entries and a phasor
	// RHS value: ac_magnitude at the given phase angle (converted from
	// degrees to radians).
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

	// Stamps the voltage source for transient analysis, using the
	// instantaneous source value at time t. The matrix entries are
	// identical to those of stampDC.
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

	// Updates the RHS after a time step by adding the difference between
	// the new and old source values. The matrix entries are unchanged, so
	// only the RHS row needs touching. This avoids re-stamping the whole
	// source on every transient step.
	void updateTRAN_stamp(VectorXd& RHS_Ivector, const int& external_nodes, const int& internal_nodes, const double& time, const double& dt)
	{
		int& i = source_idx;
		int source_row = external_nodes + internal_nodes + i;
		RHS_Ivector(source_row) += -source_at(time - dt) + source_at(time);
	}
};



// The current source class mirrors the voltage source class: same fields,
// same waveform options, same structure. The only difference is the stamp:
// a current source has no branch current unknown, so it contributes only
// to the RHS (current injection at its two nodes).
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

	// SINE parameters, only meaningful when is_tran is true.
	double sine_offset = 0.0;
	double sine_amplitude = 0.0;
	double sine_frequency = 0.0;
	double sine_delay = 0.0;
	double sine_theta = 0.0;
	double sine_phase = 0.0;


	Current() { is_ac = false; dc_magnitude = ac_magnitude = 0; phase = 0; node_plus = node_minus = -1; }

	// Same waveform evaluation as Voltage::source_at.
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

	// Stamps the DC current: +I into node_plus, -I out of node_minus.
	// The sign convention matches SPICE (current flows from n+ to n-
	// through the source internally, so it is injected into n- and
	// extracted from n+).
	void stampDC(VectorXd& RHS_Ivector)
	{
		int& node1 = node_plus;
		int& node2 = node_minus;
		if (node1 != 0)
			RHS_Ivector(node1 - 1) -= (source_at(0));
		if (node2 != 0)
			RHS_Ivector(node2 - 1) += (source_at(0));
	}

	// Complex phasor injection for AC analysis.
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

	// Instantaneous injection at time t for transient analysis.
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

	// Incremental RHS update between time steps.
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



// Voltage-controlled voltage source (VCVS).
// Models an ideal amplifier: V(n+) - V(n-) = gain * (V(nc+) - V(nc-)).
// Like a plain voltage source, it contributes an extra matrix row/column
// for its branch current.
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

	// Stamps the VCVS into the DC system.
	//
	// Row source_row:    V(n+) - V(n-) - gain*V(nc+) + gain*V(nc-) = 0
	// Column source_row: branch current injected into n+ and drawn from n-
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

	// Complex version of stampDC. Identical structure, complex entries.
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

	// The VCVS is linear and time-invariant, so its transient stamp is
	// identical to its DC stamp.
	void stampTRAN(vector<Triplet<double>>& matrix_initializer, const int& vs_size, const int& external_nodes, const int& internal_nodes)
	{
		stampDC(matrix_initializer, vs_size, external_nodes, internal_nodes);
	}
};



// Voltage-controlled current source (VCCS).
// Output current: I = transconductance * (V(nc+) - V(nc-)), flowing from
// node_minus to node_plus internally. No extra matrix row is needed —
// the branch current is a pure function of node voltages.
class VCCS
{
public:
	string name;
	int node_plus;
	int node_minus;
	int control_node_plus;
	int control_node_minus;
	double transconductance;

	// Stamps the VCCS conductances:
	//   Row n+:  +gm * V(nc+) - gm * V(nc-)
	//   Row n-:  -gm * V(nc+) + gm * V(nc-)
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

	// Complex version of stampDC.
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

	// Same as stampDC.
	void stampTRAN(vector<Triplet<double>>& matrix_initializer)
	{
		stampDC(matrix_initializer);
	}
};



// Current-controlled voltage source (CCVS).
// Output voltage: V(n+) - V(n-) = transimpedance * I_control, where
// I_control is the current through the named voltage source.
//
// The control current is the branch-current unknown of that voltage
// source, so this element requires an extra matrix row/column of its
// own (for its branch current).
class CCVS
{
public:
	string name;
	int node_plus;
	int node_minus;
	string Vcontrol; // name of the voltage source whose current controls this CCVS
	double transimpedence;
	int source_idx;

	// Stamps the CCVS into the DC system.
	//
	// Row source_row:    V(n+) - V(n-) - Zt * I_control = 0
	// Column source_row: branch current injected into n+ and drawn from n-
	//
	// The control row is located by scanning the vs vector for a matching
	// name. If no match is found the element is skipped and an error is
	// printed.
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

	// Complex version of stampDC.
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

	// Same as stampDC.
	void stampTRAN(vector<Triplet<double>>& matrix_initializer, const vector<Voltage>& vs
		, const int& vs_size, const int& vcvs_size, const int& external_nodes, const int& internal_nodes)
	{
		stampDC(matrix_initializer, vs, vs_size, vcvs_size, external_nodes, internal_nodes);
	}
};



// Current-controlled current source (CCCS).
// Output current: I = current_gain * I_control, where I_control is the
// current through the named voltage source. No extra matrix row needed.
class CCCS
{
public:
	string name;
	int node_plus;
	int node_minus;
	string Vcontrol; // name of the voltage source whose current controls this CCCS
	double current_gain;

	// Stamps the CCCS contribution to the KCL rows at n+ and n-.
	// The current is expressed in terms of the control source's branch
	// current unknown.
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

	// Complex version of stampDC.
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

	// Same as stampDC.
	void stampTRAN(vector<Triplet<double>>& matrix_initializer, const vector<Voltage>& vs, const int& external_nodes, const int& internal_nodes)
	{
		stampDC(matrix_initializer, vs, external_nodes, internal_nodes);
	}
};


//////////////*** PASSIVES ***//////////////



// Resistor: linear, time-invariant, purely real.
// Contributes the standard 1/R conductance stamp to all four quadrants
// of its two-node submatrix.
class Resistor
{
public:
	string name;
	int node_plus;
	int node_minus;
	double value;

	// DC stamp: standard 2x2 conductance matrix.
	//   G(n+,n+) = +1/R   G(n-,n-) = +1/R
	//   G(n+,n-) = -1/R   G(n-,n+) = -1/R
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

	// Complex version of stampDC. A resistor is real, so the entries
	// are the same numbers but stored as complex.
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

	// Same as stampDC.
	void stampTRAN(vector<Triplet<double>>& matrix_initializer)
	{
		stampDC(matrix_initializer);
	}
};



// Inductor.
//
//   DC    : short circuit. The branch current becomes an extra unknown,
//           so an inductor contributes a row/column like a 0 V source.
//   AC    : complex admittance 1 / (j*omega*L).
//   TRAN  : trapezoidal companion model — conductance G_eq in parallel
//           with a current source I_eq. Both are derived from the
//           previous time step's voltage and current.
class Inductor
{
public:
	string name;
	int node_plus;
	int node_minus;
	double value;

	// Reserved for a future initial-condition feature. Currently unused
	// (the OP solution initializes inductors instead).
	double initial_condition;

	// State from the previous time step, needed to update the transient
	// companion model.
	double previous_voltage;
	double previous_current;

	// Position of this inductor's extra matrix row (used only in DC).
	// Assigned in Circuit::set_indices().
	int element_idx;
	Inductor() { initial_condition = 0; previous_voltage = previous_current = 0; }

	// DC stamp: behaves as a short circuit, contributing an extra row
	// enforcing V(n+) - V(n-) = 0, and an extra current unknown.
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

	// AC stamp: Y = 1 / (j*omega*L).
	// Computed via polar form so that the phase (-90 degrees) is explicit.
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

	// Removes the old frequency's admittance and adds the new one in place.
	// Called between AC sweep steps to avoid rebuilding the matrix.
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

	// Transient stamp: trapezoidal companion model.
	//   G_eq = dt / (2L)
	//   I_eq = G_eq * v_prev + i_prev
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
				RHS_Ivector(node2 - 1) += I_eq;
		}
	}

	// Recomputes the companion model from the previous step's solution
	// and updates the RHS accordingly. The matrix (G_eq) is unchanged
	// between steps since dt and L are constant, so only the RHS moves.
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
				RHS_Ivector(node2 - 1) -= I_eq;
		}

		I_eq = (G_eq * V_new + I_new);

		// stamping new I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) -= I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) += I_eq;
		}
	}
};



// Capacitor.
//
//   DC    : open circuit — no stamp at all.
//   AC    : complex admittance j*omega*C.
//   TRAN  : trapezoidal companion model — conductance G_eq in parallel
//           with a current source I_eq. Both are derived from the
//           previous time step's voltage and current.
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

	// No stamp for DC. A capacitor is an open circuit, contributing
	// nothing to the conductance matrix or the RHS.
	void stampDC()
	{
		// a capacitor is an open circuit in DC, so it has no stamp
	}

	// AC stamp: Y = j*omega*C.
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

	// In-place update of the AC stamp between frequency steps.
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

	// Transient stamp: trapezoidal companion model.
	//   G_eq = 2C / dt
	//   I_eq = -(G_eq * v_prev + i_prev)
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

	// Recomputes the companion model from the previous solution and
	// updates the RHS. The matrix (G_eq) does not change between steps.
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
				RHS_Ivector(node2 - 1) -= I_eq;
		}

		I_eq = -(G_eq * V_new + I_new);

		// stamping new I_eq
		{
			if (node1 != 0)
				RHS_Ivector(node1 - 1) -= I_eq;
			if (node2 != 0)
				RHS_Ivector(node2 - 1) += I_eq;
		}
	}
};

//////////////*** ACTIVES ***//////////////


// pnjlim — the classic SPICE PN-junction voltage limiter.
//
// Prevents the Newton-Raphson iteration from taking huge jumps in the
// diode voltage, which would make exp(Vd/Vt) explode and destroy the
// Jacobian. A hard ceiling (VMAX = 40 * Vt) is applied on top of the
// log-damping step-reduction.
//
//   vnew  : raw solver voltage for this iteration
//   vold  : voltage from the previous iteration
//   vt    : thermal voltage n * kT/q
//   vcrit : critical voltage above which limiting is applied
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


// Diode — the basic Shockley model.
//
//   I(Vd) = Is * (exp(Vd / (n * Vt)) - 1)
//
// Limitations (by design, to keep the model simple):
//   - no series resistance
//   - no junction capacitance
//   - no reverse breakdown
//   - no temperature dependence
//
// The linearized companion model (conductance G_eq in parallel with a
// current source I_eq) is derived by first-order Taylor expansion around
// the current estimate of the diode voltage:
//
//   G_eq = Is * exp(Vd / (n*Vt)) / (n*Vt)
//   I_eq = Is * (exp(Vd / (n*Vt)) - 1) - G_eq * Vd
//
// During Newton-Raphson, updateDC_stamp() recomputes these values from
// the latest solution and re-stamps the matrix in place.
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
		// Default parameters match a typical silicon small-signal diode.
		reverse_saturation_current = 2.52E-15;
		n = 1;
		thermal_voltage = 25e-3;
		previous_voltage = 0;
		current_voltage = 0;
		vt = n * thermal_voltage;
		vcrit = vt * std::log(vt / (std::sqrt(2.0) * reverse_saturation_current));
	}

	// Stamps the initial diode linearization at Vd = current_voltage
	// (usually 0 on the first pass). Subsequent NR iterations update
	// these values via updateDC_stamp().
	void stampDC(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector)
	{
		// The stamp is the small-signal equivalent of the linearized
		// diode: a conductance G_eq in parallel with a current source
		// I_eq. G_eq is placed in the four quadrants of the 2x2 node
		// submatrix; I_eq is injected into the RHS at the two nodes.

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

	// Advances the diode linearization by one Newton-Raphson step.
	//
	// Pattern:
	//   1. De-stamp the old G_eq / I_eq.
	//   2. Read the raw junction voltage from the solver result and pass
	//      it through pnjlim to get the next estimate.
	//   3. Re-stamp the new G_eq / I_eq.
	//
	// Modifying the matrix in place (instead of rebuilding from triplets)
	// is much cheaper per NR iteration. This is safe because the pattern
	// does not change.
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

		// Read the diode voltage out of the last solution.
		double raw_voltage = (node1 != 0 ? results(node1 - 1) : 0) - (node2 != 0 ? results(node2 - 1) : 0);

		// Limit the step to keep the exponential from blowing up.
		previous_voltage = current_voltage;
		current_voltage = pnjlim(raw_voltage, previous_voltage, vt, vcrit);

		G_eq = reverse_saturation_current * exp(current_voltage / (n * thermal_voltage)) / (n * thermal_voltage);
		I_eq = reverse_saturation_current * (exp(current_voltage / (n * thermal_voltage)) - 1) - G_eq * current_voltage;

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

	// Small-signal stamp for AC analysis: the diode is replaced by its
	// linearized conductance G_eq, evaluated at the DC operating point.
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

	// Transient initial stamp is the same as the DC stamp; the NR loop
	// updates it identically in every time step.
	void stampTRAN(vector<Triplet<double>>& matrix_initializer, VectorXd& RHS_Ivector)
	{
		stampDC(matrix_initializer, RHS_Ivector);
	}
};



// BJT — expanded Ebers-Moll model.
//
// Rather than stamping a compact nonlinear model directly, the BJT is
// expanded into a subcircuit of simpler primitives (diodes, CCCSs, and
// ideal voltage sources acting as ammeters) that are then handled by the
// existing element classes.
//
// For each transistor, the following are added to the Circuit:
//   - 3 external voltage sources (0 V) that measure Ic, Ib, Ie
//   - 2 internal voltage sources (0 V) that measure the currents through
//     the two junction diodes
//   - 2 CCCSs that model the forward- and reverse-active collector current
//   - 2 junction diodes for the base-emitter and base-collector junctions
//
// 5 internal nodes are allocated per BJT and wired as:
//   B1 = base side of the external Ib monitor
//   B2 = base side of the BC junction
//   B3 = base side of the BE junction
//   C1 = collector side of the Ic monitor and BC junction
//   E1 = emitter side of the Ie monitor and BE junction
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
		// Forward and reverse common-base current gains.
		// Typical silicon values for the Ebers-Moll model.
		alphaF = 0.971;
		alphaR = 0.748;
	}

	// Expands this BJT into its primitive subcircuit and appends the
	// new elements to the vs, diodes, and cccs vectors. Must be called
	// after parsing and before stamping.
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

		// External monitoring sources (0 V ammeters).
		// These measure the terminal currents so they appear in the
		// output as I(Qname_collector_current) etc.
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

		// Internal monitoring sources. These provide the control currents
		// for the two CCCSs, so they are prefixed with "__" to signal
		// that they should not appear in the variable list.
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

		// Current-controlled current sources: alpha_F * I_BE and
		// alpha_R * I_BC, injected between C1/B1 and E1/B1 respectively.
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

		// Junction diodes for the two PN junctions of the BJT.
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

// Stores the analysis command parsed from a dot-command (.OP, .AC, .TRAN).
// Which fields are meaningful depends on the type:
//   OP:    no additional parameters
//   AC_*:  points, freq_start, freq_stop
//   TRAN:  tstep, tstop
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

// Pairs a solution-vector row index with a human-readable label such as
// "V(3)" or "I(V1)". Used by the result-printing and export functions.
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

// The Circuit class holds everything the solver needs to know about the
// netlist: the parsed elements, node counts, analysis command, and the
// helper methods that build and update the MNA matrix.
//
// MNA matrix layout:
//   Rows [0 .. external_nodes-1]                      node voltages
//   Rows [external_nodes .. +internal_nodes-1]        internal node voltages
//   Rows [.. + vs.size()]                             vs branch currents
//   Rows [.. + vcvs.size()]                           vcvs branch currents
//   Rows [.. + ccvs.size()]                           ccvs branch currents
//   Rows [.. + inductors.size()] (DC only)            inductor branch currents
class Circuit
{
public:
	// Independent and dependent sources.
	vector<Voltage> vs;
	vector<Current> cs;
	vector<VCVS> vcvs;
	vector<VCCS> vccs;
	vector<CCVS> ccvs;
	vector<CCCS> cccs;

	// Linear passives.
	vector<Resistor> resistors;
	vector<Capacitor> caps;
	vector<Inductor> inductors;

	// Nonlinear active devices.
	vector<Diode> diodes;
	vector<BJT> bipolar_junction_transistors;

	AnalysisCommand analysis;

	// Number of external (netlist) nodes and internal (BJT expansion)
	// nodes. Used to size the MNA matrix.
	int external_nodes = 0;
	int internal_nodes = 0;

	// Labels for each row of the solution vector, built by
	// set_variable_names().
	vector<row_label_pair> variable_names;

	// Expands any multi-element devices (currently just BJTs) into
	// their primitive subcircuits. Must be called before set_indices()
	// and before any stamp method.
	void resolve_multi_element_devices()
	{
		for (int i = 0; i < bipolar_junction_transistors.size(); ++i)
			bipolar_junction_transistors[i].element_idx = i;

		for (auto& x : bipolar_junction_transistors)
			x.expand_device_model(vs, diodes, cccs, external_nodes, internal_nodes);
	}

	// Assigns a unique index to every element that needs an extra matrix
	// row of its own (voltage sources, VCVS, CCVS, inductors). The index
	// is used to locate the corresponding row in the matrix.
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

	// Builds the real MNA matrix and RHS for DC / OP analysis.
	//
	// The matrix size includes one row per external node, internal node,
	// voltage source, VCVS branch, CCVS branch, and inductor branch. The
	// matrix is assembled from triplets and handed to Eigen.
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

	// Builds the complex MNA matrix and RHS for AC analysis at the
	// starting frequency.
	//
	// The AC matrix has no inductor-current rows (inductors are stamped
	// as admittances), so its size differs from the DC matrix. The size
	// is consistent across the whole sweep.
	void buildAC(SparseMatrix<complex<double>>& G_matrix, VectorXcd& RHS_Ivector, double& fstart)
	{
		int matrix_size = external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size();

		G_matrix.resize(matrix_size, matrix_size);
		RHS_Ivector = VectorXcd::Zero(matrix_size);

		vector<Triplet <complex< double >> > matrix_initializer;
		matrix_initializer.reserve(5 * matrix_size);

		// AC stamping sources
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

		// AC stamping passives
		for (auto& x : resistors)
			x.stampAC(matrix_initializer);
		for (auto& x : caps)
			x.stampAC(matrix_initializer, fstart);
		for (auto& x : inductors)
			x.stampAC(matrix_initializer, fstart);

		// AC stamping actives
		for (auto& x : diodes)
			x.stampAC(matrix_initializer);

		G_matrix.setFromTriplets(matrix_initializer.begin(), matrix_initializer.end());
	}

	// Builds the real MNA matrix and RHS for the first transient step.
	//
	// Before stamping, every capacitor and inductor is initialized with
	// its voltage and current from the DC operating point, so the first
	// time step continues smoothly from the OP solution.
	void buildTRAN(SparseMatrix<double>& G_matrix, VectorXd& RHS_Ivector, VectorXd& dc_results)
	{
		int matrix_size = external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size();

		G_matrix.resize(matrix_size, matrix_size);
		RHS_Ivector = VectorXd::Zero(matrix_size);

		vector<Triplet<double>> matrix_initializer;
		matrix_initializer.reserve(5 * matrix_size);

		// Initialize cap state from the OP.
		for (auto& x : caps)
		{
			double V_plus = x.node_plus ? dc_results(x.node_plus - 1) : 0.0;
			double V_minus = x.node_minus ? dc_results(x.node_minus - 1) : 0.0;
			x.previous_voltage = V_plus - V_minus;
			x.previous_current = 0.0;
		}

		// Initialize inductor state from the OP branch current.
		for (auto& x : inductors)
		{
			x.previous_voltage = 0.0;
			x.previous_current = dc_results(external_nodes + internal_nodes + (int)vs.size() + (int)vcvs.size() + (int)ccvs.size() + x.element_idx);
		}

		// TRAN stamping sources
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

		// TRAN stamping passives
		for (auto& x : resistors)
			x.stampTRAN(matrix_initializer);
		for (auto& x : caps)
			x.stampTRAN(matrix_initializer, RHS_Ivector, analysis.tstep);
		for (auto& x : inductors)
			x.stampTRAN(matrix_initializer, RHS_Ivector, analysis.tstep);

		// TRAN stamping actives
		for (auto& x : diodes)
			x.stampTRAN(matrix_initializer, RHS_Ivector);

		G_matrix.setFromTriplets(matrix_initializer.begin(), matrix_initializer.end());
	}

	// Re-linearizes every diode around the latest solution and returns
	// the number of diodes whose junction voltage has converged. The
	// caller stops NR when this equals diodes.size().
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

	// In-place update of the frequency-dependent stamps (C and L) between
	// AC sweep steps.
	void updateAC_stamp(SparseMatrix< complex<double> >& G_matrix, VectorXcd& RHS_Ivector, const double& new_freq, const double& old_freq)
	{
		for (auto& x : caps)
			x.updateAC_stamp(G_matrix, RHS_Ivector, new_freq, old_freq);
		for (auto& x : inductors)
			x.updateAC_stamp(G_matrix, RHS_Ivector, new_freq, old_freq);
	}

	// Updates the RHS of the transient system between steps:
	//   - independent sources: incremental source value
	//   - capacitors / inductors: new companion-model I_eq
	//
	// The matrix itself does not change between steps, so nothing is
	// done to G_matrix here.
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

	// Builds the human-readable labels for every row in the solution
	// vector. Node voltages are "V(n)". Branch currents of named voltage
	// sources, VCVS, CCVS, and inductors are "I(name)". Elements whose
	// names start with "__" (internal BJT primitives) are skipped.
	void set_variable_names()
	{
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
