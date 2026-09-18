#include <iostream>
#include <string>
#include <fstream>
#include "CIRCUIT_ELEMENT.h"

using namespace std;

// Converts a string to upper case in place.
// Used to make keywords like "DC", "AC", "SIN", ".OP" case-insensitive.
void str_toUpper(string& input)
{
	for (auto& x : input)
		x = toupper(x);
}

// Forward declarations of the individual parse_* functions.
// Each one reads a single netlist line and appends the resulting element
// to the appropriate vector in the Circuit.
double read_value(const string& mag);
void parse_voltage(vector<Voltage>& vs, string& line, int& external_nodes);
void parse_current(vector<Current>& cs, string& line, int& external_nodes);
void parse_VCVS(vector<VCVS>& vcvs, const string& line, int& external_nodes);
void parse_VCCS(vector<VCCS>& vccs, const string& line, int& external_nodes);
void parse_CCVS(vector<CCVS>& ccvs, const string& line, int& external_nodes);
void parse_CCCS(vector<CCCS>& cccs, const string& line, int& external_nodes);
void parse_resistor(vector<Resistor>& resistors, const string& line, int& external_nodes);
void parse_capacitor(vector<Capacitor>& caps, const string& line, int& external_nodes);
void parse_inductor(vector<Inductor>& inductors, const string& line, int& external_nodes);
void parse_diode(vector<Diode>& diodes, const string& line, int& external_nodes);
void parse_BJT(vector<BJT>& bipolar_junction_transistors, const string& line, int& external_nodes, int& internal_nodes);
void parse_command(AnalysisCommand& analysis, const string& line);


// Reads a netlist file line by line and dispatches each line to the
// appropriate parse_* function based on its first character.
//
// The very first line is always skipped (treated as the netlist title,
// matching SPICE convention). Lines starting with '*' are treated as
// comments and ignored.
//
// After parsing, three post-processing steps are run:
//   - resolve_multi_element_devices(): expands BJTs into diodes, CCCSs,
//     and internal monitoring sources.
//   - set_indices(): assigns source_idx / element_idx to every element
//     that needs a matrix row of its own.
//   - set_variable_names(): builds human-readable labels for the solution
//     vector so results can be printed and exported.
void parser(const string& file_name, Circuit& circuit)
{
	ifstream netlist(file_name);
	if (netlist.fail())
	{
		cout << "failed to open file\n";
		return;
	}
	// skipping the first line in a netlist
	netlist.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	string line;
	while (getline(netlist, line))
	{
		// Dispatch based on the first character of the line.
		// SPICE uses the first letter of the element name as its type
		// designator (R, C, L, V, I, D, Q, E, G, H, F), and '.' for
		// dot-commands.
		switch (line[0])
		{
		case'v':
		case'V':
			parse_voltage(circuit.vs, line, circuit.external_nodes);
			break;

		case'i':
		case'I':
			parse_current(circuit.cs, line, circuit.external_nodes);
			break;

		case'e':
		case'E':
			parse_VCVS(circuit.vcvs, line, circuit.external_nodes);
			break;

		case'g':
		case'G':
			parse_VCCS(circuit.vccs, line, circuit.external_nodes);
			break;

		case'h':
		case'H':
			parse_CCVS(circuit.ccvs, line, circuit.external_nodes);
			break;

		case'f':
		case'F':
			parse_CCCS(circuit.cccs, line, circuit.external_nodes);
			break;

		case'r':
		case'R':
			parse_resistor(circuit.resistors, line, circuit.external_nodes);
			break;

		case'c':
		case'C':
			parse_capacitor(circuit.caps, line, circuit.external_nodes);
			break;

		case'l':
		case'L':
			parse_inductor(circuit.inductors, line, circuit.external_nodes);
			break;

		case'd':
		case'D':
			parse_diode(circuit.diodes, line, circuit.external_nodes);
			break;

		case'q':
		case'Q':
			parse_BJT(circuit.bipolar_junction_transistors, line, circuit.external_nodes, circuit.internal_nodes);
			break;

		case'.':
			// Dot-commands (.OP, .AC, .TRAN) set the analysis mode.
			parse_command(circuit.analysis, line);
			break;

		case'*':
			// Full-line comment; skip.
			continue;

		default:
			break;
		}
	}
	netlist.close();

	circuit.resolve_multi_element_devices();
	circuit.set_indices();
	circuit.set_variable_names();
	cout << "\nParsing completed successfuly.\n\n";
}


// Converts a numeric token with an optional single-character suffix into a
// double. The suffix is a SPICE-style scale factor:
//
//   T = 1e+12   G = 1e+9   M = 1e+6   k/K = 1e+3
//   m = 1e-3    u = 1e-6   n = 1e-9   p = 1e-12   f = 1e-15
//
// Values without a suffix are parsed directly. Suffixes are case-sensitive
// and single-character only, so "1M" is 1e6 while "1m" is 1e-3. Multi-char
// suffixes such as "1Meg" and scientific notation like "1e3" are NOT
// supported and will be misinterpreted.
double read_value(const string& mag)
{
	double result;
	if (!isdigit(mag[mag.size() - 1]))
	{
		char units[] = { 'T', 'G', 'M','k','K' , 'm', 'u', 'n', 'p', 'f' };
		double multiplier[] = { 1e+12, 1e+9, 1e+6, 1e+3, 1e+3, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15 };
		string temp = mag.substr(0, mag.length() - 1);
		int i;
		for (i = 0; i < 9; ++i)
			if (units[i] == mag[mag.size() - 1])
				break;
		result = stod(temp) * multiplier[i];
	}
	else
		result = stod(mag);
	return result;
}


///////////////////////////////////////////////////
///////////////*** PARSE SOURCES ***///////////////
///////////////////////////////////////////////////

// Parses a voltage source line and appends it to the vs vector.
//
// Supported forms (angle brackets denote optional tokens):
//   Vname n+ n- DC value
//   Vname n+ n- AC magnitude [phase]
//   Vname n+ n- value                        (treated as DC)
//   Vname n+ n- SIN(offset amplitude freq [delay [theta [phase]]])
//
// Note: the current implementation treats DC and SINE fields on the same
// source as mutually exclusive choices at stamp time. Mixing them on one
// line is legal syntactically but the SINE fields take precedence during
// transient analysis and the DC value is used during OP.
//
// Every '(' is replaced with a space before tokenizing so that the SINE
// parameters can be read as plain whitespace-separated tokens.
void parse_voltage(vector<Voltage>& vs, string& line, int& external_nodes)
{
	// remove any '(' if exist
	for (auto& x : line)
		if (x == '(') x = ' ';

	Voltage new_voltage_source;
	istringstream component_line(line);

	component_line >> new_voltage_source.name >> new_voltage_source.node_plus
		>> new_voltage_source.node_minus;
	string temp;

	while (component_line >> temp)
	{
		str_toUpper(temp);
		if (temp == "SIN" || temp == "SINE")
		{
			// SINE(...) waveform. Read everything up to ')' as a
			// whitespace-separated list of parameters. Missing
			// trailing parameters keep their default value of 0.
			new_voltage_source.is_tran = true;

			string vals;
			getline(component_line, vals, ')');
			istringstream param_values(vals);
			string val;

			param_values >> vals;
			new_voltage_source.sine_offset = read_value(vals);

			param_values >> val;
			new_voltage_source.sine_amplitude = read_value(val);

			param_values >> val;
			new_voltage_source.sine_frequency = read_value(val);

			param_values >> val;
			new_voltage_source.sine_delay = read_value(val);

			param_values >> val;
			new_voltage_source.sine_theta = read_value(val);

			param_values >> val;
			new_voltage_source.sine_phase = read_value(val);

			continue;
		}

		if (temp == "DC")
		{
			component_line >> temp;
			new_voltage_source.dc_magnitude = read_value(temp);
		}
		else if (temp == "AC")
		{
			// "AC mag phase" — phase is read directly as a double.
			// If phase is missing, cin-style extraction leaves it
			// unchanged (default 0).
			new_voltage_source.is_ac = true;
			component_line >> temp;
			new_voltage_source.ac_magnitude = read_value(temp);
			component_line >> new_voltage_source.phase;
		}
		else
		{
			// Bare value -> treat as DC magnitude.
			new_voltage_source.dc_magnitude = read_value(temp);
		}
	}

	vs.push_back(new_voltage_source);

	// Track the highest node number seen so far. The Circuit uses this
	// to size the MNA matrix.
	int pot_max = max(new_voltage_source.node_plus, new_voltage_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}


// Same as parse_voltage, but for current sources. The grammar is identical
// except that node roles are current-injection points rather than KVL
// branch terminals, so there is no extra matrix row for current sources.
void parse_current(vector<Current>& cs, string& line, int& external_nodes)
{
	for (auto& x : line)
		if (x == '(') x = ' ';

	Current new_current_source;
	istringstream component_line(line);

	component_line >> new_current_source.name >> new_current_source.node_plus
		>> new_current_source.node_minus;

	string temp;
	while (component_line >> temp)
	{
		str_toUpper(temp);
		if (temp == "SIN" || temp == "SINE")
		{
			// Same parameter list as the voltage source's SINE().
			new_current_source.is_tran = true;

			string vals;
			getline(component_line, vals, ')');
			istringstream param_values(vals);
			string val;

			param_values >> vals;
			new_current_source.sine_offset = read_value(vals);

			param_values >> val;
			new_current_source.sine_amplitude = read_value(val);

			param_values >> val;
			new_current_source.sine_frequency = read_value(val);

			param_values >> val;
			new_current_source.sine_delay = read_value(val);

			param_values >> val;
			new_current_source.sine_theta = read_value(val);

			param_values >> val;
			new_current_source.sine_phase = read_value(val);

			continue;
		}

		if (temp == "DC")
		{
			component_line >> temp;
			new_current_source.dc_magnitude = read_value(temp);
		}
		else if (temp == "AC")
		{
			new_current_source.is_ac = true;
			component_line >> temp;
			new_current_source.ac_magnitude = read_value(temp);
			component_line >> new_current_source.phase;
		}
		else
		{
			new_current_source.dc_magnitude = read_value(temp);
		}
	}
	cs.push_back(new_current_source);

	int pot_max = max(new_current_source.node_plus, new_current_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a voltage-controlled voltage source (VCVS):
//   Ename n+ n- nc+ nc- gain
// "n+ n-" are the output nodes; "nc+ nc-" are the controlling nodes.
// An extra matrix row/column is required (like a voltage source).
void parse_VCVS(vector<VCVS>& vcvs, const string& line, int& external_nodes)
{
	VCVS new_source;
	istringstream component_line(line);
	component_line >> new_source.name >> new_source.node_plus
		>> new_source.node_minus;
	component_line >> new_source.control_node_plus
		>> new_source.control_node_minus;
	string temp;
	component_line >> temp;
	new_source.voltage_gain = read_value(temp);
	vcvs.push_back(new_source);

	int pot_max = max(new_source.node_plus, new_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a voltage-controlled current source (VCCS):
//   Gname n+ n- nc+ nc- transconductance
// The output current is transconductance * V(nc+) - V(nc-).
// No extra matrix row is needed (it's a pure branch equation).
void parse_VCCS(vector<VCCS>& vccs, const string& line, int& external_nodes)
{
	VCCS new_source;
	istringstream component_line(line);
	component_line >> new_source.name >> new_source.node_plus
		>> new_source.node_minus;
	component_line >> new_source.control_node_plus
		>> new_source.control_node_minus;
	string temp;
	component_line >> temp;
	new_source.transconductance = read_value(temp);
	vccs.push_back(new_source);

	int pot_max = max(new_source.node_plus, new_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a current-controlled voltage source (CCVS):
//   Hname n+ n- Vcontrol transimpedance
// "Vcontrol" is the name of a voltage source through which the controlling
// current flows (the current variable of that source is used).
// Requires an extra matrix row, like a plain voltage source.
void parse_CCVS(vector<CCVS>& ccvs, const string& line, int& external_nodes)
{
	CCVS new_source;
	istringstream component_line(line);
	component_line >> new_source.name >> new_source.node_plus
		>> new_source.node_minus;
	component_line >> new_source.Vcontrol;
	string temp;
	component_line >> temp;
	new_source.transimpedence = read_value(temp);
	ccvs.push_back(new_source);

	int pot_max = max(new_source.node_plus, new_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a current-controlled current source (CCCS):
//   Fname n+ n- Vcontrol current_gain
// The output current is current_gain times the current through the named
// control voltage source.
void parse_CCCS(vector<CCCS>& cccs, const string& line, int& external_nodes)
{
	CCCS new_source;
	istringstream component_line(line);
	component_line >> new_source.name >> new_source.node_plus
		>> new_source.node_minus;
	component_line >> new_source.Vcontrol;
	string temp;
	component_line >> temp;
	new_source.current_gain = read_value(temp);
	cccs.push_back(new_source);

	int pot_max = max(new_source.node_plus, new_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}


//////////////////////////////////////////////////
//////////////*** PARSE PASSIVES ***//////////////
//////////////////////////////////////////////////

// Parses a resistor: Rname n+ n- value
void parse_resistor(vector<Resistor>& resistors, const string& line, int& external_nodes)
{
	Resistor new_passive;
	istringstream component_line(line);
	component_line >> new_passive.name >> new_passive.node_plus >> new_passive.node_minus;
	string temp;
	component_line >> temp;
	new_passive.value = read_value(temp);
	resistors.push_back(new_passive);

	int pot_max = max(new_passive.node_plus, new_passive.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a capacitor: Cname n+ n- value
// The cap's value is only used in AC and TRAN analyses; in DC it is an
// open circuit and has no stamp.
void parse_capacitor(vector<Capacitor>& caps, const string& line, int& external_nodes)
{
	Capacitor new_passive;
	istringstream component_line(line);
	component_line >> new_passive.name >> new_passive.node_plus >> new_passive.node_minus;
	string temp;
	component_line >> temp;
	new_passive.value = read_value(temp);
	caps.push_back(new_passive);

	int pot_max = max(new_passive.node_plus, new_passive.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses an inductor: Lname n+ n- value
// In DC the inductor is a short, so its stamp adds a branch-current row
// (like a zero-volt voltage source). In AC and TRAN it has its own stamps.
void parse_inductor(vector<Inductor>& inductors, const string& line, int& external_nodes)
{
	Inductor new_passive;
	istringstream component_line(line);
	component_line >> new_passive.name >> new_passive.node_plus >> new_passive.node_minus;
	string temp;
	component_line >> temp;
	new_passive.value = read_value(temp);
	inductors.push_back(new_passive);

	int pot_max = max(new_passive.node_plus, new_passive.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}


//////////////////////////////////////////////////
///////////////*** PARSE ACTIVES ***//////////////
//////////////////////////////////////////////////

// Parses a diode: Dname n+ n- model
// The model name is read but currently ignored — every diode uses the
// same Shockley parameters.
void parse_diode(vector<Diode>& diodes, const string& line, int& external_nodes)
{
	Diode new_active;
	istringstream component_line(line);
	component_line >> new_active.name >> new_active.node_plus
		>> new_active.node_minus >> new_active.model_name;

	diodes.push_back(new_active);

	int pot_max = max(new_active.node_plus, new_active.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a BJT: Qname C B E [model]
// The model name is optional and ignored. Each BJT contributes 5 internal
// nodes (used by the Ebers-Moll expansion in Circuit::resolve_multi_element_devices),
// which the caller accounts for by incrementing internal_nodes.
void parse_BJT(vector<BJT>& bipolar_junction_transistors, const string& line, int& external_nodes, int& internal_nodes)
{
	BJT new_active;
	istringstream component_line(line);
	component_line >> new_active.name >> new_active.collector_node
		>> new_active.base_node >> new_active.emitter_node;

	bipolar_junction_transistors.push_back(new_active);
	internal_nodes += 5;
	int pot_max = max(new_active.base_node, max(new_active.collector_node, new_active.emitter_node));
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

// Parses a dot-command (.OP, .AC, .TRAN) and fills the AnalysisCommand
// struct with the requested analysis and its parameters.
//
//   .OP
//   .AC LIN|DEC|OCT np fstart fstop
//   .TRAN tstep tstop
//
// The type_line field stores a short human-readable description used in
// the exported result files' header.
void parse_command(AnalysisCommand& analysis, const string& line)
{
	istringstream command_line(line);
	string temp;
	command_line >> temp;
	str_toUpper(temp);

	if (temp == ".OP")
	{
		analysis.type = AnalysisType::OP;
		analysis.type_line = string() + "OP";
	}
	else if (temp == ".AC")
	{
		// Second token selects the sweep kind.
		string sweep_kind;
		command_line >> sweep_kind;
		str_toUpper(sweep_kind);
		analysis.type_line = string() + "AC";

		if (sweep_kind == "LIN") { analysis.type = AnalysisType::AC_LIN; analysis.type_line += " LIN"; }
		else if (sweep_kind == "DEC") { analysis.type = AnalysisType::AC_DEC; analysis.type_line += " DEC"; }
		else if (sweep_kind == "OCT") { analysis.type = AnalysisType::AC_OCT; analysis.type_line += " OCT"; }

		// np = points, fstart = start frequency, fstop = stop frequency.
		string points_str, fstart_str, fstop_str;
		command_line >> points_str >> fstart_str >> fstop_str;
		analysis.points = (int)read_value(points_str);
		analysis.freq_start = read_value(fstart_str);
		analysis.freq_stop = read_value(fstop_str);
	}
	else if (temp == ".TRAN")
	{
		analysis.type = AnalysisType::TRAN;
		string tstep_str, tstop_str;
		command_line >> tstep_str >> tstop_str;
		analysis.tstep = read_value(tstep_str);
		analysis.tstop = read_value(tstop_str);

		analysis.type_line = string() + "TRAN";
	}
}
