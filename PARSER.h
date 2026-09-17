#include <iostream>
#include <string>
#include <fstream>
#include "CIRCUIT_ELEMENT.h"

using namespace std;

void str_toUpper(string& input)
{
	for (auto& x : input)
		x = toupper(x);
}

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
			parse_command(circuit.analysis, line);
			break;
		case'*':
			continue;
			break;

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


double read_value(const string& mag)
{
	double result;
	//cout << mag << endl;
	if (!isdigit(mag[mag.size() - 1]))
	{
		char units[] = { 'T', 'G', 'M','k','K' , 'm', 'u', 'n', 'p', 'f' };
		double multiplier[] = { 1e+12, 1e+9, 1e+6, 1e+3, 1e+3, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15 };
		string temp = mag.substr(0, mag.length() - 1);
		int i;
		for (i = 0; i < 9; ++i)
			if (units[i] == mag[mag.length() - 1])
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

// the function only parses AC and DC voltage sources
// still doesn't support PULSE, SINE, EXP or PWL that are used in transient anlysis
void parse_voltage(vector<Voltage>& vs, string& line, int& external_nodes)
{
	// remove any '(' if exist
	for (auto& x : line)
		if (x == '(') x = ' ';
	//cout << line << endl;

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
			new_voltage_source.is_tran = true;
			//if (false)
			//{
			//	string val = temp.substr(4, temp.size() - 4);
			//	new_voltage_source.sine_offset = read_value(val);
			//}
			//else
			//{
			//	component_line >> temp;
			//	new_voltage_source.sine_offset = read_value(temp);
			//}
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
			new_voltage_source.is_ac = true;
			component_line >> temp;
			new_voltage_source.ac_magnitude = read_value(temp);
			component_line >> new_voltage_source.phase;
		}
		else
		{
			new_voltage_source.dc_magnitude = read_value(temp);
		}
	}

	vs.push_back(new_voltage_source);

	int pot_max = max(new_voltage_source.node_plus, new_voltage_source.node_minus);
	if (external_nodes < pot_max)
		external_nodes = pot_max;
}

void parse_current(vector<Current>& cs, string& line, int& external_nodes)
{
	// remove any '(' if exist
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
			new_current_source.is_tran = true;
			//if (false)
			//{
			//	string val = temp.substr(4, temp.size() - 4);
			//	new_current_source.sine_offset = read_value(val);
			//}
			//else
			//{
			//	component_line >> temp;
			//	new_current_source.sine_offset = read_value(temp);
			//}
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

// Ename n+ n- nc+ nc- gain
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

// Gname n+ n- nc+ nc- transconductance
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

// Hname n+ n- Vcontrol transimpedence
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

// Fname n+ n- Vcontrol current_gain
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

// Rname n+ n- value
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

// Cname n+ n- value
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

// Lname n+ n- value
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

// Dname n+ n- model_name
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

// Qname C B E model_name
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

// The parse command function is used to parse the analysis mode
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
		string sweep_kind;
		command_line >> sweep_kind;
		str_toUpper(sweep_kind);
		analysis.type_line = string() + "AC";

		if (sweep_kind == "LIN") { analysis.type = AnalysisType::AC_LIN; analysis.type_line += " LIN"; }
		else if (sweep_kind == "DEC") { analysis.type = AnalysisType::AC_DEC; analysis.type_line += " DEC"; }
		else if (sweep_kind == "OCT") { analysis.type = AnalysisType::AC_OCT; analysis.type_line += " OCT"; }

		string points_str, fstart_str, fstop_str;
		command_line >> points_str >> fstart_str >> fstop_str;
		analysis.points = (int)read_value(points_str);
		analysis.freq_start = read_value(fstart_str);
		analysis.freq_stop = read_value(fstop_str);

		//analysis.type_line += " " + to_string(analysis.points) + " " + to_string(analysis.freq_start) + " " + to_string(analysis.freq_stop);
	}
	else if (temp == ".TRAN")
	{
		analysis.type = AnalysisType::TRAN;
		string tstep_str, tstop_str;
		command_line >> tstep_str >> tstop_str;
		analysis.tstep = read_value(tstep_str);
		analysis.tstop = read_value(tstop_str);

		analysis.type_line = string() + "TRAN";
		//analysis.type_line += " " + to_string(analysis.tstep) + " " + to_string(analysis.tstop);
	}
}