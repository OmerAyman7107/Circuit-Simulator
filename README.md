# Circuit-Simulator
This is a simple circuit simulation program that is able to perform DC, AC and transient analysis of linear and non linear circuits.

## How to use the program
The program prompts the user to enter a file address containg the netlist the user wants to simulate.
The syntax of the netlist is the same as a standerd spice netlist but with slight modifications.
### Netlist syntax 
The program has support for dependent and independent sources, resistors, capacitors, iductors, diodes and BJTs, the sytax for these elements are as follows:
- ***Resistor:*** Rname n+ n- value
- ***Capacitor:*** Cname n+ n- value
- ***Inductor:*** Lname n+ n- value
- ***Voltage:*** Vname n+ n- <[DC] value> <AC magnitude [phase]>
- ***Transient Voltage:*** Vname n+ n- SIN(VO VA fo TD a phase)
- ***Voltage controlled voltage source"*** Ename n+ n- nc1 nc2 value
- ***Voltage controlled current source:*** Gname n+ n- nc1 nv2 value
- ***Current controlled voltage source:*** Hname n+ n- Vcontrol value
- ***Current controlled current source:*** Fname n+ n- Vcontrol value
- ***Diode:*** Dname n+ n- model
- ***BJT:*** Qname C B E model

Some syntax isn't utilized as it should be like device models, because device parameters are already set to a default value and their isn't an option to choose a different device model other than changing the default values in the program.

### Analysis commands
The simulator is capable of performing DC, AC and Transient analysis, their syntax is listed below:
- ***DC analysis:*** .OP
- ***AC analysis:*** .AC <LIN/DEC/OCT> np fstart fend
- ***Transient analysis:*** .TRAN tstep tstop

In the next section we will simulate some example netlists.

## Simulation examples

