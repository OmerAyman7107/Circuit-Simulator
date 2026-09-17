# Circuit-Simulator
This is a simple circuit simulation program that is able to perform DC, AC, and transient analysis of linear and non linear circuits.

## How to use the program
The program prompts the user to enter a file address containing the netlist the user wants to simulate.
The syntax of the netlist is the same as a standard spice netlist but it doesn't provide the same flexibility and customizability of a standard spice netlist parser.
The output is written to a text file and can easily be plotted using a simple script as will be seen in the examples section. 
### Netlist syntax 
The program has support for dependent and independent sources, resistors, capacitors, inductors, diodes, and BJTs, the syntax for these elements are as follows:
- ***Resistor:*** `Rname n+ n- value`
- ***Capacitor:*** `Cname n+ n- value`
- ***Inductor:*** `Lname n+ n- value`
- ***Voltage:*** `Vname n+ n- <[DC] value> <AC magnitude [phase]>`
- ***Transient Voltage:*** `Vname n+ n- SIN(VO VA fo TD a phase)`

  `VO` = offset, `VA` = amplitude, `fo` = frequency, `TD` = delay, `a` = damping, `phase` = phase in degrees.
- ***Voltage controlled voltage source*** `Ename n+ n- nc1 nc2 value`
- ***Voltage controlled current source:*** `Gname n+ n- nc1 nc2 value`
- ***Current controlled voltage source:*** `Hname n+ n- Vcontrol value`
- ***Current controlled current source:*** `Fname n+ n- Vcontrol value`
- ***Diode:*** `Dname n+ n- model`
- ***BJT:*** `Qname C B E model`

Some syntax isn't utilized as it should be like device models, because device parameters are already set to a default value and there isn't an option to choose a different device model other than changing the default values in the program, in fact you can even not write the dveice model and it will work just fine as it has no actual use in the program.

### Analysis commands
The simulator is capable of performing DC, AC and Transient analysis, their syntax is listed below:
- ***DC analysis:*** `.OP`
- ***AC analysis:*** `.AC <LIN/DEC/OCT> np fstart fend`
- ***Transient analysis:*** `.TRAN tstep tstop`

In the next section we will simulate some example netlists.

## Simulation examples
Example 1:

```
  Full bridge rectifier circuit with filter capacitor
  V1 1 3 SIN(0 100 1k 1m 0 90)
  D1 1 2
  D2 0 1
  D3 3 2
  D4 0 3
  R1 0 2 100k
  C1 0 2 10u
  .tran 10u 10m
```
<p align="center">
<img width="551" height="316" alt="image" src="https://github.com/user-attachments/assets/448ae850-042d-4524-bef4-105456b91067" />
</p>
The output is exported to a text file and by plotting the result using a simple python script we get this waveform.

<img width="1920" height="974" alt="Figure_1" src="https://github.com/user-attachments/assets/b7209697-a76f-437e-bb94-97e2e4add681" />

Example 2:
```
Simple bandpass filter circuit
V1 1 0 AC 1
L1 1 2 10u
C1 3 2 25p
R1 3 0 50
.ac dec 100 1 10G
```
<p align="center">
<img width="550" height="310" alt="image" src="https://github.com/user-attachments/assets/31ed8ef9-87ac-488a-803d-efb5d013e9d8" />
</p>

The output is exported to a text file and by plotting the result using a simple python script we get this waveform.

<img width="1920" height="974" alt="Figure_2" src="https://github.com/user-attachments/assets/69a7c30f-1236-4ff8-98b7-674268ef7893" />

Example 3:
```
A simple BJT biasing circuit
R1 1 2 10k
R2 5 0 1k
R3 4 3 100k
Q1 2 4 5 
V1 1 0 12
V2 3 0 5
.op
```
<p align="center">
<img width="674" height="413" alt="image" src="https://github.com/user-attachments/assets/76cf6165-0205-4089-b94a-1c260d016f95" />
</p>


## How the simulator works
The simulator performs 3 main functions:
<ol>
  <li>Parse the netlist</li>
  <li>Build the circuit data structure</li>
  <li>Solve the circuit</li>
</ol>

In the following sections we will discuss the these three in detail.

### Parser
The netlist parser is a function that works by reading space separated tokens and storing them to build the circuit data structure.
Some compromises were made when reading scale factors, they are case sensitive and can only be single characters. Unlike standard spice syntax, the scaling factors are provided below to avoid confusion and errors when using the program:
 - T -> E+12
 - G -> E+9
 - M -> E+6
 - K or k -> E+3
 - m -> E-3
 - u -> E-6
 - n -> E-9
 - p -> E-12
 - f -> E-15

### Circuit data structure
The circuit data structure is a class containing all the information and helper methods that the solver might need to solve the circuit, it stores the elements of the circuit and their positions in the circuit and their values, it also stores the number of nodes of the circuit and information about the mode of analysis. The circuit class also provides helper that build and update the modified nodal analysis matrix using element stamps, the element stamps may change according to the type of analysis like capacitors and inductors or diodes.

### Solver
The solver class is used to solve the matrix by determining the mode of analysis and then solving the circuit. To do this, it calls the circuit's helper method that builds the modified nodal analysis matrix then solves it using sparse LU factorization.

## Numerical methods
This section discusses the numerical methods that were used in the program.

### Newton-Raphson method
The Newton-Raphson method is an iterative method used to solve nonlinear equations by linearizing around an initial guess then solving iterativly until the solution converges, it is used for linearizing nonlinear circuit elements and finding a linearized companion model to be replaced with the element then solve the system until the solution converges.

### Trapezoidal method
The Trapezoidal method is a numerical integration method that approximates the area under the graph of the function as a trapezoid, this method can be used in transient analysis to derive a linear companion model for capacitors and inductors  to solve the system of differential equations iteratively at every time step.

### Sparse matrices
Sparse matrices are matrices that are mostly populated with zeros, this property can be used to reduce the memory taken by large matrices and reduce the time taken to perform operations on it.

## Known limitations

This simulator was written as a learning project, and its scope and capabilities are deliberately smaller than a standard SPICE simulator. The following limitations are worth knowing before you use it.

### Analyses
- Only one analysis command per netlist. If you include more than one (`.OP`, `.AC`, `.TRAN`), only the last one is used.
- AC analysis is small‑signal only, and only sources marked with the `AC` keyword contribute. There is no noise analysis, distortion analysis, or parameter sweep.

### Netlist syntax
- Only `*` at the start of a line is treated as a comment. Inline comments are not supported.
- Line continuation is not supported.
- Node names must be integers. Node `0` is ground.
- Scale factors are single‑character only and case‑sensitive. Multi‑character suffixes such as `1Meg` or scientific notation such as `1e3` are not accepted; write `1M` and `1k` instead.
- There is no support for `.MODEL` cards, `.SUBCKT`, `.PARAM`, `.INCLUDE`, or any other dot‑command other than `.OP`, `.AC`, and `.TRAN`.
- Transient sources support only the `SIN`/`SINE` waveform. `PULSE`, `EXP`, and `PWL` are not implemented.

### Device models
- Device parameters are hardcoded and cannot be changed from the netlist. `Dname n+ n- model` and `Qname C B E model` parse a model name, but it is ignored.
- The diode model is the ideal Shockley equation. It has no series resistance, no junction capacitance, no breakdown voltage, and no temperature dependence.
- The BJT model is a basic Ebers‑Moll model with fixed forward and reverse alphas. It has no Early effect, no junction capacitance, and no temperature dependence.
- There is no support for MOSFETs, switches, transmission lines, coupled inductors, or any other element not listed in the syntax section.

### Output
- The simulator writes results to a single text file chosen by the user, one variable at a time. There is no batch export or multi‑variable output.
- There is no built‑in plotting. Waveforms are exported as plain text and must be plotted separately (for example, with Python and Matplotlib, as shown in the examples).

As you can see there are many limitations to the simulator, and I was lazy enough that I had to let an Ai list them all for me instead of doing it myself.

At the end I would like to thank my friend who helped me throughout this project for the fruitful discussions and brainstorming sessions we had, and I would also like to thank rofessor Hesham Omran who's lectures were really extremly useful in understanding spice simultors. None of this could have been possible without them (after god's blessing firstly).
