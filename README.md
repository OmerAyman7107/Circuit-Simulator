# Circuit-Simulator

A simple circuit simulation program that performs DC, AC, and transient
analysis of linear and nonlinear circuits.

## How to use the program

The program prompts the user to enter a file path containing the netlist to
simulate. The netlist syntax follows standard SPICE conventions, though it
does not offer the same flexibility and customizability as a full SPICE
parser.

Results are written to a text file and can be plotted with a simple script,
as shown in the examples section.

### Netlist syntax

The program supports dependent and independent sources, resistors,
capacitors, inductors, diodes, and BJTs. The syntax for each element is as
follows:

- **Resistor:** `Rname n+ n- value`
- **Capacitor:** `Cname n+ n- value`
- **Inductor:** `Lname n+ n- value`
- **Voltage source:** `Vname n+ n- <[DC] value> <AC magnitude [phase]>`
- **Transient voltage source:** `Vname n+ n- SIN(VO VA fo TD a phase)`

  `VO` = offset, `VA` = amplitude, `fo` = frequency, `TD` = delay,
  `a` = damping, `phase` = phase in degrees.

- **Voltage-controlled voltage source:** `Ename n+ n- nc1 nc2 value`
- **Voltage-controlled current source:** `Gname n+ n- nc1 nc2 value`
- **Current-controlled voltage source:** `Hname n+ n- Vcontrol value`
- **Current-controlled current source:** `Fname n+ n- Vcontrol value`
- **Diode:** `Dname n+ n- model`
- **BJT:** `Qname C B E model`

Some fields are parsed but not used. Device models, for example, are read
and ignored: every device uses a fixed default model whose parameters can
only be changed in the source code, and the model name itself has no effect
on the simulation.

### Analysis commands

- **DC analysis:** `.OP`
- **AC analysis:** `.AC <LIN/DEC/OCT> np fstart fend`
- **Transient analysis:** `.TRAN tstep tstop`

The next section walks through a few example netlists.

## Simulation examples

Example 1: full bridge rectifier with filter capacitor:

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

The output is exported to a text file. Plotting it with a short Python script produces the waveform below.

<img width="1920" height="974" alt="Figure_1" src="https://github.com/user-attachments/assets/b7209697-a76f-437e-bb94-97e2e4add681" />

Example 2: simple band-pass filter:
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

Band-pass filter frequency response:

<img width="1920" height="974" alt="Figure_2" src="https://github.com/user-attachments/assets/69a7c30f-1236-4ff8-98b7-674268ef7893" />

Example 3: simple BJT biasing circuit:
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

The simulator performs three main functions:

1. Parse the netlist
2. Build the circuit data structure
3. Solve the circuit

Each is discussed below.

### Parser

The netlist parser reads space-separated tokens and uses them to build the
circuit data structure.

Scale factors are case-sensitive and limited to a single character. Unlike
standard SPICE syntax, the following suffixes are recognized:

- `T` → 1e+12
- `G` → 1e+9
- `M` → 1e+6
- `K` or `k` → 1e+3
- `m` → 1e-3
- `u` → 1e-6
- `n` → 1e-9
- `p` → 1e-12
- `f` → 1e-15

### Circuit data structure

The circuit data structure is a class that holds everything the solver
needs: the elements of the circuit, their positions and values, the number
of nodes, and the mode of analysis. It also provides helper methods that
build and update the Modified Nodal Analysis (MNA) matrix using element
stamps. The stamps themselves change with the analysis type; capacitors
and inductors contribute different stamps in transient versus AC, and
diodes contribute nonlinear stamps that are updated on every
Newton-Raphson iteration.

### Solver

The solver class determines the analysis type and calls the circuit's
helper method to build the MNA matrix, then solves it using sparse LU
factorization.

## Numerical methods

### Newton-Raphson method

The Newton-Raphson method solves nonlinear equations iteratively by
linearizing around an initial guess and refining the solution until it
converges. Here it is used to linearize nonlinear elements - diodes and
the internal junctions of BJTs - into a linear companion model that is
solved as part of the linear system. The process repeats until the
solution stops changing between iterations.

### Trapezoidal method

The trapezoidal method approximates the area under a curve as a trapezoid.
In transient analysis, it is used to derive linear companion models for
capacitors and inductors, turning the circuit's differential equations
into a system of linear equations that can be solved step by step.

### Sparse matrices

Sparse matrices are matrices in which most entries are zero. Storing only
the nonzero entries reduces both the memory footprint and the time
required for operations such as factorization, an advantage that grows
with the size of the circuit.

## Known limitations

This simulator was written as a learning project, and its scope is
deliberately smaller than that of a standard SPICE simulator. The
following limitations are worth knowing before you use it.

### Analyses

- Only one analysis command per netlist. If more than one (`.OP`, `.AC`,
  `.TRAN`) is included, only the last one is used.
- AC analysis is small-signal only, and only sources marked with the `AC`
  keyword contribute. Noise, distortion, and parameter sweeps are not
  supported.

### Netlist syntax

- Only `*` at the start of a line is treated as a comment. Inline comments
  are not supported.
- Line continuation is not supported.
- Node names must be integers. Node `0` is ground.
- Scale factors are single-character and case-sensitive. Multi-character
  suffixes such as `1Meg`, and scientific notation such as `1e3`, are not
  accepted; write `1M` and `1k` instead.
- `.MODEL`, `.SUBCKT`, `.PARAM`, `.INCLUDE`, and other dot-commands are
  not supported. Only `.OP`, `.AC`, and `.TRAN` are recognized.
- Transient sources support only the `SIN`/`SINE` waveform. `PULSE`, `EXP`,
  and `PWL` are not implemented.

### Device models

- Device parameters are hardcoded and cannot be changed from the netlist.
  The model name in `Dname n+ n- model` and `Qname C B E model` is parsed
  but ignored.
- The diode model is the ideal Shockley equation, with no series
  resistance, junction capacitance, breakdown voltage, or temperature
  dependence.
- The BJT model is a basic Ebers-Moll model with fixed forward and reverse
  alphas. It has no Early effect, junction capacitance, or temperature
  dependence.
- MOSFETs, switches, transmission lines, coupled inductors, and other
  elements beyond those listed above are not supported.

### Output

- Results are written to a single text file, one variable at a time. There
  is no batch export or multi-variable output.
- There is no built-in plotting. Waveforms are exported as plain text and
  must be plotted separately, for example, with Python and Matplotlib, as
  shown in the examples.

## Acknowledgements

I would like to thank my friend, whose fruitful discussions and
brainstorming sessions helped shape this project from the beginning. I
would also like to thank Professor Hesham Omran, whose lectures on SPICE
simulators were invaluable in helping me understand how these tools work
under the hood. First and foremost, thanks be to God; this would not have
been possible without His blessing.
