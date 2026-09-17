# Circuit-Simulator
This is a simple circuit simulation program that is able to perform DC, AC and transient analysis of linear and non linear circuits.

## How to use the program
The program prompts the user to enter a file address containg the netlist the user wants to simulate.
The syntax of the netlist is the same as a standerd spice netlist but with slight modifications.
The output is written to a file 
### Netlist syntax 
The program has support for dependent and independent sources, resistors, capacitors, iductors, diodes and BJTs, the sytax for these elements are as follows:
- ***Resistor:*** `Rname n+ n- value`
- ***Capacitor:*** `Cname n+ n- value`
- ***Inductor:*** `Lname n+ n- value`
- ***Voltage:*** `Vname n+ n- <[DC] value> <AC magnitude [phase]>`
- ***Transient Voltage:*** `Vname n+ n- SIN(VO VA fo TD a phase)`
- ***Voltage controlled voltage source"*** `Ename n+ n- nc1 nc2 value`
- ***Voltage controlled current source:*** `Gname n+ n- nc1 nv2 value`
- ***Current controlled voltage source:*** `Hname n+ n- Vcontrol value`
- ***Current controlled current source:*** `Fname n+ n- Vcontrol value`
- ***Diode:*** `Dname n+ n- model`
- ***BJT:*** `Qname C B E model`

Some syntax isn't utilized as it should be like device models, because device parameters are already set to a default value and their isn't an option to choose a different device model other than changing the default values in the program.

### Analysis commands
The simulator is capable of performing DC, AC and Transient analysis, their syntax is listed below:
- ***DC analysis:*** `.OP`
- ***AC analysis:*** `.AC <LIN/DEC/OCT> np fstart fend`
- ***Transient analysis:*** `.TRAN tstep tstop`

In the next section we will simulate some example netlists.

## Simulation examples
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

