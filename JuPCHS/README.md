# JuPCHS
# A package for simulation, analysis and control of Port-Controlled Hamiltonian systems (PCHS)

## Introduction

Port-controlled Hamiltonian systems (PCHS) are a specialized class of dynamical systems that integrate concepts from Hamiltonian mechanics with port-based modeling, allowing for the analysis and control of complex physical systems.  Many electromechanical systems can be described by this paradigm, incluing, electric circuits and power systems, electric machines, hidraulic machines, power electronic converters and renewable energy sources. 

The dynamics of a PCHS can be expressed as follows:

$\dot{x} = (J(x)-R(x))\nabla H(x) + G(x) u(x)$

$y = G(x)^\top \nabla H(x)$

where 
* $x\in\mathbb{R}^{n}$ represents the state variables
* $u\in\mathbb{R}^{m}$ represents the control inputs
* $y\in\mathbb{R}^{m}$ represents the outputs
* $H:\mathbb{R}^n\rightarrow\mathbb{R}$ is the Hamiltonian function
* $\nabla H(x)$ represents the gradient of $H$
* $J(x)$ is a skew-symmetric interconection matrix
* $R(x)$ is a positive semidefinite disipation matrix
* $G(x)$ account for external influences

## What can JuPCHS do ?

JuPCHS is a package for simulation, analysis and control port-controlled Hamiltonian systems. JuPCHS make easy to create PCHS considering structural properties, namely:
* It checks positivity of $R$ and skew-symmetry of $J$. 
* It allows operations such as interconection and feedback.  
* It is easy to simulate the system using different methods for solving ordinary differential equations such as forward and backward Euler, Runge Kutta and MidPoint.
* Allows to easily find an equilibrium point
* It has specialized functions for sensitivity analysis
* It is easy to create passive-PID
* It is easy to develop model-predictive controls (MPC) 
  

## Geting started


### Installation

Just clone the project from github

    XXX

Open the folder in Visual-Studio Code.

### A simple example

A single sinchronous machine connected to an infinite bus can be represented as follows:

$M \Delta\dot{\omega} = p_m-p_e-\xi \Delta\omega$

$\dot{\delta} = \omega_B\Delta\omega$

where: $M$ is the inertia, $\omega_B$ the nominal frequency, $p_m$ the mechanical power, $p_e=p_\text{max}\sin(\delta)$ the electrical power, $\xi$ the damping, $\omega_b$ the nominal frequency, $\Delta\omega$ the frequency deviation and, $\delta$ the angle.

The dynamics of this system can be written under the framework of JuPCHS by defining state variables: $x=[x_\omega,x_\delta]$ with $x_\omega=M\omega$ and $x_\delta = \delta$. We only require to define $J$,$R$,$H$ and $\nabla H$ as given below:

    include("scr/JuPCHS.jl")
    using .JuPCHS
    Ms = 30
    wB = 2*pi*60
    pmax = 10
    Js = [0 -wB; wB 0]
    Rs = [188.5 0; 0 0]
    Hs(x) = 1/(2*Ms)*x[1]^2-pmax/wB*cos(x[2])
    dHs(x) = [x[1]/Ms; pmax/wB*sin(x[2])]

    SMIB = BuildPCHS(nx=2,nu=1,dt=1/60/10,
                     J=Js,R=Rs,H=Hs,dH=dHs)
    println(SMIB)

In this case, we have assumed damping $\xi=188.5$, $M=30$ and, $p_\text{max}=10$. Also, we assume mechanical power equal to zero (this may represent a synchronous condenser). JuPCHS returns some warnings indicating these assumptions. The function BuildPCHS returns a PCHS object named SMIB with all the information required for the simulation.  The println command returns the following output with the most important information:

    📌 build_PHS (Port Controlled Hamiltonian System)
        x' = (J(x)-R(x))∇H(x) + G(x)u(x)
    State variables:        [x1, x2]
    Inputs:                 [u1]
    Outputs:                [y1]
    Discretization time:    0.00833


Simulating is as simple as calling the function Simulate, namely:

    results = Simulate(PCHS=SMIB,nt=1000,xini=[0.1;0.1])
    
In this case, JuPCHS excecuted a simulation with 1000 points starting from xini = [0.1;0.1], and return a DataFrame named results.  The function PlotResults shows the state variables, the inputs and, the Hamiltonian, as given below:

    PlotResults(PCHS=SMIB,data=results)

![image](docs/PlotExample01.svg)
    
JuPCHS utilizes the backward-Euler as by default.  However, it is possible to change the solver to Forward Euler, Runge Kutta, Implicit Runge Kutta and, MidPoint, among others.  More details from this and other examples can be found in the tutorial (below).


     
### Tutorial

* [Introduction](docs/INTRO.md)
* [Defining Port-Controlled Hamiltonian Systems](docs/CH01.md)
* [Operations with PCHS that preserves structural properties](docs/CH02.md)
* [Simulation methods for PCHS](docs/CH03.md)
* [Analysis methods for PCHS](docs/CH04.md)

## Citing JuPCHS

If you find JuPCHS usefull for your research, we kindly request that you cite the following paper ()

 XXXXXXXXXXXXXXXXXX

---
## Contact

**Alejandro Garcés Ruiz**
Department of Electric Power Engineering
Universidad Tecnológica de Pereira
alejandro.garces@utp.edu.co
(https://github.com/alejandrogarces/)

## License

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC_BY--NC--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
