# Port-Controlled Hamiltonian Systems

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

## Examples of PCHS

Several electromechanical systems can be respresented as PCHS. A short list of these type of systems is presented below (click in the link to see the description):

* [A sincronous machine connected to an infinite bus](EX01.md)
* [RLC circuits](EX02.md)
* [Lossless power-system](EX03.md)
* [Second-order DC/DC converter](EX03.md)

More examples in the following link:

* [Examples](CH05.md)

## Passivity of PCHS

Passivity is a key property shared by many physical sistems (see [[Van-der-Schaft & Jeltsema 2014]](https://ieeexplore.ieee.org/document/8187102) for more details):  

**Definition:** A system $\dot{x} = f(x, u)$, $y = h(x, u)$, where $x \in \mathcal{X}\subseteq \mathbb{R}^n$ and $u, y \in \mathbb{R}^m$, is called *passive* if there exists a differentiable storage function $S : \mathcal{R} \rightarrow \mathbb{R}$ with $S(x) \geq 0, x \in \mathcal{X}$ , satisfying the differential dissipation inequality:

$\dot{S} \leq u^\top y$

This property is closely related to Lyapunov stability.  If $S$ has a strict minimum at a certain state $x^\star$, then $x^\star$ is a stable equlibrium point for the unforced dynamics (i.e., $u=0).    

A storage function for a PCHS can be defined as $S(x)=H(x)-H(x^\star)$ with $x^\star \in \operatorname{argmin} H(x)$. Therefore, the system is passive.

Passivity is preserved after parallel and feedback interconnection. These relevant properties are maintained in a natural manner using JuPCHS.  The manual below shows the use of the package.

## Manual

* [Defining PCHS](CH01.md)
* [Interconnecting PHCS](CH02.md)
* [Simulating for PCHS](CH03.md)
* [Analysing for PCHS](CH04.md)
* [Examples](CH05.md)
* [List of functions](CH06.md)
---

Next: [Defining PCHS](CH01.md)

