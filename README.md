# Yamada Model - Numerical Continuation Scripts

This repository contains the code used to calculate the bifurcation diagram and phase resetting results for the Yamada model of a Q-switching laser [1,2]

$$ \dot{G} = \gamma ( A - G - G I ) ,$$

$$ \dot{Q} = \gamma ( B - Q - a Q I ) ,$$

$$ \dot{I} = ( G - Q - 1 ) I ,$$

where $G$ is the gain, $Q$ the absorption, and $I$ the intensity.

This repository essentially contains two programs: `bifurcations` and `phase_resetting`. 
- `bifurcations` contains the code which calculates and generates the two-parameter bifurcation diagram in $A$ and $\gamma$. 
- `phase_resetting` contains code which finds a stable periodic orbit, and then computes some "phase resetting" calculations. 

The code for both of these programs is written in either [AUTO-07P](https://www.github.com/auto-07p/auto-07p/) (in Python and Fortran) or [Continuation Core and Toolboxes (COCO)](https://sourceforge.net/projects/cocotools/) (in Matlab). 

As I learnt COCO before learning AUTO, I have written the AUTO code to resemble that of COCO.

## References
[1] M. Yamada. "A theoretical analysis of self-sustained pulsation phenomena in narrow-stripe semiconductor lasers". *IEEE J. Quantum Electron* **29**, 1330 (1993).

[2] J. L. A. Dubbeldam and B. Krauskopf. "Self-pulsations of lasers with saturable absorbers: Dynamics and bifurcations". *Opt. Commun.* **159**, 325 (1999).