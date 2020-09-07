# Ultracold_Atoms_src
## Author:
Tomás Sánchez Sánchez-Pastor, Grupo de Sistemas Complejos UPM

## 1. Getting started
Firstly, you have to add an ssh key to your github account. To do so, just execute on the terminal:
```
$ ssh-keygen
```
And press enter to everything.

Then, copy the content on .ssh/id_rsa.pub into your ssh keys on GitHub.

Once you did that you only need to run:
```
$ git clone git@github.com:TsspGit/CodVid19_CIEMAT.git
```
And the project would be copied to your current folder.

## 2. Folders explanation
### src/: 
In this folder can be found the more general notebooks on ultracold atoms theory such as quasi-1D Harmonic Oscillator + Delta potential spectrum or a toy model of Landau Zener transitions. 

### utils/: 
- **atomic_units.py** contains the values of the fundamental constants of the Hartee atomic units.
- **VjtoIj.py** python script to input the data as written in the papers and output the data as must be written in the codes. Run: ./VjtoIj.py and the code start to ask the parameters. Output: Intensities in W/cm2 and the trap length.
- **Units.ipynb** a python notebook seemly to VjtoIj.py. It contains the equations in LaTeX format.
- **General_Figures** Figures of the src/ notebooks.

### convergence/:
- **Convergence_With_Potential.ipynb** this notebook uses two useful functions to read and plot the eigenvalues eassily, furthermore it presents an example of a well-convergenced execution.

### out/:
Saved logs of otaold_2hm.csh and the scatlength vs potential parameter.

### potentials/:
- **potutils/** contains a function to read the potentials eassily.
- **Li7Li7_potential071.ipynb** an example of Li7-Li7 potential with parameter 071.