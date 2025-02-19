# Mathematica-RL
Packages for running RL within mathematica

This repository contains the code for Reinforcement Learning (RL), 
implemented in Mathematica. This was intended for use in the paper:
https://doi.org/10.1002/prop.202200034.

Both REINFROCE and A2C are implemented

## Installation

It's possible to use the package without installing, if you have access to the 
internet. Simply run
```
Get["https://raw.githubusercontent.com/harveyThomas4692/Mathematica-RL/\
main/ActorCritic.m"]
```
or 
```
Get["https://raw.githubusercontent.com/harveyThomas4692/Mathematica-RL/\
main/REINFORCE.m"]
```
in the Mathematica Kernel.

Note that importing both will lead to name conflicts.

Alternatively, if you wish to use the package without internet access, you will need
to install it. To install the package, please copy "Genetic.m" into your 
Mathematica "Applications" folder. This can be found in 
"Home_Directory/.Mathematica/Applications", where your home directory can
be found by running
```
$HomeDirectory

```
in the Mathematica kernel. The package can then be imported by running the line
```
Get["ActorCritic`"]
```
```
Get["REINFORCE`"]
```
Again, note that importing both will lead to name conflicts.
## Examples
There is currently one example notebook for eac method contained in the repository.

## Use
Documentation is built into the package. This can be accessed by running the
following lines in Mathematica:
```
?ActorCritic
```
```
?ActorCriticModules
```
for actor critic. Or 
```
?REINFORCE
```
```
?REINFORCEModulesModules
```
for REINFORCE.

Similarly, any function listed in this last command has their own documentation.
This can be found by running 
```
?"FunctionName"
```
I also suggest looking at the examples, as these demonstrate how one can connect
custom functions to this package.

## Citation
If you use this code, please cite the following bib entries:

```
@article{RLPackage,
  author = "Harvey, Thomas R. and Lukas, Andre",
  title = "{Mathematica package for Reinforcement Learning}",
  url="https://github.com/harveyThomas4692/Mathematica-RL",
  note="https://github.com/harveyThomas4692/Mathematica-RL"
}
@article{Harvey:2021oue,
    author = "Harvey, T. R. and Lukas, A.",
    title = "{Quark Mass Models and Reinforcement Learning}",
    eprint = "2103.04759",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    doi = "10.1007/JHEP08(2021)161",
    journal = "JHEP",
    volume = "08",
    pages = "161",
    year = "2021"
}
@article{Constantin:2021for,
    author = "Constantin, Andrei and Harvey, Thomas R. and Lukas, Andre",
    title = "{Heterotic String Model Building with Monad Bundles and Reinforcement Learning}",
    eprint = "2108.07316",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    doi = "10.1002/prop.202100186",
    journal = "Fortsch. Phys.",
    volume = "70",
    number = "2-3",
    pages = "2100186",
    year = "2022"
}
```
