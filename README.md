# mesoCells repository
By Jack Treado, Yale University

## Getting started

The code for this project is all contained in the `src` directory. It is written completely in `MATLAB`, so everything should work out of the box. 

The main file that runs the simulation is the function `runMesoCellExtension.m`. This is a `MATLAB` function that you can run from the command
line in the MATLAB editor screen. 

The function is defined as follows:
`function [hList, xList, yList, shapeList, calAList] = runMesoCellExtension(NV,NPINS,Kl,Kb,lambdaA,lambdaB,plotIt,varargin)`
* `NV`: integer number of vertices that make up the deformable particle (DP)
* `NPINS`: integer number of vertices that make up the DP. 
	* **NOTE**: MUST be at least $\leq$ than `NV`
