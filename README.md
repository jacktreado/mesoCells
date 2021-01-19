# mesoCells repository
By Jack Treado, Yale University

## Downloading code

You can download this whole repository to get access to the code.

The code is all located in the `src` directory. You can only run the code while the `src` directory is your working directory in MATLAB.  

## Getting started

The main file that runs the simulation is the function `runMesoCellExtension.m`. This is a MATLAB function that you can run from the command
line in the MATLAB editor screen. 

There is another function, `vertexFIREPinnedVerts.m`, included in the `src` directory. This is needed during the simulation, so make sure it always stays in the same working directory as
`runMesoCellExtension.m`. 

The main function is defined as follows:

`function [hList, xList, yList, shapeList, calAList] = runMesoCellExtension(NV,NPINS,Kl,Kb,lambdaA,lambdaB,plotIt,movieFileStr)`

### Inputs:
* `NV`: integer number of vertices that make up the deformable particle (DP)
* `NPINS`: integer number of pins to force the DP to change shape during extension. 
	* **NOTE**: MUST be at least _less than_ than `NV`
* `Kl`: mechanical constant for perimeter
* `Kb`: mechanical constant for curvature
* `lambdaA`: rate constant for growth/decay of the preferred shape parameter `calA0` based on the instantaneous shape parameter `calA`
* `lambdaB`: rate constant for growth of `Kb` during stretching simulation
	* **NOTE**: Subject to change, form of Kb growth is current pure exponential, may want to change
* `plotIt`: Binary var. to either draw (`1`) or not draw (`0`) cells during simulation

#### Optional Input:
* `movieFileStr`: If included, save animation to movie file. 
	* **NOTE**: File will be saved as an .mp4 file, currently does not support other formats.

### Outputs:
* `hList`: List of pin stretch steps parametrized by variable <img src="https://render.githubusercontent.com/render/math?math=h">. 
	* **NOTE**: Negative entries mean pins extended out from original cell boundary.
* `xList`: _x_-coordinates of all cell vertices at each stretch step
* `yList`: _y_-coordinates of all cell vertices at each stretch step
	* **NOTE**: Coordinates are stored as MATLAB `cell` variables. To access at step `ii`, use syntax `x{ii}` or `y{ii}`.
* `shapeList`: List of relevant shape information during simulation. Column 1 is `calA0`, column 2 is `Kb.`
* `calAList`: Instantaneous shape parameter <img src="https://render.githubusercontent.com/render/math?math=\mathcal{A} = p^2/4\pi a"> during stretching simulation.

## Running the code

To run the code on the command line, simply use

`>> [hList, xList, yList, shapeList, calAList] = runMesoCellExtension(NV,NPINS,Kl,Kb,lambdaA,lambdaB,plotIt,movieFileStr);`

You can use values instead of variable names for the inputs. If you want to use variables, make sure they are defined in the MATLAB workspace. Remember that you do not need to include `movieFileStr`.

 
