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

`function [hList, xList, yList, shapeList, calAList] = runMesoCellExtension(NV,NPINS,calA0Init,Kl,KbInit,cL,cB,cKb,plotIt,movieFileStr)`

### Inputs:
* `NV`: integer number of vertices that make up the deformable particle (DP)
* `NPINS`: integer number of pins to force the DP to change shape during extension 
	* **NOTE**: MUST be at least _less than_ `NV`
* `calA0Init`: initial **preferred** shape parameter, will change if you make `cL` greater than 0. 
	* Note that `calA0Init` should be greater than or equal to 1 always. 
* `Kl`: mechanical constant for perimeter
* `KbInit`: mechanical constant for curvature
	* **NOTE**: `KbInit` will change if parameter `cKb > 0`. 
* `cL`: rate constant for growth/decay of the preferred shape parameter `calA0` based on the instantaneous shape parameter `calA`
* `cB`: rate constant for growth of preferred angles `th0` during stretching simulation
* `cKb`: rate constant for change in bending energy, currently updates based on magnitude of angle
* `plotIt`: Binary var. to either draw (`1`) or not draw (`0`) cells during simulation

#### Optional Input:
* `movieFileStr`: If included, save animation to movie file. 
	* **NOTE**: File string must end with the `.mp4` extension, currently does not support other formats.

### Outputs:
* `hList`: List of pin stretch steps parametrized by variable <img src="https://render.githubusercontent.com/render/math?math=h">. 
	* **NOTE**: h is unitless, defined as the distance travelled outside of the original perimeter divided by the vertex radius.
* `xList`: _x_-coordinates of all cell vertices at each stretch step
* `yList`: _y_-coordinates of all cell vertices at each stretch step
	* **NOTE**: Coordinates are stored as MATLAB `cell` variables. To access at step `ii`, use syntax `x{ii}` or `y{ii}`.
* `shapeList`: List of relevant shape information during simulation. Column 1 is `calA0` (_double_), column 2 is `Kb` (_array_), column 3 is preferred angles `th0` (_array_)
* `calAList`: Instantaneous shape parameter <img src="https://render.githubusercontent.com/render/math?math=\mathcal{A} = p^2/4\pi a"> during stretching simulation.

## Running the code

To run the code on the command line, simply use

`>> [hList, xList, yList, shapeList, calAList] = runMesoCellExtension(NV,NPINS,calA0Init,Kl,KbInit,cL,cB,cKb,plotIt,movieFileStr);`

You can use values instead of variable names for the inputs. If you want to use variables, make sure they are defined in the MATLAB workspace. Remember that you do not need to include `movieFileStr`.

End of README.md
