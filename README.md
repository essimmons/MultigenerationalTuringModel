# A multigenerational Turing model reproduces transgressive petal spot phenotypes in hybrid *Mimulus*: Supplementary Material

Emily S. G. Simmons, Arielle M. Cooley, Joshua R. Puzey, and Gregory D. Conradi Smith

## Description
The code included in this project reproduces simulations for a population of F2 hybrids as submitted to bioRxiv at https://doi.org/10.1101/2023.08.10.552669. 

### Included files
- Main file for running simulations.
	- MultigenTuringModel.m
- Supporting function files (listed alphabetically).
	- buildF2Fig.m
    - config.m   
	- findCriticalk2.m
	- findJacobian.m
	- findSteadyState.m
	- implicitRxnDiff2D.m
	- interpretZygosity.m
	- revertFromF1ToGenes.m
	- runF2Simulations.m
- 6 different parameter files titled "FigXXParams.mat" which reproduce the corresponding published figures.
	- Fig04Params.mat
	- Fig5AParams.mat
	- Fig5BParams.mat
	- Fig5CParams.mat
	- Fig5DParams.mat
	- Fig6Params.mat

## Requirements
This project is built on MATLAB_2023a.

## Usage
- Open MatLab
- Load parameters from one of the 6 corresponding parameter files, or make your own.
	- If making your own parameter set, reference one of the parameter sets provided for formatting.
- Set up your simulation configurations in `config.m`
	- `configs.N`: the number of mesh points, $N$, you wish to use in your simulation.
	    - To reproduce the figures as they appear in the paper, choose `configs.N=100`.
	- `configs.axisSize`: Choose the spatial domain as multiples of the anticipated F1 wavelength, $\omega$.
	- `configs.makeF2Fig`: Choose whether or not you want to build and save the figure. 
		- If yes, set `configs.makeF2Fig=1`. 
		- If no, set `configs.makeF2Fig=0`.
	- `configs.showProgress`: Choose whether or not you want to see the progression of the simulations as they are calculated.
		- If yes, set `configs.showProgress=1`.
		- If no, set `configs.showProgress=0`.
- Run the file `MultigenTuringModel.m`.
- When complete, the program will save the steady-state activator concentration results of all 12 simulations as
		"UdataSz"+ `N` + ".mat"
		"UdataSz"+ `N` + ".dat"
	- If you opted to build the figure, the figure will be saved as 
		  "F2SimulationSz" +  `N`  + ".fig"
