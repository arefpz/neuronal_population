# Neuronal_population
The code runs a population of N LIF neurons.

# Developers
The code was developed by Aref Pariz and Farhad Daei in 2021.

# Code description
This repository includes several codes written in FORTRAN to simulate a population(s) of neurons.  
The code can be run either on a single core or multiple cores using the OpenMP flag. 

Specifically, these codes were used to generate the results for the following paper:

Selective control of synaptic plasticity in heterogeneous networks through transcranial alternating current stimulation (tACS)  
Aref Pariz, Daniel Trotter, Axel Hutt, Jeremie Lefebvre  
https://www.biorxiv.org/content/10.1101/2022.11.15.516556v1

Although here we have used the LIF neuron model with Hebbian Spike Timing dependant plasticity (STDP), it is possible to change the code to work with any neuron model by adding it in modu.f90 file. 

For more information about the code, please read the Materials and Methods section in the aforementioned paper.

Further information and a video explaining the code will be available ASAP.

# Description of the files
vars.f90: the parameters and variables used in the code.  
modu.f90: the code contains modules that are being used to generate the connectivity, synaptic weights, delay matrices, and membrane potential evolution.  
writer.f90: The file contains the subroutines that can be used to write the data as a .txt file on disk.  
log.f90: The code uses vars.f90 to write the used parameters in simulation onto the disk.  
random.f90: The publicly available code is written by Alan Miller. For more information, please read the description inside the code.  
main_p.f90: this is the main file that simulates the activity of neurons.  

# Running the code
The code needs to be compiled first using gfortran. The following command will compile and run the code.  

- on a single core:  
gfortran -m64 -O3 vars.f90 log.f90 writer.f90 random.f90 modu.f90 main_p.f90 -o a.out  #This compiles the code 
./a.out

- on multiple cores:  
gfortran -m64 -O3 -fopenmp vars.f90 log.f90 writer.f90 random.f90 modu.f90 main_p.f90 -o a.out  
./a.out

In the next step, you must provide the stimulation frequency (variable: omega) and amplitude (variable: amp). The amplitude will then be multiplied by the damp value in vars.f90.
If you are simulating a population without stimulation, set the variable IS_SIgnal to 0 and comment the line in main_p.f90, which asks you to enter the values.  
Note that to run the code on several cores, you must define the number of CPUs in "ompnum" variable in vars.f90. Default value is 4.

# Converting the results into .MAT file
Use MATLAB software to read and save the data as MAT files. Use the Matlab codes (*.m) to analyze and plot the figures.  
data_reader.m: By running this code, the code will try to load the .txt files and save them as a .MAT files with appropriate names for each part of simulations and ensemble.  
analyser_taum_less_more.m: The code uses the data saved by "data_reader.m" and moved to "data/" folder. This code finds the distribution of the synaptic weights among neurons with different membrane time constants and saves the result.  
g_plotter.m: the code uses the result from the previous step and plots the distribution of synaptic weights for different stimulation frequencies.  


# Coupled neurons
stdp_coupled_LIF_neuron.m simulates the coupled neurons stimulated by tACS and save the results in 'data/.'

analyzer.m will analyze the data saved in the 'data/' folder. The outcome contains the average synaptic weight between coupled neurons at different stimulation frequencies and membrane time constant of second neurons, while the MTC of neuron 1 is 10ms.

dT_spikes.m will find the distribution of spike timing differences among neuron spikes.

phase_diagram_taum.m finds the phase of spikes for both neurons.

stdp_sum_omega.m finds the convolution of the STDP time window and the spike timing differences of neurons.
