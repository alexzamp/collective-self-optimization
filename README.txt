This is a repository containing the relevant simulation code and data for our manuscript entitled: 
"Collective self-optimization of communicating active particles",
by A. Zampetaki, B. Liebchen, A. V. Ivlev and H. Löwen.

This code is  is made available under the GPL license. Over and above the legal restrictions imposed by this license, if you use any of the contents for academic research then you are obliged to provide proper attribution to the related paper:

A. V. Zampetaki, B. Liebchen, A. V. Ivlev and H. Löwen, 
"Collective self-optimization of communicating active particles", 
(publication details to be provided soon).

## Codes
Two codes are provided "ham1.c" and "non_ham1.c", which were used to extract the simulation data for the Hamiltonian approximation and the full non-Hamiltonian system,studied and discussed in detail in the aforementioned paper. The codes were written for running in Windows 10, but can be easily adapted to run in Linux by changing appropriately the random number generator functions.

## Data
The enclosed raw data refer to the Figures 3 and 5 of the manuscript and include the final configurations used for the phase diagram of Fig.3 and the 14 particle trajectories used to construct Fig. 5. Due to size limitations they are provided in a compressed form (".zip" files). The ".zip" files for Figure 3 contain folders with names such as "Ph1Dk1.000_T0.010_Lx24_Np256". It is thus implied that k=1.000, T=0.010, Lx=24 and Np=256. Here k is the screening paramete, T the optimal value of the field, Lx the x-dimension of the simulation box and Np the number of particles. Some folder names include also the value of the diffusion contant Dt used in the simulation.  In turn, the folders contain ".dat" files whose names correspond to  the different  values of the time step.  Each ".dat" file has 5 space-separated columns  representing:
1. time step,
2. x-coordinate,
3. y-coordinate,
4. total force projection in x,
5. total force projection in y.
Each row then accounts for a different particle. More information can be provided in the relevant code "non_ham1.c" and in the paper.

## Contact
Any questions regarding the enclosed code and data should be addressed to:
alexandra.zampetaki@hotmail.com
