# PREPARE_MalawiRisk
Data and MATLAB codes for seismic risk assessment in Malawi

The depository includes the following items: 
- Data_EQrisk_Malawi.mat
- open_code_simu_Malawi.m
- ModelFit_MalawiSeismicHazard.m.

The MATLAB data file (Item 1) contains information of the grid location information (latitude, longitude), hazard information (peak ground acceleration at the return periods of 100, 200, 500, 750, 1,000, 2,000, 2,500, 5,000, and 10,000 years), upper tail approximation results (probability distribution type, slope and intercept of the best-fitting model) and building collapse probability results. There are four sets of the building collapse probability results. The first set is for the equally weighted collapse fragility models in terms of ultimate behavior. For this result set, results for three census building classes (i.e., permanent, semi-permanent, and traditional) are included for the nine return periods. In other words, the first result set consists of 27 columns (3 building types times 9 return periods). The same information by considering individual ultimate behavior (i.e., strength degradation, geometrical instability, and limited ductility) is included as the second, third, and fourth result sets, respectively. In the MATLAB data file, the results for the 2,905 grids (0.1-degree spacing) and the results for the 18,714 enumeration areas are included. 

The MATLAB script file (Item 2) is the main MATLAB code: (i) to perform the upper tail approximations of seismic hazard curves, (ii) to simulate the annual maximum peak ground acceleration values for the specified duration, and (iii) to calculate the probability of collapse of a specified building type. The building type can be specified from nine cases (i.e., seismic vulnerability classes A, B, and C, and ultimate behavior of a static pushover curve, i.e., strength degradation, geometrical instability, and limited ductility). Building collapse risk simulations can be carried out for the 2,905 grids or for the 18,714 enumeration areas.

The MATLAB script file (Item 3) is the subfunction to carry out the upper tail approximation of a seismic hazard curve by fitting a different probability distribution.
