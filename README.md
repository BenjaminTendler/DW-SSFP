# DW-SSFP
This repository contains code to simulate and process Diffusion-Weighted Steady-State Free Precession (DW-SSFP) data.

Details are as follows:

**DWSSFP_simulation_plot.py**

Code will simulate the expected DW-SSFP signal as a function of diffusion gradient amplitude/diffusion gradient time/TR/T1/T2/flip angle/diffusion coefficient. Simulations are based on the Buxton model of DW-SSFP (Buxton, Richard B. "*The diffusion sensitivity of fast steady‐state free precession imaging.*" Magnetic resonance in medicine 29.2 (1993): 235-243).

Requires *numpy*, *sys* and *matplotlib*

**DiffusionTensor_SSFP_2TP.py**

Code will fit DW-SSFP data to a diffusion tensor model, under the conditions of the two transverse-period approximation (defined in Eq. [3] of Buxton, Richard B. "*The diffusion sensitivity of fast steady‐state free precession imaging.*" Magnetic resonance in medicine 29.2 (1993): 235-243), valid when TR >= ~1.5 * T2. Deviations from this condition will lead to errors in estimated eigenvalues. 

This code provides rapid estimation of the tensor eigenvalue and eigenvectors, as ADCs are estimated using an analytical solution for the two transverse-period approximation (Appendix Eq. [A2] in Tendler, Benjamin C., et al. "*Modeling an equivalent b‐value in diffusion‐weighted steady‐state free precession.*" Magnetic resonance in medicine 84.2 (2020): 873-884.)

Requires *numpy*, *sys*, and *nibabel*

**Ball2Sticks_DWSSFP** and **DiffusionTensor_DWSSFP**

Code to fit a ball & 2 sticks model and diffusion tensor model to DW-SSFP data aquired at one and two flip angles, based on the process described in Tendler, Benjamin C., et al. "*Use of multi-flip angle measurements to account for transmit inhomogeneity and non-Gaussian diffusion in DW-SSFP*", NeuroImage 220 (2020): 117113. Fitting performed using the cuDIMOT toolbox (https://users.fmrib.ox.ac.uk/~moisesf/cudimot/). 

**GammaFitting**

Code to fit a gamma distribution of diffusivisities to multi-flip angle DW-SSFP data, as described in Benjamin C., et al. "*Use of multi-flip angle measurements to account for transmit inhomogeneity and non-Gaussian diffusion in DW-SSFP*", NeuroImage 220 (2020): 117113. Current version is implemented with signal weighting based on the standard deviation of the diffusivity estimates at each flip angle. 

Any questions please contact benjamin.tendler@ndcn.ox.ac.uk

