# DW-SSFP
This repository contains code to simulate the Diffusion-Weighted Steady-State Free Precession (DW-SSSFP) signal. Details are as follows:

**DWSSFP_simulation_plot.py**

Code will simulate the expected DW-SSFP signal as a function of diffusion gradient amplitude/diffusion gradient time/TR/T1/T2/flip angle/diffusion coefficient. Simulations are based on the Buxon model of DW-SSFP (Buxton, Richard B. "*The diffusion sensitivity of fast steady‐state free precession imaging.*" Magnetic resonance in medicine 29.2 (1993): 235-243).

Requires *numpy*, *sys* and *matplotlib*

**DiffusionTensor_SSFP_2TP.py**

Code will fit DW-SSFP data to a diffusion tensor model, under the conditions of the two transverse-period approximation (defined in Eq. [3] of Buxton, Richard B. "*The diffusion sensitivity of fast steady‐state free precession imaging.*" Magnetic resonance in medicine 29.2 (1993): 235-243), valid when TR >= ~1.5 * T2. Deviations from this condtion will lead to errors in estimated eigenvalues. 

This code provides rapid estimation of the tensor eigenvalue and eigenvectors, as ADCs are estimated using an analytical solution for the two transverse-period approximation (Appendix Eq. [A2] in Tendler, Benjamin C., et al. "*Modeling an equivalent b‐value in diffusion‐weighted steady‐state free precession.*" Magnetic resonance in medicine 84.2 (2020): 873-884.

Requires *numpy*, *sys*, and *nibabel*

Any questions regarding this code and DW-SSFP, please contact benjamin.tendler@ndcn.ox.ac.uk

