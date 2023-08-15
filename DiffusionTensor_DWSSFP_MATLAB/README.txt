Software estimates the Diffusion Tensor from experimental DW-SSFP data. 
Implementation based on the DW-SSFP model described in Freed et al., Steady-state free precession experiments and exact treatment of diffusion in a uniform gradient,J. Chem. Phys. 115, 4249 (2001); https://doi.org/10.1063/1.1389859

INPUT DATA:

The input data directory (FilePath) should contain the following files
- data.nii.gz (b0 & dwi data concatenated into a single file)
- T1map (T1 map coregistered to diffusion data)
- T2map (T2 map coregistered to diffusion data)
- B1map (B1 map coregistered to diffusion data)
- mask  (mask coregistered to diffusion data)
- bvecs (bvectors text file with dimensions 3xN)
- flipAngles (flip angle text file in degrees with dimensions 1xN)
- diffGradDurs (diffusion gradient duration file in seconds with dimensions 1xN)
- diffGradAmps (diffusion gradient amplitude file in G/cm with dimensions 1xN. Divide mT/m by 10 to arrive at G/cm. SET AMPLITUDE EQUAL TO ZERO FOR b0 VOLUMES)
- TRs (TR file in seconds with dimensions 1xN)
- noisefloor (noise floor estimate file in with dimensions 1xN. Create 1xN files of all zeros if no noise floor estimate)

OUTPUT DATA:

FSL is used to convert the estimated Tensor components into Tensor outputs (V1/2/3, L1/2/3, FA, MD, MO). Please set the paths in "setenvFSL.m" for your machine, or edit the output file (TensorOut) for your own conversion software.

NOTES:
The DW-SSFP data must be fit voxelwise to estimate Tensors, leading to lengthy processing times. This code can be readily parallelised on a computing cluster. For MATLAB parallelisation, uncomment the 'parfor' loop in TensorAnalyticalWrap.m 

