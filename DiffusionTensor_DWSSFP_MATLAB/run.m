%%
%Define directory path
FilePath='ExampleData/';
OutputPath='Output/';
%%
%Define data paths
[DataPath,T1Path,T2Path,B1Path,MaskPath,bvecsPath,FlipAnglePath,tauPath,GPath,TRsPath,noisefloorPath]=DatasetPaths(FilePath);
%%
%Import files
[Data,T1,T2,B1,Mask,noisefloor,bvecs,FlipAngle,tau,G,TRs] = ImportDataAnalytical(DataPath,T1Path,T2Path,B1Path,MaskPath,noisefloorPath,bvecsPath,FlipAnglePath,tauPath,GPath,TRsPath);
%%
%Run Fitting
[D,SigReconTensor]=TensorAnalyticalWrap(Data,T1,T2,B1,Mask,FlipAngle,noisefloor,tau,G,TRs,bvecs);
%Output Tensor
TensorOut(D,SigReconTensor,FilePath,OutputPath,'');