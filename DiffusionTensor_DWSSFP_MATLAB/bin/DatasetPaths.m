function [DataPath,T1Path,T2Path,B1Path,MaskPath,bvecsPath,FlipAnglePath,tauPath,GPath,TRsPath,noisefloorPath]=DatasetPaths(FilePath)
%%
%Define Paths
DataPath=[FilePath,'data.nii.gz'];
T1Path=[FilePath,'T1map.nii.gz'];
T2Path=[FilePath,'T2map.nii.gz'];
B1Path=[FilePath,'B1map.nii.gz'];
MaskPath=[FilePath,'mask.nii.gz'];
bvecsPath=[FilePath,'bvecs'];
FlipAnglePath=[FilePath,'flipAngles'];
tauPath=[FilePath,'diffGradDurs'];
GPath=[FilePath,'diffGradAmps'];
TRsPath=[FilePath,'TRs'];
noisefloorPath=[FilePath,'noisefloor'];
