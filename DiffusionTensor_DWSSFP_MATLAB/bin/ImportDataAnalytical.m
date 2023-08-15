function [Data,T1,T2,B1,mask,noisefloor,bvecs,FlipAngleVector,tau,G,TRs] = ImportDataAnalytical(DataPath,T1Path,T2Path,B1Path,maskPath,noisefloorPath,bvecsPath,FlipAngleVectorPath,tauPath,GPath,TRsPath)
%%
%Import data
Data=niftiread(DataPath);
T1=niftiread(T1Path);
T2=niftiread(T2Path);
B1=niftiread(B1Path);
mask=niftiread(maskPath);
noisefloor=importdata(noisefloorPath);
bvecs=importdata(bvecsPath);
FlipAngleVector=importdata(FlipAngleVectorPath);
tau=importdata(tauPath);
G=importdata(GPath);
TRs=importdata(TRsPath);
%%
%Distinguish b0 and diffusion-weighted volumes
bvecs(:,tau==0)=0;
