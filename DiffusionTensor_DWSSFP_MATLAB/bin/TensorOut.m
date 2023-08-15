function TensorOut(D,SigRecon,FilePath,OutputPath,PathAppend)
%%
%Initialisation message
disp('Outputting Tensor data...');
%Set FSL Environment
setenvFSL();
%%
%Add slashes
FilePath=[FilePath,'/'];
OutputPath=[OutputPath,'/'];
PathAppend=[PathAppend,'/'];
%%
system(['mkdir -p ',OutputPath,PathAppend]);
%%
%Output Data
niftiwrite(D,[OutputPath,PathAppend,'dti'],'Compressed',true);
niftiwrite(SigRecon,[OutputPath,PathAppend,'dti_SigRecon'],'Compressed',true)
%%
%Correct header information
system([getenv('FSLDIR'),'bin/fslcpgeom ',FilePath,'T1map ',OutputPath,PathAppend,'/dti ', '-d']);
system([getenv('FSLDIR'),'bin/fslcpgeom ',FilePath,'T1map ',OutputPath,PathAppend,'/dti_SigRecon ', '-d']);
%Create diffusion tensor components
system([getenv('FSLDIR'),'bin/fslmaths ',OutputPath,PathAppend,'dti',' -tensor_decomp ',OutputPath,PathAppend, 'dti']);
%%
%Completion message
disp('Complete');
