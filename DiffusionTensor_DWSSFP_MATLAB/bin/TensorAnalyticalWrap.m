function [D,SigRecon] = TensorAnalyticalWrap(Data,T1,T2,B1,mask,FlipAngleVector,noisefloor,tau,G,TR,bvecs,sliceNo)
%%
%Create Output Arrays
D=zeros([size(Data,1),size(Data,2),size(Data,3),6]);
SigRecon=zeros(size(Data));
%%
%Define slice to fit
if nargin==12
    loopInit=sliceNo;
    loopEnd=sliceNo;
else
    loopInit=1;
    loopEnd=size(Data,3);
end
%%
%Create Flip Angle file
for k=1:length(FlipAngleVector)
    FlipAngle(:,:,:,k)=round(B1.*FlipAngleVector(k)).*mask;
end
%%
%Perform Fitting
if sum(T1(:))==0
    ['Slice Complete']
else
    for k=loopInit:loopEnd
        %Slice data to reduce dimensionality if parallelising
        [DataSlice,T1Slice,T2Slice,FlipAngleSlice,row,col] = pullSlice(Data,T1,T2,FlipAngle,k);
        %Perform fitting
        %parfor l=1:length(row)
        for l=1:length(row)
            [fit_out(:,l),Sig(:,l)] = TensorFittingWrapAnalytical(noisefloor,bvecs,tau,G,TR,squeeze(FlipAngleSlice(row(l),col(l),:)),squeeze(T1Slice(row(l),col(l))),squeeze(T2Slice(row(l),col(l))),squeeze(DataSlice(row(l),col(l),:)));
        end
        %Assign Variables
        for l=1:length(row)
            D(row(l),col(l),k,:)=fit_out(1:6,l);
            SigRecon(row(l),col(l),k,:)=Sig(:,l);
        end
        clear fit_out Sig
        ['Slice ',num2str(k),' of ',num2str(size(Data,3)), ' Complete']
    end
end
