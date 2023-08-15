function [DataSlice,T1Slice,T2Slice,FlipAngleSlice,row,col] = pullSlice(Data,T1,T2,FlipAngle,sliceNo,mask)
%%
%Extract slice from data
DataSlice=(Data(:,:,sliceNo,:));
T1Slice=(T1(:,:,sliceNo));
T2Slice=(T2(:,:,sliceNo));
FlipAngleSlice=(FlipAngle(:,:,sliceNo,:));
%%
%Obtain row and column indices
if nargin==5
    [row,col]=find(T1Slice);
else
    [row,col]=find(T1Slice.*mask(:,:,sliceNo));
end