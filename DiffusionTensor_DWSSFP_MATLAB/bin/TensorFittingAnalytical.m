function [S,J] = TensorFittingAnalytical(x,Data,S0,alpha,TR,tau,T1,T2,G,bvecs,noisefloor,res)
%%
%Define parameters - eigenvalues and vector angles
Dxx=x(1);
Dxy=x(2);
Dxz=x(3);
Dyy=x(4);
Dyz=x(5);
Dzz=x(6);
%%
%%
%Define ADC
ADC=bvecs(1,:).*bvecs(1,:)*Dxx+bvecs(2,:).*bvecs(2,:)*Dyy+bvecs(3,:).*bvecs(3,:)*Dzz+2*bvecs(1,:).*bvecs(2,:)*Dxy+2*bvecs(1,:).*bvecs(3,:)*Dxz+2*bvecs(2,:).*bvecs(3,:)*Dyz;
%%
%Estimate signal
Sig=zeros(size(ADC));
for k=1:length(ADC)
    Sig(k)=S0.*ssfp_diff_signal_Freed_Pulsed(G(k),tau(k)*1000,TR(k)*1000,alpha(k),ADC(k),T1,T2,0);
end
%%
%Add noise floor
S=(Sig.^2+noisefloor.^2).^0.5;
%%
%Calculate residual
if res==1
    S=double((abs(S)'-abs(Data)));
end
