function [fit_out,Sig] = TensorFittingWrapAnalytical(noisefloor,bvecs,tau,G,TR,alpha,T1,T2,Data)
%%
%Estimate S0 from b0 data
DataS0=squeeze(abs(Data(sum(abs(bvecs),1)==0).^2-noisefloor(sum(abs(bvecs),1)==0)'.^2).^0.5);
S0=abs(mean(DataS0)./ssfp_diff_signal_Freed_Pulsed(0,0,TR(1)*1000,alpha(1),0,T1,T2,0));
if isnan(S0) || isinf(S0)
    fit_out=zeros([1,6]);
    Sig=zeros(size(Data'));
else
    %%
    %Define fitting options
    options = optimoptions(@lsqnonlin,'Display','off');
    %%
    %Perform fitting
    f=@(x)TensorFittingAnalytical(x,Data,S0,alpha,TR,tau,T1,T2,G,bvecs,noisefloor,1);
    [fit_out]=lsqnonlin(f,[0.0001,0,0,0.0001,0,0.0001],[0,-0.01,-0.01,0,-0.01,0],[0.01,0.01,0.01,0.01,0.01,0.01],options);
    Sig=TensorFittingAnalytical(fit_out,Data,S0,alpha,TR,tau,T1,T2,G,bvecs,noisefloor,0);
end

