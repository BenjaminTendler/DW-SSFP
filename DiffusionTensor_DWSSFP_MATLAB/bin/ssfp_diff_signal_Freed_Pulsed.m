function [M] = ssfp_diff_signal_Freed_Pulsed(G,tau,TR,alpha,D,T1,T2,shortpulse)
% function [M] = ssfp_diff_signal_Freed(G,tau,TR,alpha[,D,T1,T2])
%
% calculates steady-state magnetization with diffusion weighting
% for monopolar diffusion gradient based on the Freed model. 
%
% input:
%	G   = strength of diffusion gradient (G/cm)
%	tau = duration of diffusion gradient (ms)
%	TR  = repetition time (ms)
%	alpha = flip angle (degrees)
%	D   = diffusion coefficient (mm^2/s)
%	T1  = longitudinal relaxation time (ms)
%	T2  = transverse relaxation time (ms)
%   shortpulse = shortpulse approximation (1 = yes, 0 = no)
%   output:
%	M = steady-state magnetization
%

  gam = 4258*2*pi;   % Hz/G

  TR = TR*10^-3;     % convert to s
  tau = tau*10^-3;   % convert to s
  T1 = T1*10^-3;     % convert to s
  T2 = T2*10^-3;     % convert to s
  G = G*10^-1;       % convert to G/mm
  alpha = alpha*pi/180;

%%
cosa = cos(alpha);
sina = sin(alpha);
E1p = @(p)(exp(-TR/T1-D*gam^2*G^2*tau^2*TR.*p.^2));
%For short pulse approximation
if shortpulse==1
    E2p = @(p)(exp(-TR/T2-D*gam^2*G^2*tau^2.*TR.*(p^2)));
else
    E2p = @(p)(exp(-TR/T2-D*gamma^2*G^2*tau^2.*((p^2+p+1/3)*(tau)+(p+1)^2*(TR-tau))));
end
Ap=@(p)(1/2.*(E1p(p)-1).*(1+cosa));
Bp=@(p)(1/2.*(E1p(p)+1).*(1-cosa));
Cp=@(p)(E1p(p)-cosa);
np=@(p)(-E2p(-p).*E2p(p-1).*(Ap(p).^2).*Bp(p-1)./Bp(p));
dp=@(p)((Ap(p)-Bp(p))+E2p(-p-1).*E2p(p).*Bp(p).*Cp(p+1)./Bp(p+1));
ep=@(p)(-E2p(p).*E2p(-p-1).*Bp(p).*Cp(p+1)./Bp(p+1));

x1=0;
for k=10:-1:1
    if k==10
        x1=np(k)./(dp(k)+ep(k));
    else
        x1=np(k)./(dp(k)+x1);
    end
end  
r1=x1./(E2p(-1).*Bp(0))+(E2p(0).*Cp(1))./Bp(1);

 M=r1.*sina.*(1-E1p(0)).*E2p(-1)./(Ap(0)-Bp(0)+E2p(-1).*Cp(0).*r1);
