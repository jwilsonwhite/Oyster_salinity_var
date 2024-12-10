function kxy = mkkern(x,y,Sal,Param,M_master,Type)

% Construct kernels for IPM
% NOTE: at present there is no salinity effect on growth/mortality, so
% input argument 'Sal' goes unused here

isjuv = 1 - normcdf(x,Param.Lf,Param.Lfs); % which sizes get fished
isjuv2 = 1 - normcdf(x,Param.j_length,1);  % which sizes experience juvnile mortality rate

%SURVIVAL PART OF KERNEL
switch Type
    case {'mortality','growth','both','fecundity'}
  % natural mortality only
Ma = Param.M; 
Mj = Param.Mj; 
m = ones(size(x)).*Ma.*(1-isjuv2) + ones(size(x)).*Mj.*(isjuv2);
    case {'fishing'}
M = M_master;
m = (1-isjuv).*M; 
end
%this is a matrix size x,
%mortality for each size 
                                              
p1 = exp(-m); % convert mortality rate to survivorship

%GROWTH PART OF KERNEL
Linf = Param.Linf;
k = Param.k; 
%growth
pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)
pmean1 = min(pmean1,Linf);

%add variability around von Bertalanffy growth
psig1 = pmean1*Param.vBs; % make this a parameter; %.*varparm(3); 
psig1 = min(psig1,Param.vBs.*Linf); % cap the variability to prevent leakage into very large sizes

%evaluate growth part of kernel
p2 = normpdf(y, pmean1, psig1);
%keyboard
% FECUNDITY PART OF KERNEL

%define pr(reproductive) using a maturity ogive
ismat = normcdf(x,Param.Mat,diff(x(1,1:2))/2);
pmean2 = repmat(Param.Fec(:)',[length(x),1]); 
pmean2 = ismat.*pmean2; 

pmean3 = repmat(Param.Rvec(:),[1,length(x)]);
p3 = pmean2 .* pmean3;

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives



%added LSS - enforcing that each col of the kernel sums to 1/dx
for i = 1:length(p2)
    p2(:,i) = (p2(:,i)/sum(p2(:,i)))*(1/diff(x(1,1:2)));
end


if k == 0 %when k = 0 the kernel goes "funny" so making a special case zero growth kernel where the chances of staying in your current size class are 100%
    p2 = eye(size(p2))*(1/diff(x(1,1:2)));
end
%end added LSS     


switch Type
    case 'growth'
        kxy = p2;
    case {'mortality','fishing'}
        kxy = diag(diag(p1));
    case 'both'
kxy = p1.*p2;
    case 'fecundity'
        kxy = p3;
end

assignin('base',"p1",p1)
assignin('base',"p2",p2)
assignin('base',"p3",p3)


end