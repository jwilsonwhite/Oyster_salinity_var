%edited by Laura July 2022 - using new mock salinity time series, salinity
%is an input now 

%function Oyster_IPM_Disturbance(Meta,T,TS_temp,Rn,F,Fp,Y,Runs,savename)
function Oyster_IPM_Disturbance(Meta,T,TS_temp,TS_sal,Rn,F,Fp,savename,saveopt,salopt)

%for i = 1:Runs % trial run: 10 simulations; full run: 10000 sims

warning off all

if ~exist('saveopt','var')
    saveopt = false;
end

if ~exist('salopt','var')
    salopt = 'none';
end

%salsave_tmp = strcat('Salinity_ts/',salsave)
%load(salsave_tmp)

% Implement Oyster Predator-Prey IPM w/ Varying Disturbance Regimes

% Inputs:
% Meta = structure of parameters & metadata from 'oyster_PP_params.mat'
% T = length of simulation [T = 40 semi-years]
% TS_temp = temperature timeseries. 26 weeks * T semiyears [T = 40]
% Rn = recruitment; used for 'Open' and 'Mixed' recruitment only
% F = oyster fishing mortality (e.g., 0.1)
% Fp = drill fishing/culling mort.
% savename = save as filename (e.g., 'y0_2000sims.mat')

% Note that predation & other parameters is on a weekly time scale
% Spatial units = 0.25m^-2

% TS_sal = salinity timeseries. 26 weeks * T semiyears [T = 40];
% y# (e.g., y0) = return time factor scaled from 0-3 by 0.5
% see ClimatologyForDisturbanceOysterIPM.m for details on generation of
% salinity data. "i" corresponds to i-th time series for i-th simulation

% 2) initialize model starting conditions
N = nan(Meta.IPM.Prey.meshsize,T);
P = nan(Meta.IPM.Pred.meshsize,T);
%P_master = Meta.Params.Pred.Sizedist;
N(:,1) = 1;
P(:,1) = 5/Meta.IPM.Pred.meshsize/Meta.IPM.Pred.dy; % starting density of ~5 per 0.25 m^-2

% other state variables:
S = nan(T,1); % reef structure
R = nan(T,1); % number of recruits
S(1) = 1;
R(1) = 1;

if saveopt
Ntss = nan(Meta.IPM.Prey.meshsize,26,T);
end
RunningMean = TS_sal(1,1); % initial value (1st salinity value)

for t = 2:T % each timestep is 6 months. 
    % ****** For data comparison model will have to deal with date-alignment
    % issues.
    
    % Beginning of timestep: oyster & drill reproduction & recruitment
    A = S(t-1);% + Meta.IPM.Prey.Sy*(Meta.Params.Prey.LW(:).*N(:,t-1)); % total area available for settlement in a timestep
    A = A.^2/3; % convert biomass to area
    switch Meta.Rect
        case 'Open'
    RR = Rn(t);  % external recruitment
    RR = RR.*Meta.Params.Prey.Rvec*Meta.IPM.Prey.dy;
        case 'Closed' % USE 'Closed' Recruitment -> Meta.Rect = 'Closed'
   kmatNr = kernmatSimp(Meta.IPM.Prey.x,TS_sal(1,t),Meta.Params.Prey,0,'fecundity'); % get kernel        
   RR = kmatNr*N(:,t-1)*Meta.IPM.Prey.dy;         
        case 'Mixed'
    kmatNr = kernmatSimp(Meta.IPM.Prey.x,TS_sal(1,t),Meta.Params.Prey,0,'fecundity'); % get kernel        
   RR = kmatNr*N(:,t-1)*Meta.IPM.Prey.dy + Rn(t)*Meta.IPM.Prey.dy; %sum(kmatNr*N(:,t-1)*Meta.IPM.Prey.dy + Rn(t));         
    end
    
% Apply salinity-dependent larval mortality:
    if t == 1 %currently this is never triggered because t starts at 2
        Savg = TS_sal(1,1); %-M_rs2;
    else
       % Savg = nanmean(TS_sal(23:end,t-1)); %mean of prior 4 weeks
       % Tavg = nanmean(TS_temp(23:end,t-1)); %mean of prior 4 weeks
       Savg = (TS_sal(23:end,t-1)); %mean of prior 4 weeks
       Tavg = (TS_temp(23:end,t-1)); %mean of prior 4 weeks
    end
    
    switch salopt
        case 'NoLarval'
            M_r = 0.76; % approximate maximum of the penalty function 

            disp(M_r)
        otherwise
% New approach: combination of early & later survival
% 1st week
 M_r1 = Meta.Params.Prey.b0 + Meta.Params.Prey.b1.*Tavg(1) + Meta.Params.Prey.b2.*Savg(1) + Meta.Params.Prey.b3.*Tavg(1).^2 + Meta.Params.Prey.b4.*Savg(1).^2 + Meta.Params.Prey.b5.*(Tavg(1).*Savg(1));
 M_r1 = max(0,M_r1);
 %M_r2 = -104.389 + 10.221.*Tavg(2:end) + -0.2006.*Tavg(2:end).^2 + 0.0996 .* Tavg(2:end) .* Savg(2:end) + -0.1455 .* Savg(2:end).^2 + 2.6621 .* Savg(2:end);
 M_r2 = Meta.Params.Prey.b0a + Meta.Params.Prey.b1a.*Tavg(2:end) + Meta.Params.Prey.b2a.*Tavg(2:end).^2 + Meta.Params.Prey.b3a.*Tavg(2:end).*Savg(2:end) + Meta.Params.Prey.b4a.*Savg(2:end).^2 + Meta.Params.Prey.b4a.*Savg(2:end)
 M_r2 = max(0,M_r2);
   % M_r = prod(M_r).^(1/4); % geometric mean...this caused issues with
   % zeros so use arithmetic instead
    M_r2 = mean(M_r2);

    M_r1 = M_r1/100; %it's a % survival but the number is on O(10^1) instead of decimal
    M_r2 = M_r2/100;
    M_r = M_r1.*M_r2.^3;  
    end % end switch salopt


    RR = RR * M_r * exp(-Meta.Params.Prey.M_r); 
    RR_surv = 0.7./(1 + sum(RR).*4.*6e-4);
    RR = RR * RR_surv;
  

    N(:,t) = N(:,t-1) + RR(:);

    
   
   kmatPr = kernmatSimp(Meta.IPM.Pred.x,TS_sal(1,t),Meta.Params.Pred,0,'fecundity'); % pred. fec. kernel        
   RR_P = kmatPr*P(:,t-1)*Meta.IPM.Pred.dy;
   
   LEP = Meta.Params.Pred.LEP; % drill lifetime egg production (R0)
   a = 1/(0.1*LEP); % a = slope at the origin of the Beverton Holt curve
   %b = 15; % max density drills = 15/m^2 (from Kimbro 2017 paper)
   b = 0.25;
   RR_P = (a*RR_P)./(1+(a*RR_P)/b); % pred. DD, Beverton Holt like equation
   P(:,t) = P(:,t-1) + RR_P(:);
    
    % temporary state variables within this loop
    Nt = N(:,t); 
    Pt = P(:,t);
    St = S(t-1);
    
    for tt = 1:26 % sub-loop on a weekly timescale for growth & predation dynamics

    
    switch salopt
        case 'NoAdult'
            kmatNg = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey,0,'growth'); % use regular growth kernel
            kmatNm = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey,0,'mortality');  % use regular natural mortality kernel
   
        otherwise
    % Get salinity-dependent kernels:
    % oyster salinity-dependent conditions based on Kennedy et al. 1996 (Ch. 13,
    % pgs. 468-469). Additional citations: growth- La Peyre, Eberline, Soniat, & La Peyre, 2013; 
    % Lowe, Sehlinger, Soniat, & La Peyre, 2017. mortality- Schlesselman, 1955; Turner, 2006; 
    % Volety et al., 2009; La Peyre et al., 2013
    
    if TS_sal(tt,t) <= 5 || TS_sal(tt,t) >=30
  Meta.Params.Prey2 = Meta.Params.Prey; % new params = old params
  Meta.Params.Prey2.k = 0; % but no growth
  Meta.Params.Prey2.M = 0.0104; % higher mortality M (100% increase)
  Meta.Params.Prey2.Mj = 0.0372; % higher juv. mort. Mj (100% increase)
  kmatNg = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey2,0,'growth'); % use growth kernel with new params
  kmatNm = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey2,0,'mortality'); % use natural mortality kernel with new params
    else
  kmatNg = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey,0,'growth'); % use regular growth kernel
  kmatNm = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey,0,'mortality');  % use regular natural mortality kernel
    end
    end % end salopt
    
  kmatNf = kernmatSimp(Meta.IPM.Prey.x,TS_sal(tt,t),Meta.Params.Prey,F(t),'fishing'); % oyster harvest
  
  kmatPg = kernmatSimp(Meta.IPM.Pred.x,TS_sal(tt,t),Meta.Params.Pred,0,'growth'); % use regular growth kernel
  
  % drill salinity-dependent condition based on Butler 1985 (pg. 5)
  % 10ppt reasonable/consensus cutoff for higher mortality
  
  switch salopt
      case 'NoDrill'
                  kmatPm = kernmatSimp(Meta.IPM.Pred.x,TS_sal(tt,t),Meta.Params.Pred,0,'mortality');  % use regular natural mortality kernel
      otherwise
    if TS_sal(tt,t) <= 10
        Meta.Params.Pred2 = Meta.Params.Pred; % new params = old params
        Meta.Params.Pred2.M = 0.005; % but higher mortality M (100% increase)
        Meta.Params.Pred2.Mj = 0.005; % and higher juv. mort. Mj (100% increase)
        kmatPm = kernmatSimp(Meta.IPM.Pred.x,TS_sal(tt,t),Meta.Params.Pred2,0,'mortality'); % use natural mortality kernel with new params
    else
        kmatPm = kernmatSimp(Meta.IPM.Pred.x,TS_sal(tt,t),Meta.Params.Pred,0,'mortality');  % use regular natural mortality kernel
    end
  end % end switch salopt

    kmatPf = kernmatSimp(Meta.IPM.Pred.x,TS_sal(tt,t),Meta.Params.Pred,Fp(t),'fishing'); % drill harvest/culling
    
  %  Nt2 = (Meta.IPM.Prey.Symat.*kmatNm)*Nt; % oysters dead due to natural causes
    Nt2 = kmatNm*Nt;
    DeadM = Nt - Nt2;
    
   % Nt = (Meta.IPM.Prey.Symat.*kmatNg.*kmatNm.*kmatNf)*Nt; 
    Nt = (kmatNm.*kmatNf)*Nt;
    Nt = kmatNg*Nt.*Meta.IPM.Prey.dy; 
    % Oysters grow & experience natural mortality & fishing.a
    % Natural mortality leaves shell behind; fishing does not
    
    % Predators may change behavior or die due to salinity
    RunningMean = nanmean([RunningMean,TS_sal(tt,t)]); % update running mean salinity


    switch salopt
        case 'NoDrill' % just apply the temperature penalty
            if TS_temp(tt,t)<20
                salP = max(0,(TS_temp(tt,t)-15)/20);
            else
            salP = 1;
            end
        otherwise
    salP = predator_salinity_penalty(TS_sal(tt:end,t),TS_temp(tt,t),RunningMean);
    end

    
    % Predation depends on salinity, habitat complexity, prey density
   
    aP = Meta.Params.Pred.aP ; % attack rate
    
    hP = Meta.Params.Pred.hP ; % handling time
    % hP is a nxp matrix giving the handling time for pred size p on prey size
    % n
    
    cP = Meta.Params.Pred.cP; % Predator interference:
    
    Tp = Meta.Params.Timestep.*salP; % penalize total predation time
    
   % Crawley-Martin function:
   % Crowley-Martin functional response is best fit to mesocosm data.
   % Equation: dN/dt = a*N*P/(1 + b*N + c*P + b*c*N*P)
   AttackRate_t = (( Tp*aP'*Pt)) ./ (1 + hP.*Nt + cP.*repmat(Meta.IPM.Pred.Sy*Pt,[size(Nt,1),1]) + hP.*cP.*repmat(Meta.IPM.Pred.Sy*Pt,[size(Nt,1),1]).*Nt);
   DeadP_t = AttackRate_t.*Nt;
   Nt = Nt - DeadP_t;
    
if RunningMean > 30
  %  keyboard
end

    
    % Predators grow & die (depends on salinity)
    Pt = (kmatPm.*kmatPf)*Pt;
	Pt = kmatPg*Pt.*Meta.IPM.Pred.dy;

    AttackRate(tt,t) = sum(AttackRate_t);

    if saveopt
        Ntss(:,tt,t) = Nt;
    end
    
    end % end inner weekly loop over tt
        
    % store results for semiannual timestep:
    N(:,t) = Nt;
    P(:,t) = Pt;
    S(t) = St;

    DeadP(:,t) = DeadP_t;
    RR_s(:,t) = RR;
  
    

end % end loop over T


save(savename)
