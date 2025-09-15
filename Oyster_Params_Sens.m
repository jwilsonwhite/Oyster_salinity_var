%Parameters for the oyster predator-prey model
function Xvar = Oyster_Params_Sens(CV)

% This version varies parameters randomly in order to conduct global
% sensitivity analysis. The randomly varied variables are output in Xvar

% 'CV' specifies how much variability to add

Meta.Params.Timestep = 1;

%%%%%%%%%%%%%%%%%%%%%%%
% IPM parameters - do not alter
Meta.IPM.Prey.meshmin = 0;
Meta.IPM.Prey.meshmax = 240; % approx 2x Linf
Meta.IPM.Prey.meshsize = 250;
Meta.IPM.Prey.x = linspace(Meta.IPM.Prey.meshmin,Meta.IPM.Prey.meshmax,Meta.IPM.Prey.meshsize);
Meta.IPM.Prey.meshdiff = diff(Meta.IPM.Prey.x(1:2));
Meta.IPM.Prey.dy = Meta.IPM.Prey.meshdiff;
Meta.IPM.Prey.Sy = makeSimpVec(Meta.IPM.Prey.dy,Meta.IPM.Prey.meshsize);
Meta.IPM.Prey.Symat = repmat(Meta.IPM.Prey.Sy(:)',[Meta.IPM.Prey.meshsize,1]);

% Predator:
Meta.IPM.Pred.meshmin = 0;
Meta.IPM.Pred.meshmax = 120;
Meta.IPM.Pred.meshsize = 100;
Meta.IPM.Pred.meshmin = 0;
Meta.IPM.Pred.x = linspace(Meta.IPM.Pred.meshmin,Meta.IPM.Pred.meshmax,Meta.IPM.Pred.meshsize);
Meta.IPM.Pred.meshdiff = diff(Meta.IPM.Pred.x(1:2));
Meta.IPM.Pred.dy = Meta.IPM.Pred.meshdiff;
Meta.IPM.Pred.Sy = makeSimpVec(Meta.IPM.Pred.dy,Meta.IPM.Pred.meshsize);
Meta.IPM.Pred.Symat = repmat(Meta.IPM.Pred.Sy(:)',[Meta.IPM.Pred.meshsize,1]);




%%%%%%%%%%%%%%%%%%%%%%%
%Oyster PARAMETERS - a lot of these were taken from oyster_PP_params or EPR_NERR or Christian's code
Xvar(1) = max(realmin,normrnd(0.0122,0.0122*CV));
Meta.Params.Prey.k = Xvar(1); %taken from Christian's code - weekly 
Xvar(2) = max(realmin,normrnd(120.3500,120.3500*CV));
Meta.Params.Prey.Linf = Xvar(2); %taken from Christian's code - mms
Meta.Params.Prey.vBs = 0.1; % cv of length-at-age
Xvar(3) = max(realmin,normrnd(0.0186,0.0186*CV));
Meta.Params.Prey.Mj = Xvar(3); %taken from Christian's code - Based on 80% survival over 12 weeks in cages from Kimbro experiment A2ii
Xvar(4) = max(realmin,normrnd( 0.0052, 0.0052*CV));
Meta.Params.Prey.M = Xvar(4); %taken from Christian's code - Based on 94% survival over 12 weeks in cages from Kimbro experiment A1v
Meta.Params.Prey.j_length = 15; %juvenile mortality size limit - taken from Christian's 
Xvar(5) = max(realmin,normrnd( 35, 35*CV));
Meta.Params.Prey.Mat = Xvar(5); %maturity parameter, taken from oyster_PP_params

Meta.Params.Prey.LW = 5.77e-4*Meta.IPM.Prey.x.^2.234; % this is WET SHELL WEIGHT; units = g ; derived from Kimbro field samples in salinity zone 2
Meta.Params.Prey.LW_afdw = 5.09e-5*Meta.IPM.Prey.x.^2.365; % Ash Free Dry Weight; units = g ; derived from Kimbro field samples in salinity zone 2 (taken from oyster_PP_params)
Meta.Params.Prey.Lf = 76; % legal fishing size limit (mm)
Meta.Params.Prey.Lfs = 1; % SD of selectivity ogive

Meta.Params.Prey.Rmean = 2.5; % mean spat size in week 1 - taken from Christian's code 
Meta.Params.Prey.Rstd = 0.5; % std in spat size - taken from Christian's code 
Meta.Params.Prey.Rvec = normpdf(Meta.IPM.Prey.x,Meta.Params.Prey.Rmean,Meta.Params.Prey.Rstd); %taken from EPR_NERR, initializing recruitment vector
Meta.Params.Prey.isR = normcdf(Meta.IPM.Prey.x,Meta.Params.Prey.Rmean,Meta.Params.Prey.Rstd)<0.99;

Meta.Params.Prey.Rnum = [4.75 11.2]; % % mean # recruits (oysters < 15) per quad. For use in model with external recruitment

Xvar(6) = max(realmin,normrnd( 19.86e4, 19.86e4*CV));
Meta.Params.Prey.Fec = Xvar(6) .* Meta.Params.Prey.LW_afdw.^(1.17); %fecundity parameter, taken from oyster_PP_params  ##### Adjusted downward two orders of magnitude (in the proportionality coefficient) to match Cox & Mann 1992, then further reduced to account for egg to veliger transition
Meta.Params.Prey.Fec = Meta.Params.Prey.Fec .* (Meta.IPM.Prey.x>=35); %Impose size at maturity of 35 mm
Meta.Params.Prey.SelfR = 1; % probability of self recruitment

Xvar(7) = max(realmin,normrnd( 0.75, 0.75*CV));
Xvar(8) = max(realmin,normrnd( 1.5e-3, 1.5e-3*CV));
Meta.Params.Prey.DDa = Xvar(7);
Meta.Params.Prey.DD = Xvar(8)*4; % based on Puckett & Eggleston (2012), Fig 7. (taken from oyster_PP_params), rescaled to 0.25 m2 scale by factor of 4

%salinity based larval mortality parameters - (Deprecated)
Meta.Params.Prey.M_r = 7.8; % 0.26 per day * 30 d (Rumrill 1990)
Meta.Params.Prey.M_rs = 2.934; % fitted param from SSIPM
Meta.Params.Prey.M_rs2 = 15; % optimum salinity (Davis 1958 Biol Bull)

%new salinity based mortality coefficients from Lough 1975, 2 day survival
Xvar(9) = max(realmin,normrnd( 32.8617, 32.8617*CV));
Xvar(10) = min(-realmin,normrnd( -0.6234, 0.6234*CV));
%(part 1 of the salinity penalty)
Meta.Params.Prey.b0 = -643.9149; %constant
Meta.Params.Prey.b1 = 27.7755; %T
Meta.Params.Prey.b2 = Xvar(9); %S
Meta.Params.Prey.b3 = -0.5195; %T^2
Meta.Params.Prey.b4 = Xvar(10); %S^2
Meta.Params.Prey.b5 = -0.0971; %T*S

Xvar(11) = max(realmin,normrnd( 2.6621, 2.6621*CV));
Xvar(12) = min(-realmin,normrnd(  -0.1455,  0.1455*CV));
%(part 2 of the salinity penalty)
Meta.Params.Prey.b0a = -104.389; %constant
Meta.Params.Prey.b1a = 10.221; %T
Meta.Params.Prey.b2a = -0.2006; %T^2
Meta.Params.Prey.b3a = 0.0996; %T*S
Meta.Params.Prey.b4a = Xvar(12); %S^2
Meta.Params.Prey.b5a = Xvar(11); %S




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%predator parameters
Meta.Params.Pred.k = 0.1; % units = 1 week; based on info in Butler (1985), see notebook p. 10
Meta.Params.Pred.Linf = 60;  % units = mm ****
Meta.Params.Pred.vBs = 0.1; % CV spread for growth
Xvar(13) = normrnd(  0.0025,  0.0025*CV);
Meta.Params.Pred.M = Xvar(13); % units = 1 week. Based on 12% annual loss rate (Butler 1985)
Meta.Params.Pred.Mj = Xvar(13);
Meta.Params.Pred.j_length = 21; %this doesn't matter because Mj and M are the same - source: Butler, P. A. (1985). Synoptic review...southern oyster drill
Meta.Params.Pred.Lf = 0; % size at first harvest
Meta.Params.Pred.Lfs = 1; % sd of selectivity ogive

Meta.Params.Pred.LW = NaN; % deprecated

Meta.Params.Pred.Mat = 25; % 2.5 cm, smallest size observed producing viable embryos, (Butler 1985)

Meta.Params.Pred.Rmean = 2; % mean drill size in week 1
Meta.Params.Pred.Rstd = 0.5; % std in spat size
Meta.Params.Pred.Rvec = normpdf(Meta.IPM.Pred.x,Meta.Params.Pred.Rmean,Meta.Params.Pred.Rstd); %changed March 2021, testing 

%Pred.Fec = 43350 * Params.x.^(0.782); % maybe 3e5 .* Meta.IPM.Pred.x.^(0.782); ???
Xvar(14) = max(realmin,normrnd(  43350, 43350*CV));
Meta.Params.Pred.Fec = Xvar(14) * Meta.IPM.Pred.x.^(0.782); %changed March 2021, testing 

Xvar(15) = max(realmin,normrnd( 2332400, 2332400*CV));
Meta.Params.Pred.LEP = Xvar(15); % lifetime egg production

% Values from Pusack et al. 2018 - mesocosm experiment with Crowley-Martin
% functional response
Xvar(16) = max(realmin,normrnd( 0.2068, 0.2068*CV));
aP = Xvar(16) * 7; % units in paper are per day, convert to weekly
Meta.Params.Pred.cP = 0.2162;
Meta.Params.Pred.hP = 0.5142;

% Adjust attack rate based on size preference of drills vs oysters.
% (Mesocosm, Report Fig. 23A)
PrefMat = normpdf(repmat(Meta.IPM.Prey.x(:)',[length(Meta.IPM.Pred.x),1]),...
    (39.33 + 0.367.*repmat(Meta.IPM.Pred.x(:),[1,length(Meta.IPM.Prey.x)])),235^0.5);
% Should sum to one for each predator size class:
PrefMat = PrefMat./repmat(sum(PrefMat,2),[1,size(PrefMat,2)]);
%Meta.Params.Pred.aP = 3.8041.*PrefMat; % attack rate (size dependent)
Meta.Params.Pred.aP = aP.*PrefMat;

Meta.Params.lambdaTAF = 0.1/52;  % Based on annual rate of 0.1 from Powell et al. (2012)


save('oyster_PP_params_sens.mat')

