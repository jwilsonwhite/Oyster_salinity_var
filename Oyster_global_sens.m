% Run a global sensitivity analysis different scenarios for the oyster predator-prey model
% 'CV' specifies the level of variability to use
function Oyster_global_sens(CV)

if ~exist('salopt','var')
    salopt = 'none';
end

F_pred = 0; % no fishing on drills for this analysis

% 20yr timeseries of historical 
% weekly temps; reflects historical temperature regime of estuary
load("CatPt_2002_to_2021_AvgWeeklyTemp.mat")
% Convert to climatology, and then reset to align with April-Oct timeframe
TS_temp = reshape(WeeklyTemp,[52,20]);
TS_temp_mean = nanmean(TS_temp,2);
TS_temp_mean = [TS_temp_mean(14:end); TS_temp_mean(1:13)];
TS_temp = repmat(TS_temp_mean,[1,20]);
TS_temp = reshape(TS_temp,[26,40]);
save('Temp_ts.mat','TS_temp')

T = 40; % length of simulation. T = 40 --> 40 seasons/20 years (each column in TS_sal and TS_temp is 26wks/6mos)


Rn = ones(T,1)*1; % recruitment seed? For 'Open' and 'Mixed' recruitment only
F = ones(T,1)*(0/52); % level of fishing mortality (harvest) for oysters. 0 for Ch. 2 runs
Fp = ones(T,1)*(F_pred/52); % level of fishing mortality (culling/harvest) for drills. 0 for Ch. 2 runs


N_sens_runs = 1e4; % how many runs to do for the global sensitivity
Results_mat = nan(N_sens_runs,17); % one row per set of simulations, one column for parameters + result

for n = 1:N_sens_runs

if mod(n,100)==0
    disp(n)
end


%tic;Salinity
Xvar = Oyster_Params_Sens(CV); % Get a randomly varied structure of params & metadata. The varying values are stored in Xvar
load('oyster_PP_params_sens.mat') % load params & metadata ('Meta')
Meta.Rect = 'Closed'; % use 'Closed' recruitment for model runs



runs = 1e1; % just take the average of 10 runs
Y= 0; %[0 0.1 0.5]; 
name = 'Salinity_stdev_multiplied';

%Xvar

    parfor j = 1:runs
    
    load_name = strcat('Salinity_ts/', name, num2str(Y), '_run_num_', num2str(j));
    Tmp = load(load_name,'TS_sal');
    TS_sal = Tmp.TS_sal;

    TS_sal(TS_sal>=36) = 36; % ceiling on salinity levels

  %  savename = strcat('Global_Sensitivity_31July2025/Results_16Oct2024_F',num2str(F_pred),'_stdev_multiplied_',num2str(Y(i)),'_salopt_',salopt,'_run_num_',num2str(j), ...
  %      '.mat');
  savename = strcat('Global_Sensitivity_31July2025/dummy.mat');
  Ntmp = Oyster_IPM_Disturbance(Meta,T,TS_temp,TS_sal,Rn,F,Fp,savename,salopt);


  Ntmp2(j) = mean(sum(Ntmp(:,20:T)));



    end

    % Summarize the results in terms of N, and store in matrix for later
% analysis

Results_mat(n,1) = mean(Ntmp2);
Results_mat(n,2:end) = Xvar;

end % %end loop over global sensitivity runs

%keyboard

% Variable order in Xvar
% Prey.k
% Prey.Linf
% Prey.Mj
% Prey.M
% Prey.Mat
% Prey.Fec.coeff
% Prey.DD
% prey.b2 (larval salinity dependence)
% prey.b4
% prey.b5a
% prey.b4a
% pred.M
% pred.fec.coeff
% pred.LEP
% pred.aP

savename = strcat('Global_Sensitivity_31July2025/Global_sens_results.mat');
save(savename,'Results_mat');
