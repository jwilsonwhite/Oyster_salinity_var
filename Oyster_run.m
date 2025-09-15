% Run the different scenarios for the oyster predator-prey model
function Oyster_run(F_pred,salopt)

if ~exist('salopt','var')
    salopt = 'none';
end

saveopt = 'none';

%tic;Salinity
Oyster_Params; % run, if needed, to get up-to-date structure of params & metadata
load('oyster_PP_params.mat') % load params & metadata ('Meta')
Meta.Rect = 'Closed'; % use 'Closed' recruitment for model runs
T = 40; % length of simulation. T = 40 --> 40 seasons/20 years (each column in TS_sal and TS_temp is 26wks/6mos)


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

Rn = ones(T,1)*1; % recruitment seed? For 'Open' and 'Mixed' recruitment only
F = ones(T,1)*(0/52); % level of fishing mortality (harvest) for oysters. 0 for Ch. 2 runs
Fp = ones(T,1)*(F_pred/52); % level of fishing mortality (culling/harvest) for drills. 0 for Ch. 2 runs


%runs = 2e3;
runs = 10;%
Y= 0; %[0 0.1 0.5]; 
name = 'Salinity_stdev_multiplied';

dosizeplot = true;
if dosizeplot
    figure (1)
    clf
    Cm = parula(runs);
end

%range of multipliers for standard deviation (0 means the standard deviation stays the same, 0.5 means 50% increase)
% Based on Pendergrass et al, showing a 0-5% increase in variability per
% degree K change for this region.

for i = 1:length(Y)
    for j = 1:runs
    
    load_name = strcat('Salinity_ts/', name, num2str(Y(i)), '_run_num_', num2str(j));
    load(load_name)

 %   TS_sal(TS_sal>36) = 36; % ceiling on salinity levels

    savename = strcat('Results_31July2025/Results_31July2025_F',num2str(F_pred),'_stdev_multiplied_',num2str(Y(i)),'_salopt_',salopt,'_run_num_',num2str(j), ...
        '.mat');
    N = Oyster_IPM_Disturbance(Meta,T,TS_temp,TS_sal,Rn,F,Fp,savename,saveopt,salopt);

if dosizeplot
    figure(1)
    hold on

    Ntmp = N(:,35:end);
    Ntmp = Ntmp./repmat(sum(Ntmp),[250,1]); % rescale to relative abundance

    plot(Meta.IPM.Prey.x,Ntmp*5,'color',Cm(j,:)) % rescale to match integration scale of the Apalachicola data
    xlim([0 150])
    ylim([0 0.3])
    %keyboard
end % end doplot

    end

end
 
