function Oyster_IPM_Disturbance_Analysis
% Oyster_IPM_Disturbance Results Analysis:



doFig1 = false; % salinity timeseries statistics
doFig2 = false;% example salinity effects on predator
doFig3 =  false; % example time series
doFig4 = true; % histograms of abundance, attack rate, etc. & Table of extinction probs


Months = {'J','F','M','A','M','J','J','A','S','O','N','D'};
%----------------------------
if doFig4
% Figure 1: Histograms of abundance, attack rate, etc. & Table of
% extinction probabilities
figure(4)
clf

Ys = [0 0.1 0.5]; % levels of increased variability in salinity (based on Pendergrass)
Fs = [0 0.13 0.26]; % levels of fishing on the predator
Nsims = 2e3;

for y = 1:length(Ys)
% 1) load data: e.g., load('y0_2000sims_FullRun1.mat')
% Output (columns in N_outer & P_outer cells) is density distribution of abundance,
% so to get total abundance, we have to integrate via Simpsons rule integration (multiply by Sy vector):
% Total abundance (N) = Meta.IPM.Prey.Sy * N; Total abundance (P) = Meta.IPM.Pred.Sy * P

% 2) loop through N_outer & P_outer to get total abundance for all sims

% pre-allocate:
N_TotalAbundance_mean_Last10T = zeros(1, Nsims);
N_TotalAbundance_stdev_Last10T = zeros(1, Nsims);

P_TotalAbundance_mean_Last10T = zeros(1, Nsims);
P_TotalAbundance_stdev_Last10T = zeros(1, Nsims);

for f = 1:length(Fs)
for i = 1:Nsims
    
% [change the file names to match what you are using]
fname = strcat('Results_31July2025/Results_31July2025_F',num2str(Fs(f)),'_stdev_multiplied_',num2str(Ys(y)),'_salopt_none_run_num_',num2str(i),'.mat');
load(fname,'Meta','N','P','AttackRate','RR_s')


N_TotalAbundance{i} = Meta.IPM.Prey.dy * N;

P_TotalAbundance{i} = Meta.IPM.Pred.dy * P;

N_AdultAbundance{i} = Meta.IPM.Prey.dy * (N .* (Meta.IPM.Prey.x(:) >= Meta.Params.Prey.Lf));

%P_AdultAbundance{i} = Meta.IPM.Pred.Sy * P_outer{1,i};

AR(y,f,i) = mean(mean(AttackRate(:,31:40))); % mean attack rate for final 10 years. 
% Divide by # size classes of N to get mean attack rate per oyster (instead of summed over all sizes as was calculated in the code)
% multiply by 52 to get annual rate
AR(y,f,i) = AR(y,f,i)/size(N,1)*52;

%R_total{i} = Meta.IPM.Prey.Sy * (N_outer{i}.*repmat(Meta.Params.Prey.isR(:),1,40));
R_total{i} = Meta.IPM.Prey.Sy * RR_s;


% 3) loop through total abundance cells to get mean and std. dev. across 
% last 10 semi-years (e.g., N_TotalAbundance{1,1}(:,31:40) ) across all cells/sims, 
% to determine # of sims needed: 


N_TotalAbundance_mean_Last10T(f,i) = mean(sum(N_TotalAbundance{i}(:,31:40)));
N_TotalAbundance_stdev_Last10T(f,i) = std(sum(N_TotalAbundance{i}(:,31:40)));

N_AdultAbundance_mean_Last10T(f,i) = mean(sum(N_AdultAbundance{i}(:,31:40)));
N_AdultAbundance_stdev_Last10T(f,i) = std(sum(N_AdultAbundance{i}(:,31:40)));

P_TotalAbundance_mean_Last10T(f,i) = mean(sum(P_TotalAbundance{i}(:,31:40)));
P_TotalAbundance_stdev_Last10T(f,i) = std(sum(P_TotalAbundance{i}(:,31:40)));

RR_total(y,f,i) = mean(R_total{i}(31:40)) * 2; % note x2 to get yearly rate


end


% get total mean and stdev across all sims:

N_TotalAbundance_TotalMean_Last10T(y,f) = mean(N_TotalAbundance_mean_Last10T(f,:));
N_TotalAbundance_TotalStDev_Last10T(y,f) = mean(N_TotalAbundance_stdev_Last10T(f,:)); % average St. Dev. ???

N_AdultAbundance_TotalMean_Last10T(y,f) = mean(N_AdultAbundance_mean_Last10T(f,:));
N_AdultAbundance_TotalStDev_Last10T(y,f) = mean(N_AdultAbundance_stdev_Last10T(f,:)); % average St. Dev. ???



P_TotalAbundance_TotalMean_Last10T(y,f) = mean(P_TotalAbundance_mean_Last10T(f,:));
P_TotalAbundance_TotalStDev_Last10T(y,f) = mean(P_TotalAbundance_stdev_Last10T(f,:)); % average St. Dev. ???



figure(4)
set(gcf,'units','centimeters','position',[10,10,18,20])

SPs = reshape(1:(4*length(Ys)),4,length(Ys))';
subplot(length(Ys),4,SPs(y,1))
hold on
histogram(N_AdultAbundance_mean_Last10T(f,:),0:10:400,'normalization','pdf','edgecolor','none')
%xlim([0 440])
xlabel('Adult abundance (oysters 0.25 m–2)')
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:100:500)

subplot(length(Ys),4,SPs(y,2))
hold on
histogram(N_AdultAbundance_stdev_Last10T(f,:)./N_AdultAbundance_mean_Last10T(f,:),0:0.01:0.5,'normalization','pdf','edgecolor','none')
%xlim([0 0.45])
xlabel('SD adult abundance (oysters 0.25 m–2)')
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:0.1:0.5)

subplot(length(Ys),4,SPs(y,3))
hold on
histogram(AR(y,f,:),0:0.005:0.2,'normalization','pdf','edgecolor','none')
%xlim([0 0.15])
xlabel('Predator attack rate (oyster-1 y-1)')
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])


subplot(length(Ys),4,SPs(y,4))
hold on
histogram(RR_total(y,f,:),0:10:550,'normalization','pdf','edgecolor','none') 
%xlim([100 550])
xlabel('Recruitment (spat 0.25 m-2 y-1')
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',100:100:500)

%----------------------------

% Probability of extinction (i.e., population collapse) calculated as the percentage 
% of simulations in which populations fall below 20% of the initial population abundance

% Prey:
    Extinctions = nan(length(N_TotalAbundance),1);
    
    if Ys(y) == 0 && Fs(f) == 0 % y0 case
        Equilib = N_AdultAbundance_TotalMean_Last10T;
    end

    Threshold = 0.5*Equilib;
    
for k = 1:length(N_AdultAbundance)
    X = sum(N_AdultAbundance{k});
 %   Equilib = 5.20; % equilibrium pop. size, i.e., avg. end pop size for y0 (it was initial oyster population size, 150)
    
    Extinctions(k,1) = any(X((end-10):end) <= Threshold);
end

    ExtinctProb(y,f) = mean(Extinctions);

    end % end loop over F
end % end loop over Y
keyboard
end % end if doFig2


%-----------------------------------------


%-----------------------------------------
% salinity regime characteristics (*change "y" as needed)
if doFig1

figure(1)
clf
set(gcf,'units','cent','Position',[10 10 8.5 12])

% Top panels: climatology and sd around climatology
load TimeSeries_Resampling_Output.mat MockSal MockSalStd

subplot(3,3,1:3)
hold on

Salsm = smooth(MockSal,336*2); % biweekly smoothing
Salsdsm = smooth(MockSalStd,336*2);% biweekly smoothing

plot(Salsm(:),'k-','LineWidth',1.5)
plot(Salsm(:)+Salsdsm(:),'k:','LineWidth',1.5)
plot(Salsm(:)-Salsdsm(:),'k:','LineWidth',1.5)
dm = length(MockSal(:))/12;
set(gca,'XTick',1:dm:length(MockSal(:)),'XTickLabel',Months)
set(gca,'YGrid','on','xlim',[0,length(Salsm)],'ylim',[0 39])
set(gca,'tickdir','out','xcolor','k','ycolor','k','fontsize',8)

% Add temp on the second axis
yyaxis right

load Temp_ts.mat
tempX = linspace(1,length(Salsm),52); % x-indices to plot weekly data on the same scale
Temp = [TS_temp(14:end,2); TS_temp(:,1); TS_temp(1:13,2)]; % one year
plot(tempX,Temp(:),'r-','LineWidth',1.5)
set(gca,'ycolor','r','ylim',[0 35])

subplot(3,3,4:6)
hold on

C = colororder; % default colororder matrix
Cols = C([1,3,4],:);
Cols(:,4) = 0.3; % hack to add transparency

% Middle panels: Example ts for a range of levels of SD
Num = 2;
Ys = [0 0.1 0.5];
for y = 1:length(Ys)
for j = 1:Num
Sal_filename0 = strcat('Salinity_ts/Salinity_stdev_multiplied',num2str(Ys(y)),'_run_num_',num2str(j),'.mat');
load(Sal_filename0,'TS_sal'); TS_sal = TS_sal(:);
TS_tmp = TS_sal(40:(40+104)); % 2 years, starting in January
TS_tmp(TS_tmp>=36)=36;
plot(TS_tmp,'k-','Color',Cols(y,:),'LineWidth',1.5);
end % end Num
end % end Ys

Ticks = linspace(1,104,25);
set(gca,'YGrid','on','xlim',[1,104],'xtick',Ticks(1:end-1),'XTickLabel',Months)
set(gca,'tickdir','out','ticklength',[0.015 0.015],'ylim',[0 45],'fontsize',4,'ytick',0:10:40)




% Bottom panels: Histograms of distributions
CalcSalStats = true; % only need to do this once
Ys = [0 0.05 0.1 0.5];
Nsims = 2e3;
if CalcSalStats
    for y = 1:length(Ys)
for i = 1:Nsims
    Sal_filename = strcat('Salinity_ts/Salinity_stdev_multiplied',num2str(Ys(y)),'_run_num_',num2str(i),'.mat');
    load(Sal_filename,'TS_sal');
    Summer = 1:2:40;
    Sals = TS_sal(:,Summer);
    Sals = Sals(:);
    Sals(Sals>=36)=36;

    sals_max(y,i) = quantile(Sals,0.95);
  %  sals_max(y,i) = max(Sals);
    sals_min(y,i) = quantile(Sals,0.05);
    sals_mean(y,i) = mean(Sals,'omitnan');
    sals_std(y,i) = std(Sals,'omitnan');

 %   sal_min_total = mean(sals_min{1,i});
 %   sal_mean_total = mean(sals_mean{1,i});
 %   sal_max_total = mean(sals_max{1,i});

 save CalcSalStats.mat sals_max sals_min sals_mean sals_std
end % end for Nsims
    end % end Ys
    else
load CalcSalStats.mat
    end % end if CalcSalStats

for y = 1:length(Ys)
subplot(3,3,8)
hold on
histogram(sals_mean(y,:),50,'Normalization','pdf','EdgeColor','none')
set(gca,'tickdir','out','ticklength',[0.015 0.015],'ytick',[],...
    'xcolor','k','fontsize',6)
set(gca,'xtick',20:1:24)

subplot(3,3,9)
hold on
histogram(sals_max(y,:),50,'Normalization','pdf','EdgeColor','none')
set(gca,'tickdir','out','ticklength',[0.015 0.015],'ytick',[],...
    'xcolor','k','fontsize',6)
set(gca,'xtick',28:2:36)

subplot(3,3,7)
hold on
histogram(sals_min(y,:),50,'Normalization','pdf','EdgeColor','none')
set(gca,'tickdir','out','ticklength',[0.015 0.015],'ytick',[],...
    'xcolor','k','fontsize',6)
set(gca,'xtick',3:2:12)

end

keyboard

end % end if doFig1



%%%%%%-----------------------------------------------------------
% Fig 2
if doFig2
% Representative plot of salinity vs. drill attack rate
% Scenarios to use: Ys = [0, 0.5]

Panels = [1,3,2,4];
Runs = [50, 50, 1, 1];
Summers = [33,35,39,37];
Ys = [0, 0, 0.5, 0.5];
% Panel 1: 'good' summer with fresh water
% Panel 3: 'bad' summer with salty water

figure(2)
clf
set(gcf,'units','centimeters','Position',[15,10, 8.5, 7 ])

for p = 1:length(Panels)
fname = strcat('Results_31July2025/Results_31July2025_F0_stdev_multiplied_',num2str(Ys(p)),'_salopt_none_run_num_',num2str(Runs(p)),'.mat');
load(fname,'Meta','N','P','AttackRate','RR_s')
sname = strcat('Salinity_ts/Salinity_stdev_multiplied',num2str(Ys(p)),'_run_num_',num2str(Runs(p)),'.mat');
load(sname,'TS_sal')
load('Temp_ts.mat','TS_temp');

TS_sal(TS_sal>36)=36;

% emulate the runningmean variable
RunningMeanVec1 = TS_sal(:);
RunningMeanVec1 = [RunningMeanVec1(1); RunningMeanVec1(1:end-1)];
RunningMeanVec2 = TS_sal(:);
RunningMean = mean([RunningMeanVec1, RunningMeanVec2],2);

for s = 1:length(RunningMean)
Penalty(s) = predator_salinity_penalty(TS_sal(s),TS_temp(s),RunningMean(s));
end
Penalty = reshape(Penalty,[26,40]);

subplot(2,2,Panels(p));
hold on

TS_plot = TS_sal(:,Summers(p):(Summers(p)+1));
Weekvals = 14:(14+51);
Weekticks = linspace(14,14+51,13);
Weeklab = [Months(4:12) Months(1:3)];
plot(Weekvals,TS_plot(:),'k-','LineWidth',1.5)

set(gca,'tickdir','out','ticklength',[0.015,0.015],'xcolor','k','ycolor','k','color','none')
set(gca,'ylim',[0 39],'xlim',[14,65],'xtick',Weekticks(1:end-1),'xticklabels',Weeklab,'fontsize',6)

Pen_tmp = Penalty(:,Summers(p):(Summers(p)+1));

Pen_store(:,p) = Pen_tmp(:); % save for later use

yyaxis right
plot(Weekvals,Pen_tmp(:),'r-','LineWidth',1.5,'Color',[0.8 0.1 0.1 0.5])
plot([14,39],[quantile(Pen_tmp(1:26),0.25),quantile(Pen_tmp(1:26),0.25)],'r--','LineWidth',0.5)
plot([14,39],[quantile(Pen_tmp(1:26),0.75),quantile(Pen_tmp(1:26),0.75)],'r--','LineWidth',0.5)
plot([40,66],[quantile(Pen_tmp(27:52),0.25),quantile(Pen_tmp(27:52),0.25)],'r--','LineWidth',0.5)
plot([40,66],[quantile(Pen_tmp(27:52),0.75),quantile(Pen_tmp(27:52),0.75)],'r--','LineWidth',0.5)

set(gca,'ycolor','k','YTick',0:0.2:1,'ylim',[-0.03 1.1])

end

end


%%%%%%%%%%
% Figure 3 
% Timeseries of oyster population under particular scenarios
if doFig3   
doRuns = true;

figure(3)
clf
set(gcf,'units','centimeters','Position',[20, 20, 16, 12])

figure(4)
clf
set(gcf,'units','centimeters','Position',[10, 15, 12, 10])


figure(5)
clf

        % Use the same runs as in Fig 2:
        Runs = [1, 2];
        Summers = [31, 37; 31,37];%[33,35; 39,37];
        Ys = [0, 0.5];
        Fs = [0,0.26];
        Subplots=[1,3;2,4];
    if doRuns % if we need to rerun scenarios to get the data

        load('oyster_PP_params.mat') % load params & metadata ('Meta')
        Meta.Rect = 'Closed';
        T = 40;
        Rn = ones(T,1)*1; % recruitment seed? For 'Open' and 'Mixed' recruitment only

        load('Temp_ts.mat','TS_temp');

        for j = 1:length(Runs)
            for f = 1:length(Fs)

            F = ones(T,1)*(0/52); % level of fishing mortality (harvest) for oysters. 0 for Ch. 2 runs
            Fp = ones(T,1)*(Fs(f)/52); % level of fishing mortality (culling/harvest) for drills. 0 for Ch. 2 runs

        
        sname = strcat('Salinity_ts/Salinity_stdev_multiplied',num2str(Ys(j)),'_run_num_',num2str(Runs(j)),'.mat');
        load(sname,'TS_sal')
        TmpSaveName = strcat('Exrun_NoLarval_sd',num2str(Ys(j)),'_Fp',num2str(Fs(f)),'_run',num2str(Runs(j)),'.mat');

        Oyster_IPM_Disturbance(Meta,T,TS_temp,TS_sal,Rn,F,Fp,TmpSaveName,true,'NoLarval')
            end % end Fs
        end % end Runs
    end % end if doRuns

   % Load save file and make plots
        for j = 1:length(Runs)

            sname = strcat('Salinity_ts/Salinity_stdev_multiplied',num2str(Ys(j)),'_run_num_',num2str(Runs(j)),'.mat');
        load(sname,'TS_sal')


            for f = 1 %:length(Fs) % just plot with no predation, the difference in scale is confusing

                        TmpSaveName = strcat('Exrun_Fig3_sd',num2str(Ys(j)),'_Fp',num2str(Fs(f)),'_run',num2str(Runs(j)),'.mat');
                        load(TmpSaveName)


                % Loop over summers in a given run to make plots
                for s = 1:2

                    Nplot = Ntss(:,:,Summers(j,s));
                    sum(RR_s(:,Summers(j,s)))

                    % Filter to include only 'large adult' oysters
                    Nplot = Meta.IPM.Prey.dy * (Nplot(Meta.IPM.Prey.x(:) >= Meta.Params.Prey.Lf,:));

                    figure(3)
                    subplot(2,2,Subplots(j,s))
                    hold on

                    Np = sum(Nplot);
                    Np = Np/Np(1);
                    plot(1:26,Np)
                    set(gca,'xtick',1:4:26,'xticklabels',Months(4:10),'xcolor','k','ycolor','k',...
                    'tickdir','out','ticklength',[0.015 0.015],...
                    'xlim',[1 26],'ylim',[0.95,1.06],'ytick',0.95:0.025:1.05);
                   

                    if j == 1
                %        ylim([0.93 1])
                    end
                   

                  if f == 1
                  yyaxis right

                  plot(TS_sal(:,Summers(j,s)),'k-')
                  set(gca,'ycolor','k','ylim',[0 39])
                  end

                    %set(gca,'ygrid','on','xgrid','on')
                    %ylim([0 80])

                end % end loop over summers

                figure(4)

                subplot(2,1,j)
                hold on
                yyaxis left
                R_index = 1:26:(26*40);
                Xtick = -13:52:(52*19);
              %  Xticklabs = 1:10;
              Xticklabs = 1:2:40;
           %     keyboard
                plot(R_index,sum(RR_s))
                set(gca,'xtick',Xtick,'xticklabels',Xticklabs,'xcolor','k',...
                    'tickdir','out','ticklength',[0.015 0.015]);
               set(gca,'xlim',[503,1027])

                yyaxis right
                plot(TS_sal(:))

                figure(5)
                subplot(2,1,j)
                hold on
                plot(sum(N(Meta.IPM.Prey.x>76,:)))

                

               % keyboard

                end % end Fs
        end % end Runs

end % end if doFig3
%%%%%%%%%%

dolast = false;
if dolast
OysterBiomass = Meta.Params.Prey.LW * N_outer{1,i};
OysterBiomass = log(OysterBiomass);

DrillBiomass = Meta.Params.Pred.LW * P_outer{1,i};
DrillBiomass = log(DrillBiomass);

% change "y" as needed
sal = csvread('Sal_ExtremeSummer_HistoricalWinter_WeeklySalinity_timeseries_y0_i.csv');
sal_max = max(sal);
sal_min = min(sal);
sal_mean = mean(sal);

figure(2)
plot(OysterBiomass)
hold on
plot(DrillBiomass)
plot(sal_min)
legend('Oyster Biomass (log scale)','Drill Biomass (log scale)','Min Salinity (ppt)')
hold off

% baseline/historical salinity & temperature plots

% pop. time series
figure(3)
plot(log(N_TotalAbundance{1,i}));
hold on
plot(log(P_TotalAbundance{1,i}));
%ylim([0 10])
hline = refline(0,log(0.26)); % extinction threshold
hline.Color = 'k';
legend('Oyster Abundance (log scale)','Drill Abundance (log scale)','Extinction Threshold')
%plot(sal_min)
%plot(sal_mean)
%plot(sal_max)
hold off
end