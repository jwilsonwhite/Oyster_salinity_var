function Apalach_data_comparison

% Create histogram of Apalachicola data for comparison to model runs

% Load data
% Dataset is the FDACS Apalachicola oyster size survey
load meta_IPM_SS_CP_Historical_Intera_17May2016.mat



Size = Meta.Data.Edges; % size bins
N = Meta.Data.Nsamp;
W = Meta.Data.ModelWeek;
Minweek= 13;
D = Meta.Data.Data(Size<76,W>=Minweek); % sublegal size oysters in the latter half of the model season

D2 = Meta.Data.Data;
Dsc = D2./repmat(N(:)',[28,1]); % density
Dsc2 = Dsc./repmat(sum(Dsc),[28,1]);

histogram(sum(D./N(W>=Minweek)),15,'normalization','pdf','facecolor','none')
%set(gca,'xlim',[0,110])

figure;
plot(Size,Dsc2(:,W>=Minweek),'color',[0.5 0.5 0.5])
hold on
plot(Size,mean(Dsc2,2),'k-','linewidth',2)

keyboard
