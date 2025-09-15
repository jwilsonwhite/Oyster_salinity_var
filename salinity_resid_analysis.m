function salinity_resid_analysis

% Analysis of weekly vs higher-resolution averaging of salinity timeseries
% in Apalachicola Bay

% Load the climatological timeseries
load("TimeSeries_Resampling_Output.mat");


% Each week is 336 datapoints
X = MockSal(1:(336*52)); % trim out the residual week
X2 = reshape(X,[336,52]);
X3 = mean(X2);
X4 = repmat(X3,[336,1]);
X5 = X4(:);

Resid = X(:) - X5(:);

max(X)
quantile(X,0.95)

figure(1)
clf

subplot(3,2,1)
hold on
plot(X)
plot(X5,'linewidth',1.5)
set(gca,'xtick',1:(336*4):(336*52),'xticklabel',1:4:52,'xlim',[0 336*52])
xlabel('Weeks')
ylabel('Salinity')

subplot(3,2,3)
histogram(Resid,'edgecolor','none','normalization','probability');
set(gca,'ytick',[])
xlabel('Salinity deviation')
ylabel('Probabililty density')
xlim([-8 8])

subplot(3,2,5)
histogram(X5(:)-mean(X5),'edgecolor','none','normalization','probability');
set(gca,'ytick',[])
xlabel('Salinity deviation')
ylabel('Probabililty density')
xlim([-8 8])


% Full dataset
Y = DateMat(1:(336*52*20),6);  % trim out the residual week
Y2 = reshape(Y,[336,52*20]);
Y3 = mean(Y2);
Y4 = repmat(Y3,[336,1]);
Y5 = Y4(:);

ResidY = Y(:) - Y5(:);



subplot(3,2,2)
hold on
plot(Y)
plot(Y5,'linewidth',1.5)
set(gca,'xtick',1:(336*52*2):(336*52*20),'xticklabel',0:2:20)
xlabel('Years')
ylabel('Salinity')

subplot(3,2,4)
histogram(ResidY,'edgecolor','none','normalization','probability');
set(gca,'ytick',[])
xlabel('Salinity deviation')
ylabel('Probabililty density')
xlim([-20 20])

subplot(3,2,6)
histogram(Y5(:)-nanmean(Y5),'edgecolor','none','normalization','probability')
set(gca,'ytick',[])
xlabel('Salinity deviation')
ylabel('Probabililty density')
xlim([-20 20])

% FFT on the data to look at the time scale for highest variance


% Interpolate to remove NaNs:

X = 1:length(Y);
    warning('off','MATLAB:interp1:NaNstrip') % turn off inevitable warning 
    It = interp1(X,Y,X(isnan(Y)),'pchip'); % shape-preserving piecewise cubic interpolation
    % Insert interpolated values:
    Y(isnan(Y)) = It;


    FFT_S = detrend(Y);
    
    Fs = 48; % sampling frequency (48x per day)
    T = 1/Fs; % sampling interval
    L = length(FFT_S); % length of timeseries
    
    % find number of points for FFT
    % (FFT is most efficient for a number of points that is an integer power of 2)
    NFFT = 2^nextpow2(L);
    
    % run the FFT
    Y = fft(FFT_S,NFFT)/L;
    
    % vector for plotting
    % cannot detect frequencies greater than Fs/2
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    % plot one-sided FFT
    figure(2)
    clf

    hold on
    plot(log10(1./f),2*abs(Y(1:NFFT/2+1)));
    xlabel('Period (days)')
    ylabel('Power (variance)')
    set(gca,'tickdir','out','ticklength',[0.01 0.01])
    % x labels:
    Xl = -1:1:4;
    Xlab = 10.^Xl;
    set(gca,'xtick',Xl,'xticklabel',Xlab)
  %  end % end if Plot2
    

  [acf,lags] = autocorr(FFT_S,NumLags=1e4);
  % Calculate and plot the autocorrelation function
  laglab = lags/48;

figure(3)
clf
plot(laglab, acf, 'k-');
xlabel('Lag (days)');
ylabel('Autocorrelation');
title('Autocorrelation Function of Salinity Time Series');
xlim([0 100]);
set(gca,'xtick',0:7:140)

keyboard
  [pacf] = parcorr(FFT_S,NumLags=1e2);
  % Calculate and plot the autocorrelation function


figure(4)
clf
plot(laglab, pacf, 'k-');
xlabel('Lag (days)');
ylabel('Partial autocorrelation');
title('Autocorrelation Function of Salinity Time Series');
xlim([0 100]);
set(gca,'xtick',0:7:140)