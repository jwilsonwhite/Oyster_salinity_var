function TimeSeries_Resampling
%L Storch 2022

%WARNING - THE MAIN PART OF THIS CODE TAKES 3 HOURS TO RUN!!!! 
%ALREADY COMPUTED OUTPUTS ARE IN TIMESERIES_RESAMPLING_OUTPUTS.MAT 

%edited July 2022 to add in temp information as well for Will's extreme
%disturbances project (Christian PhD project, REU student working on it)

%Trying to recreate the resampling method from Denny et al 2009 Eco Mono

%Have to process the files in CatPt_WaterQualityData, which are csv files
%of water quality from the Cat Point station in Apalachicola Bay from 2002
%through 2021.  Just want to grab salinity and date information for now

% Output is a salinity climatology + SD 

addpath("CatPt_WaterQualityData/") %tell matlab where to find the data files
%Data obtained from the NERR Centralized Data Management Office
%http://cdmo.baruch.sc.edu/aqs/index.cfm

for i = 2002:2021
    
    filename = strcat('apacpwq',num2str(i),'.csv');
    processtable = readtable(filename); %need a table to get the datetime stuff which is in a horrible format
    
    if i == 2002 %if it's the first run-through, have to set some things up
        TimeMat = table2array(processtable(:,3));
        TimeMat = datevec(TimeMat); %FORMAT: year month day hour minute second
        SalHours = table2array(processtable(:,11));
        TempHours = table2array(processtable(:,7));
    else %otherwise, vertically concatenate
        TimeMat = [TimeMat; datevec(table2array(processtable(:,3)))]; %FORMAT: year month day hour minute second
        SalHours = [SalHours; table2array(processtable(:,11))];
        TempHours = [TempHours; table2array(processtable(:,7))];
    end
    
    
end

%DateSal cols: (1)year (2)month (3)day (4)hour (5)minute (6)second (7)salinity
DateSal = [TimeMat, SalHours, TempHours]; %combine date matrix and sal matrix and temp matrix 
keyboard
%at some point in the data set they switched to sampling every 15 mins, get
%rid of the 15 and 45 min entries so everything is evenly spaced throughout
idx15 = find(DateSal(:,5) == 15);
DateSal(idx15,:) = [];
idx45 = find(DateSal(:,5) == 45);
DateSal(idx45,:) = [];

% find missing values and correct 
%matrix of actual dates, no missing values 
T1 = datetime(2002,1,1,0,0,0);
T2 = datetime(2021,12,31,23,30,0);
V = (T1:minutes(30):T2).';
DateMat = [V.Year,V.Month,V.Day,V.Hour,V.Minute]; %no missing values of year month day hour minute 

%according to DateMat, DateSal has one missing value 
DateSal(end+1,:) = NaN; %just making DateSal the same size as DateMat
idx = DateMat(:,4) == DateSal(:,4); %where are they the same? 
wheremissing = find(~idx,1,'first'); %where is the missing data point? 
DateMat(1:wheremissing-1,6) = DateSal(1:wheremissing-1,7); %fill in everything before missing 
DateMat(wheremissing+1:end,6) = DateSal(wheremissing+1:end,7); %fill in everything after missing 
DateMat(wheremissing,6) = NaN; %the missing one gets a NaN
DateMat(1:wheremissing-1,7) = DateSal(1:wheremissing-1,8); %fill in everything before missing 
DateMat(wheremissing+1:end,7) = DateSal(wheremissing+1:end,8); %fill in everything after missing 
DateMat(wheremissing,7) = NaN; %the missing one gets a NaN
%% creating the mock time series.  for each entry we want a mean and stdev
%we use a within-day window of 6hrs and a between-day window of 2 days

%data points are every 30 min so n-12 is -6 hours and n+12 is +6 hrs
%n+48 is the same time the next day, n-48 is the same time the previous day

n = length(DateSal);
MockSal = []; %mean values for a given day
MockSalStd = []; %std dev for a given day
daywindow = 2; %starting with +/- 2 days for between-day window
hourswindow = 12; %index if 12 is +/- 6 hrs
yearidx = 17520; % n+yearidx is the same time/day the next year (DOESN'T WORK FOR LEAP YEARS!) 

totalyears = 20; %time series goes from 2002 to 2021

count = 0;
tic
for i = (1+yearidx):(yearidx*2) %using 2003 as the year we run through so no issue with going out of bounds (and not a leap year so yearidx works) 
    
    count
    
    dataselect = [];
    dataselect = DateMat(i-hourswindow:i+hourswindow,6); %within-day of +/- 6 hrs from a given hour
    monthdayhourminute = DateMat(i,2:5); %grab month day hour minute of the given i
    match = strfind(num2cell(char(DateMat(:,2:5)),2),char(monthdayhourminute)); %look for a match to monthdayhourminute in each row of DateSal
    matchidx = find(~cellfun(@isempty,match)); %where is there a match?
    yearvec = matchidx(2:end-1); %grab all the indices for years of a given month day hour minute (chopping off the last index to avoid issues with going over the end of the data set)
    
    
    %     %both of these methods are SUPER SLOW (above and below) 
    %     idx = [];
    %     monthdayhourminute = DateSal(i,2:5); %grab month day hour minute of the given i
    %     for zz = 1+yearidx:length(DateSal)
    %         compare = isequal(DateSal(zz,2:5),monthdayhourminute); %this is dumb and clunky but i didn't know how to compare a vector against an entire matrix, so testing each row individually
    %         if compare
    %             idx = [idx; zz];
    %         end
    %     end
    %     yearvec = idx(1:end-1); %grab all the indices for years of a given month day hour minute (chopping off the last index to avoid issues with going over the end of the data set)
    %
    
    for k = 1:length(yearvec)
        
        yearidxx = yearvec(k);
        hourlyblock = DateMat(yearidxx - hourswindow: yearidxx + hourswindow,6); %grab a block of hours +/- 6 hours from the selected hour
        dataselect = [dataselect; hourlyblock]; %grabbing +/-6 hour window for same day in each year, aggregating
        
        for j = 1:daywindow %for a given year, stacking across +/- 2 day window 
            dayminus = yearidxx - 48*j;
            dayplus = yearidxx + 48*j;
            dataselect = [dataselect; DateMat(dayminus-hourswindow:dayminus+hourswindow,6); DateMat(dayplus-hourswindow:dayplus+hourswindow,6)]; %aggregating all of the data we will eventually average
        end
        
        
        
    end
    
    MockSal(1+count) = nanmean(dataselect); 
    MockSalStd(1+count) = nanstd(dataselect);
    count = count+1;
    
    
    
end
toc

%idx = find(MockSal >0);
%MockSal = MockSal(idx(1):idx(end)); %in case there are any legit zero values in the middle somewhere

%% ploting the mock time series over the original
subplot(1,2,1)
plot(MockSal)
hold on
plot(DateMat(1+yearidx:2*yearidx,6))
plot(DateMat(1:yearidx,6));
plot(DateMat(1+9*yearidx:10*yearidx,6))
legend('mocksal','2003','2002','2011')
subplot(1,2,2)
plot(MockSalStd)

%% calculate residuals and standardized residuals based on Denny procedure in main manuscript text (pg 399) - this one also takes 3 hours to run! 

MockSalTimes = DateMat(17521:35040,:);


%mocksal is for one year without a leap day so need to match the days/hours
%when calculating residuals 
resid = zeros(length(DateMat),1); %residual for a given date
stdresid = zeros(length(DateMat),1); %standardized residual for a given date 
tic
for i = 1:length(MockSal) %go through each day individually
    i
    monthdayhourminute = MockSalTimes(i,2:5); %grab month day hour minute of the given i
    match = strfind(num2cell(char(DateMat(:,2:5)),2),char(monthdayhourminute)); %look for a match to monthdayhourminute in each row of DateSal
    matchidx = find(~cellfun(@isempty,match)); %where is there a match?
    resid(matchidx) = (DateMat(matchidx,6) - MockSal(i));
    stdresid(matchidx) = resid(matchidx)/MockSalStd(i);
end
toc

%% autocorrelation - it takes about 20 days (1000 lags) to get down to ~0.25 and that's about as low as it goes  

[acf,lags] = autocorr(stdresid,NumLags=3000);

subplot(1,4,1)
plot(DateMat(yearidx:yearidx*2,6))
title('example year of data')
subplot(1,4,2)
plot(MockSal)
title('climatology')
subplot(1,4,3)
plot(MockSalStd)
title('standard deviation')
subplot(1,4,4)
plot(stdresid(yearidx:yearidx*2))
title('standard residuals of 1 year')

figure(2)
%plot(acf(1:1100))
plot(acf)
set(gca,'fontsize',20)
set(gcf,'color','white')
xlabel('time in days')
ylabel('autocorrection')
%xticks([0 100 200 300 400 500 600 700 800 900 1000 1100]);
%xtickdays = round([0 100 200 300 400 500 600 700 800 900 1000 1100]/48);
%xticklabels(xtickdays)

%% the autocorrelation doesn't go to zero but takes a dive after about a week.  we are just going to use a week as the "window" for a segment
%this is convenient because in the code the salinity measurements are
%weekly and so we will convert these to a weekly as well, and the length of
%a segment will just be one point in time

%n+48 is the same time the next day, so n+48*7 is the same day/time the
%following week (n + 336)

WeeklyMat = [];
WeeklyStd = [];

count = 1;
for i = 1:336:length(MockSal)-335
    WeeklyMat(count,1:4) = DateMat(i,1:4);
    WeeklyMat(count,5) = mean(MockSal(i:i+335)); %average salinities for the week 

    %usually we take the averages of variance, not stdev, so turn it into
    %variance, do the averaging, turn back into stdev 
    std = MockSalStd(i:i+335);
    std = std.^2;
    avgstd = mean(std);
    avgstd = sqrt(avgstd);
    WeeklyStd(count) = avgstd;
    count = count + 1;
end

%added a piece to get weekly temp readings out of DateMat
WeeklyTemp = [];
count = 1;
for i = 1:336:length(DateMat)-335
    WeeklyTemp(count) = nanmean(DateMat(i:i+335,7));
    count = count+1;
end

%there are 3 extra weeks so just chopping the last 3 off
WeeklyTemp(end-2:end) = [];
WeeklyTemp = reshape(WeeklyTemp,[26,40]);

save('CatPt_2002_to_2021_AvgWeeklyTemp',"WeeklyTemp")

subplot(1,2,1)
plot(WeeklyMat(:,5),'linewidth',2);
xlabel('Weeks')
ylabel('Salinity')
xlim([1 52])
set(gca,'fontsize',20)
set(gcf,'color','white')
subplot(1,2,2)
plot(WeeklyStd,'linewidth',2);
xlabel('Weeks')
xlim([1 52])
ylabel('Standard deviation')
set(gca,'fontsize',20)
set(gcf,'color','white')



