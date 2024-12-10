%L Storch 2022 - a code to create a fake time series of salinity
%measurements, based on the methods in Denny et al 2009.  code that creates
%the climatology is TimeSeries_Resampling.  It takes a LONG TIME (~4 hrs?)
%to run, you should not need to run that code because all of the necessary
%information was saved in a .mat file which is loaded in this function

function salmat = Create_Mock_Salinity(stp)
%stp = 0; 

%input stp is how much we want to increase the standard deviation - should
%be a number between 0 and 1 where 0 means the standard dev stays the same
%and 1 is a doubling of the standard deviation (100% increase)

%output salmat is a mock salinity time series of 20 years of data in weekly
%measurements, where each column is 26 weeks of measurements (so 26 rows and
%40 cols) 

%this function should be in the same folder as the data set we are opening
%below:

%loading up the climatology year, the corresponding std dev, and the
%standardized residuals for the entire 20-year data set that was used to create
%MockSal and MockSalStd:
load("TimeSeries_Resampling_Output.mat","MockSal","MockSalStd","stdresid");

%MockSal is a climatology - a "standard" year of salinity measurements made
%by averaging over a 20 year data set (see "TimeSeries_Resampling.m" code
%on how it was done).  MockSalStd is the corresponding set of standard
%deviations for that climatology.  stdresid is the set of standardized
%residuals taken across the 20 year data set that was used to create
%MockSal

%we are using one week of data as our "segment" length, as the
%autocorrelation tends to rapidly decrease after about a week (but never
%drops to zero, hovers around 0.25).  measurements are taken every 30 min
%so a week of data is 336 data points

%n+48 is the same time the next day, so n+48*7 is the same day/time the
%following week (n + 336)

%To create a mock salinity time series: choose at random a segment of 336
%elements of stdresid.  multiply the first stdresid value by the first
%standard deviation value in MockSalStd, then add that to the first
%expected salinity value in MockSal.  Once this has been done for 336
%elements, choose another random 336 element chunk of stdresid and continue
%on through the time series:

%first, increase the standard deviation (or stay the same, if stp = 0)
NewStd = MockSalStd + MockSalStd*stp;

% Fix any stray nans
stdresid(isnan(stdresid))=nanmean(stdresid);

n = length(stdresid);
salcell = cell(1,20); %want 20 years of data

for j = 1:20 %want 20 years of data

    %need to create one year at a time of mock time series
    count = 1;
    for i = 1:52 

        if count > length(MockSal) %if we reach past the end of MockSal, stop loop
            break
        end

        %disp(strcat('year is ',num2str(j),' and week is ',num2str(i)));

        segstart = randi(n-335); %choose a random starting point for our segment of 336 elements
        stdresseg = stdresid(segstart:segstart+335); %grab a random 336 element chunk of standardized residuals

        for k = 1:length(stdresseg)
            saltimeseries(count) = MockSal(count)+NewStd(count)*stdresseg(k);

if any(isnan(saltimeseries(count)))
        keyboard
        end

            count = count +1;
        end

    end
    
    %if there are negative values, reset them to zero
    saltimeseries(saltimeseries<0) = 0;

    salcell{j} = saltimeseries;
end

%now we have to take the 40 years of 30 min increments and average them so
%it's a weekly data set with cols of 26 weeks each (so it jives with the
%format of data that the IPM model is expecting) 
count = 1; 
for ii = 1:length(salcell)
    for jj = 1:336:length(salcell{ii})-335
        salmat(count) = nanmean(salcell{ii}(jj:jj+335)); %average salinities for the week

        count = count + 1;

        
    end
end

% However we want to have half-years that correspond to Summer-Winter, so
% have to move the initial 'winter' weeks to the end of the timeseries

salmat = salmat(:);
salmat = [salmat(14:end); salmat(1:13)];

salmat = reshape(salmat,[26,40]); %reshape salmat into the correct format for the IPM




end

