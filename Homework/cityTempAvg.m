function tempAvg = cityTempAvg(city, t)
%CITYTEMPAVG: HW 3 P03.004 - 2 - daily avg temp of city
% inputs:
%   city = from a given list(see ftn file, i.e. 'Miami')
%   t = days (vector). 1 = 1st day of yr... 365 = last
% output:
%   temp = daily average temp over time range t

addpath(['C:\Users\caleb\..B - School\_A - MATLAB Scripts\_' ...
    'ME 140 Matlab\Ch3 - Programming with MATLAB'])
validateInput(city, 'char');
validateInput(t(1), 'numeric', [1 365]);
validateInput(t(2), 'numeric', [1 365]);
validateInput(t(2)-t(1), 'numeric', [0 365]);

omega = 2*pi/365; % frequency of the annual variation
tpeak = 205;

switch city
    case 'Miami'
        Tmean = 22.1; % deg C
        Tpeak = 28.3;
    case 'Yuma'
        Tmean = 23.6;
        Tpeak = 33.6;
    case 'Bismarck'
        Tmean = 5.2;
        Tpeak = 22.1;
    case 'Seattle'
        Tmean = 10.1;
        Tpeak = 17.6;
    case 'Boston'
        Tmean = 10.7;
        Tpeak = 22.9;
    otherwise
        error('Please input a listed city.');
end

tempDay = Tmean + (Tpeak - Tmean)*cos(omega*(t - tpeak));
plot(t,tempDay);
xlabel('Day');
ylabel('Temperature');

tempAvg = mean(tempDay);
end