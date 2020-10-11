%% Analyzing Neural Time Series Data
% Matlab code for Chapter 6
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 6.2

% create sine wave
srate     = 1000;
time      = 0:1/srate:1;
frequency = 3;

sinewave = sin(2*pi*frequency.*time);

figure
subplot(311)
plot(time,sinewave,'r')
set(gca,'xlim',[-.05 time(end)*1.05],'ylim',[-1.1 1.1])
hold on
sampling1 = round(linspace(1,length(time),frequency*2));
plot(time(sampling1),sinewave(sampling1),'o')
title('continuous sine wave')

sampling2 = round(linspace(1,length(time),frequency*20));
plot(time(sampling2),sinewave(sampling2),'+')


subplot(312)
plot(time(sampling1),sinewave(sampling1),'-o')
set(gca,'xlim',[-.05 time(end)*1.05],'ylim',[-1.1 1.1])
title('sampled at 2*frequency')

subplot(313)
plot(time(sampling2),sinewave(sampling2),'-+')
title('sampled at 20*frequency')
set(gca,'xlim',[-.05 time(end)*1.05],'ylim',[-1.1 1.1])

%% end.
