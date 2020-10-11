%% Analyzing Neural Time Series Data
% Matlab code for Chapter 15
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 15.1

load sampleEEGdata

timewin = 500; % in ms

% convert ms to idx
timewinidx = round(timewin/(1000/EEG.srate));


% create hann taper function
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% detrend data (useful to attentuate super-low frequency artifacts in FFT
% from sampled data)
d = detrend(EEG.data(20,:,16));

figure
subplot(311)
plot(EEG.times,d)
title('One trial of data')

[junk,stime] = min(abs(EEG.times--50));

subplot(323)
plot(EEG.times(stime:stime+timewinidx-1),d(stime:stime+timewinidx-1))
hold on
plot(EEG.times(stime:stime+timewinidx-1),d(stime:stime+timewinidx-1).*hann_win,'r')
set(gca,'xlim',[-50 -50+timewin])
title('One short-time window of data, windowed')

dfft = fft(d(stime:stime+timewinidx-1).*hann_win);
f    = linspace(0,EEG.srate/2,floor(length(hann_win)/2)+1); % frequencies of FFT
subplot(313)
plot(f(2:end),abs(dfft(2:floor(length(hann_win)/2)+1)).^2,'.-');
title('power spectrum from that time window')
set(gca,'xlim',[1 128],'ylim',[-1000 25000],'xtick',0:10:EEG.srate/2)


% create TF matrix and input column of data at selected time point
tf = zeros(floor(length(hann_win)/2),EEG.pnts);
tf(:,stime+timewinidx/2-10:stime+timewinidx/2+10) = repmat(abs(dfft(2:floor(length(hann_win)/2)+1)')*2,1,21);

figure
imagesc(EEG.times,f,log10(tf+1))
set(gca,'clim',[-1 1]*4)

%% Figure 15.2

timewin        = 400; % in ms, for stFFT
times2save     = -300:50:1000; % in ms
channel2plot   = 'p7';
frequency2plot = 15;  % in Hz
timepoint2plot = 200; % ms

% convert from ms to index
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end
timewinidx = round(timewin/(1000/EEG.srate));
chan2useidx = strcmpi(channel2plot,{EEG.chanlocs.labels});


% create hann taper
hann_win = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% define frequencies
frex = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% initialize power output matrix
tf = zeros(length(frex),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    tempdat = squeeze(EEG.data(chan2useidx,times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2),:)); % note: the 'mod' function here corrects for even or odd number of points
    
    % taper data (using bsxfun instead of repmat... note sizes of tempdat
    % and hann_win)
    taperdat = bsxfun(@times,tempdat,hann_win');
    
    fdat = fft(taperdat,[],1)/timewinidx; % 3rd input is to make sure fft is over time
    tf(:,timepointi) = mean(abs(fdat(1:floor(timewinidx/2)+1,:)).^2,2); % average over trials
end

% plot
figure
subplot(121)
[junk,freq2plotidx]=min(abs(frex-frequency2plot));
plot(times2save,mean(log10(tf(freq2plotidx-2:freq2plotidx+2,:)),1))
title([ 'Sensor ' channel2plot ', ' num2str(frequency2plot) ' Hz' ])
axis square
set(gca,'xlim',[times2save(1) times2save(end)])

subplot(122)
[junk,time2plotidx]=min(abs(times2save-timepoint2plot));
plot(frex,log10(tf(:,time2plotidx)))
title([ 'Sensor ' channel2plot ', ' num2str(timepoint2plot) ' ms' ])
axis square
set(gca,'xlim',[frex(1) 40])

figure
contourf(times2save,frex,log10(tf),40,'linecolor','none')
set(gca,'clim',[-2 1])
title([ 'Sensor ' channel2plot ', power plot (no baseline correction)' ])

disp([ 'Overlap of ' num2str(100*(1-mean(diff(times2save))/timewin)) '%' ])

%% Figure 15.3

% create hamming
hamming_win = .54 - .46*cos(2*pi*(0:timewinidx-1)/(timewinidx-1));
hann_win    = .5*(1-cos(2*pi*(0:timewinidx-1)/(timewinidx-1)));

% create gaussian
gaus_win = exp(-.5*(2.5*(-timewinidx/2:timewinidx/2-1)/(timewinidx/2)).^2);

% plot together
figure
plot(hann_win)
hold on
plot(hamming_win,'r')
plot(gaus_win,'k')
legend({'Hann';'Hamming';'Gaussian'})

set(gca,'xlim',[-5 timewinidx+5],'ylim',[-.1 1.1],'ytick',0:.2:1)

%% figure 15.6

chan2use = 'o1';
frequency2plot = 10;  % in Hz


[junk,freq2plotidx]=min(abs(frex-frequency2plot));

% initialize ITPC output matrix
itpc = zeros(length(frex),length(times2save));

% loop over time points and perform FFT
for timepointi=1:length(times2save)
    
    % extract time series data for this center time point
    % (yes, yes, it's a long line of code. Perhaps you can understand it
    % better by breaking it up into several lines to separately identify
    % channel index and time window?)
    tempdat = squeeze(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),times2saveidx(timepointi)-floor(timewinidx/2):times2saveidx(timepointi)+floor(timewinidx/2)-mod(timewinidx+1,2),:));
    
    % taper data (Note the difference between using repmat here and bsxfun
    % earlier in the code.)
    taperdat = tempdat.*repmat(hann_win',1,EEG.trials);
    
    fdat = fft(taperdat,[],1)/timewinidx; % 3rd input is to make sure fft is over time
    itpc(:,timepointi) = abs(mean(exp(1i*angle(fdat(1:floor(timewinidx/2)+1,:))),2)); % average over trials
end

figure
contourf(times2save,frex,itpc,40,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-200 1000])
title([ 'ITPC at sensor ' chan2use ])

figure
plot(times2save,mean(itpc(freq2plotidx-2:freq2plotidx+2,:),1))
title([ 'ITPC at sensor ' chan2use ', ' num2str(frequency2plot) ' Hz' ])
set(gca,'xlim',[-200 1000])

%% end.
