%% Analyzing Neural Time Series Data
% Matlab code for Chapter 16
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 16.1

load sampleEEGdata

channel2plot = 'o1';
timewin      = 400; % in ms


timewinidx = round(timewin/(1000/EEG.srate));
tapers     = dpss(timewinidx,5); % this line will crash without matlab signal processing toolbox

% extract a bit of EEG data
d = detrend(squeeze(EEG.data(strcmpi(channel2plot,{EEG.chanlocs.labels}),200:200+timewinidx-1,10)));

% plot EEG data snippet
figure
subplot(5,2,1)
plot(d)
axis tight,axis off

% plot tapers
for i=1:5
    subplot(5,2,(2*(i-1))+2)
    plot(tapers(:,i))
    axis tight,axis off
end

% plot taper.*data
figure
for i=1:5
    subplot(5,2,(2*(i-1))+1)
    plot(tapers(:,i).*d')
    axis tight,axis off
end

% plot fft of taper.*data
f=zeros(5,timewinidx);
for i=1:5
    subplot(5,2,(2*(i-1))+2)
    f(i,:)=fft(tapers(:,i).*d');
    plot(abs(f(i,1:timewinidx/2)).^2)
    axis tight,axis off
end

figure
subplot(5,2,2)
plot(mean(abs(f(:,1:timewinidx/2)).^2,1))
axis tight, axis off

subplot(5,2,3)
hann = .5*(1-cos(2*pi*(1:timewinidx)/(timewinidx-1)));
plot(hann)
axis tight, axis off

subplot(525)
plot(hann.*d)
axis tight, axis off

subplot(526)
ff=fft(hann.*d);
plot(mean(abs(ff(1:timewinidx/2)).^2,1))
axis tight, axis off

%% Figure 16.2

channel2plot    = 'p7';
frequency2plot  = 15;  % in Hz
timepoint2plot  = 200; % ms

nw_product      = 3;  % determines the frequency smoothing, given a specified time window
times2save      = -300:50:1000;
baseline_range  = [-200 -00];
timewin         = 400; % in ms

% convert time points to indices
times2saveidx = dsearchn(EEG.times',times2save'); 
timewinidx    = round(timewin/(1000/EEG.srate));

% find baselinetimepoints
baseidx = zeros(size(baseline_range));
[~,baseidx(1)] = min(abs(times2save-baseline_range(1)));
[~,baseidx(2)] = min(abs(times2save-baseline_range(2)));
% note that the following line is equivalent to the previous three
%baseidx = dsearchn(times2save',baseline_range');

% define tapers
tapers = dpss(timewinidx,nw_product); % note that in practice, you'll want to set the temporal resolution to be a function of frequency
% define frequencies for FFT
f = linspace(0,EEG.srate/2,floor(timewinidx/2)+1);

% find logical channel index
chanidx = strcmpi(channel2plot,{EEG.chanlocs.labels});

% initialize output matrix
multitaper_tf = zeros(floor(timewinidx/2)+1,length(times2save));

% loop through time bins
for ti=1:length(times2saveidx)
    
    % initialize power vector (over tapers)
    taperpow = zeros(floor(timewinidx/2)+1,1);
    
    % loop through tapers
    for tapi = 1:size(tapers,2)-1
        
        % window and taper data, and get power spectrum
        data      = bsxfun(@times,squeeze(EEG.data(chanidx,times2saveidx(ti)-floor(timewinidx/2)+1:times2saveidx(ti)+ceil(timewinidx/2),:)),tapers(:,tapi));
        pow       = fft(data,timewinidx)/timewinidx;
        pow       = pow(1:floor(timewinidx/2)+1,:);
        taperpow  = taperpow + mean(pow.*conj(pow),2);
    end
    
    % finally, get power from closest frequency
    multitaper_tf(:,ti) = taperpow/tapi;
end

% db-correct
db_multitaper_tf = 10*log10( multitaper_tf ./ repmat(mean(multitaper_tf(:,baseidx(1):baseidx(2)),2),1,length(times2save)) );


% plot time courses at one frequency band
figure
subplot(121)
[junk,freq2plotidx]=min(abs(f-frequency2plot)); % can replace "junk" with "~"
plot(times2save,mean(log10(multitaper_tf(freq2plotidx-2:freq2plotidx+2,:)),1))
title([ 'Sensor ' channel2plot ', ' num2str(frequency2plot) ' Hz' ])
axis square
set(gca,'xlim',[times2save(1) times2save(end)])

subplot(122)
[junk,time2plotidx]=min(abs(times2save-timepoint2plot));
plot(f,log10(multitaper_tf(:,time2plotidx)))
title([ 'Sensor ' channel2plot ', ' num2str(timepoint2plot) ' ms' ])
axis square
set(gca,'xlim',[f(1) 40])


% plot full TF map
figure
contourf(times2save,f,db_multitaper_tf,40,'linecolor','none')
set(gca,'clim',[-2 2])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Power via multitaper from channel ' channel2plot ])

%% end.
