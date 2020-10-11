%% Analyzing Neural Time Series Data
% Matlab code for Chapter 26
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 26.1

% load sample EEG dataset
load sampleEEGdata

% names of the channels you want to synchronize
channel1 = 'p1';
channel2 = 'pz';

% create complex Morlet wavelet
center_freq = 5; % in Hz
time        = -1:1/EEG.srate:1; % time for wavelet
wavelet     = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2))/center_freq;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts;
n_convolution = n_wavelet+n_data-1;

% FFT of wavelet
fft_wavelet = fft(wavelet,n_convolution);

% initialize output time-frequency data
phase_data = zeros(2,EEG.pnts);
real_data  = zeros(2,EEG.pnts);

% find channel indices
chanidx = zeros(1,2); % always initialize!
chanidx(1) = find(strcmpi(channel1,{EEG.chanlocs.labels}));
chanidx(2) = find(strcmpi(channel2,{EEG.chanlocs.labels}));


% run convolution and extract filtered signal (real part) and phase
for chani=1:2
    fft_data = fft(squeeze(EEG.data(chanidx(chani),:,1)),n_convolution);
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(4/(2*pi*center_freq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
 
    % collect real and phase data
    phase_data(chani,:) = angle(convolution_result_fft);
    real_data(chani,:)  = real(convolution_result_fft);
end


% open and name figure
figure, set(gcf,'Name','Movie magic minimizes the mystery.','Number','off');

% draw the filtered signals
subplot(321)
filterplotH1 = plot(EEG.times(1),real_data(1,1),'b');
hold on
filterplotH2 = plot(EEG.times(1),real_data(2,1),'m');
set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[min(real_data(:)) max(real_data(:))])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(center_freq) ' Hz' ])

% draw the phase angle time series
subplot(322)
phaseanglesH1 = plot(EEG.times(1),phase_data(1,1),'b');
hold on
phaseanglesH2 = plot(EEG.times(1),phase_data(2,1),'m');
set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[-pi pi]*1.1,'ytick',-pi:pi/2:pi)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angle differences in cartesian space
subplot(323)
filterplotDiffH1 = plot(EEG.times(1),real_data(1,1)-real_data(2,1),'b');
set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[-10 10])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(center_freq) ' Hz' ])

% draw the phase angle time series
subplot(324)
phaseanglesDiffH1 = plot(EEG.times(1),phase_data(1,1)-phase_data(2,1),'b');
set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[-pi pi]*2.2,'ytick',-2*pi:pi/2:pi*2)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angles in polar space
subplot(325)
polar2chanH1 = polar([phase_data(1,1) phase_data(1,1)]',repmat([0 1],1,1)','b');
hold on
polar2chanH2 = polar([phase_data(1,1) phase_data(2,1)]',repmat([0 1],1,1)','m');
title('Phase angles from two channels')
 
% draw phase angle differences in polar space
subplot(326)
polarAngleDiffH = polar([zeros(1,1) phase_data(2,1)-phase_data(1,1)]',repmat([0 1],1,1)','k');
title('Phase angle differences from two channels')
 
% now update plots at each timestep
% Note: in/decrease skipping by 10 to speed up/down the movie
for ti=1:10:EEG.pnts
    
    % update filtered signals
    set(filterplotH1,'XData',EEG.times(1:ti),'YData',real_data(1,1:ti))
    set(filterplotH2,'XData',EEG.times(1:ti),'YData',real_data(2,1:ti))
    
    % update cartesian plot of phase angles
    set(phaseanglesH1,'XData',EEG.times(1:ti),'YData',phase_data(1,1:ti))
    set(phaseanglesH2,'XData',EEG.times(1:ti),'YData',phase_data(2,1:ti))
    
    % update cartesian plot of phase angles differences
    set(phaseanglesDiffH1,'XData',EEG.times(1:ti),'YData',phase_data(1,1:ti)-phase_data(2,1:ti))
    set(filterplotDiffH1,'XData',EEG.times(1:ti),'YData',real_data(1,1:ti)-real_data(2,1:ti))
    
    subplot(325)
    cla
    polar(repmat(phase_data(1,1:ti),1,2)',repmat([0 1],1,ti)','b');
    hold on
    polar(repmat(phase_data(2,1:ti),1,2)',repmat([0 1],1,ti)','m');
    
    subplot(326)
    cla
    polar(repmat(phase_data(2,1:ti)-phase_data(1,1:ti),1,2)',repmat([0 1],1,ti)','k');
    
    drawnow
end

%% Figure 26.2

figure
subplot(221)
polar(repmat(phase_data(2,:)-phase_data(1,:),1,2)',repmat([0 1],1,EEG.pnts)','k');
title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(phase_data,1)))))) ])

new_phase_data = phase_data;
for i=2:4
    subplot(2,2,i)
    
    % add random phase offset
    new_phase_data(1,:) = new_phase_data(1,:)+rand*pi;
    
    % plot again
    polar(repmat(new_phase_data(2,:)-new_phase_data(1,:)+pi/2,1,2)',repmat([0 1],1,EEG.pnts)','k');
    title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(new_phase_data,1)))))) ])
end

%% Figure 26.3

% note: see commented line "time_window_idx..." below for panels C and D

channel1 = 'fz';
channel2 = 'o1';

freqs2use  = logspace(log10(4),log10(30),15); % 4-30 Hz in 15 steps
times2save = -400:20:800;
timewindow = linspace(1.5,3,length(freqs2use)); % number of cycles on either end of the center point (1.5 means a total of 3 cycles))
baselinetm = [-400 -200];

% wavelet and FFT parameters
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% time in indices
times2saveidx = dsearchn(EEG.times',times2save');
baselineidx   = dsearchn(times2save',baselinetm');

chanidx    = zeros(1,2); % always initialize!
chanidx(1) = find(strcmpi(channel1,{EEG.chanlocs.labels}));
chanidx(2) = find(strcmpi(channel2,{EEG.chanlocs.labels}));

% initialize
ispc = zeros(length(freqs2use),length(times2save));
ps   = zeros(length(freqs2use),length(times2save));

% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);


for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    phase_sig1 = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    phase_sig2 = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));
    
    % phase angle differences
    phase_diffs = phase_sig1-phase_sig2;
    
    % compute ICPS over trials
    ps(fi,:) = abs(mean(exp(1i*phase_diffs(times2saveidx,:)),2));
    
    % compute time window in indices for this frequency
    time_window_idx = round((1000/freqs2use(fi))*timewindow(fi)/(1000/EEG.srate));
%     time_window_idx = round(300/(1000/EEG.srate)); % set 300 to 100 for figure 3c/d
    
    for ti=1:length(times2save)
        
        % compute phase synchronization
        phasesynch = abs(mean(exp(1i*phase_diffs(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:)),1));
        
        % average over trials
        ispc(fi,ti) = mean(phasesynch);
    end
end % end frequency loop

figure
contourf(times2save,freqs2use,ispc-repmat(mean(ispc(:,baselineidx(1):baselineidx(2)),2),1,size(ispc,2)),20,'linecolor','none')
set(gca,'clim',[-.08 .08],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

figure
plot(freqs2use,(1000./freqs2use).*timewindow*2,'o-','markerface','k')
hold on
plot(freqs2use,(1000./freqs2use).*timewindow(1)*2,'ro-','markerface','m')
ylabel('Window width (ms)'), xlabel('Frequency (Hz)')
legend({'variable windows';'fixed 3*f window'})

%% Figure 26.4

figure
for i=1:8
    subplot(8,1,i)
    plot(phase_sig1(1:200,i)-phase_sig2(1:200,i))
end

figure
subplot(121)
polar(repmat(phase_sig1(1:200,1)-phase_sig2(1:200,1),2,1),repmat([0 1]',200,1),'k');
title('Phase angle differences over time')

subplot(122)
polar(repmat(phase_sig1(100,1:i)-phase_sig2(100,1:i),2,1),repmat([0 1]',1,i),'k');
title('Phase angle differences over trials')

%% Figure 26.5

figure
contourf(times2save,freqs2use,bsxfun(@minus,ps,mean(ps(:,baselineidx(1):baselineidx(2)),2)),20,'linecolor','none')
set(gca,'clim',[-.2 .2],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[-300 800])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%% figure 26.6

time2use = 300; % ms
niterations = 50; % you can decrease this to make the code a bit faster

% initialize
ispcByNandF = zeros(length(freqs2use),EEG.trials);
time2useidx = dsearchn(times2save',time2use);

% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);

for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    phase_sig1 = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    phase_sig2 = angle(reshape(convolution_result_fft,EEG.pnts,EEG.trials));
    
    % phase angle differences
    phase_diffs = phase_sig1-phase_sig2;
    
    for n=1:EEG.trials
        % multiple iterations to select different random sets of trials
        for iteri=1:niterations
            trials2use = randsample(EEG.trials,n);
            ispcByNandF(fi,n) = ispcByNandF(fi,n) + mean(abs(mean(exp(1i*phase_diffs(times2saveidx(time2useidx)-time_window_idx:times2saveidx(time2useidx)+time_window_idx,trials2use)),2)),1);
        end
    end
end

figure
plot(1:EEG.trials,ispcByNandF/iteri)
xlabel('Trials')
ylabel('ICPS-trials')

%% Figure 26.7

% initialize
data4test  = zeros(2,EEG.pnts,EEG.trials);
data4power = zeros(2,EEG.pnts,EEG.trials);

amp_mod = 0.00001;

for triali=1:EEG.trials
    % each trial is a random channel and trial
    trialdata1 = EEG.data(chanidx(1),:,triali);
    trialdata2 = EEG.data(chanidx(2),:,triali);
    
    % Un/comment the next line of code for band-pass filtered data.
    % This uses the eegfilt function, which is part of the eeglab toolbox.
    % You can also replace this function with your preferred filter method (chapter 14).
    trialdata1 = eegfilt(double(trialdata1),EEG.srate,10,20);
    trialdata2 = eegfilt(double(trialdata2),EEG.srate,10,20);
    
    % phase angle differences, with and without amplitude dampening
    data4test(1,:,triali) = angle(hilbert(trialdata1)) - angle(hilbert(trialdata2));
    data4test(2,:,triali) = angle(hilbert(trialdata1)) - angle(hilbert(trialdata2*amp_mod));
    
    data4power(1,:,triali) = abs(hilbert(trialdata2)).^2;
    data4power(2,:,triali) = abs(hilbert(trialdata2*amp_mod)).^2;
end

% compute ITPC
ispc_nomod = abs(mean(exp(1i*data4test(1,:,:)),3));
ispc_mod   = abs(mean(exp(1i*data4test(2,:,:)),3));

% compute power
power = squeeze(mean(data4power,3));

% plot!
figure
subplot(311)
plot(EEG.times,trialdata2)
hold on
plot(EEG.times,trialdata2*amp_mod,'r')
title('Amplitude modulator')

subplot(312)
plot(EEG.times,data4test(1,:,10))
hold on
plot(EEG.times,data4test(2,:,10),'r')
axis tight
set(gca,'ytick',-2*pi:pi:2*pi)
title('Example trials')

subplot(313)
plot(EEG.times,ispc_mod,'ro')
hold on
h=plotyy(EEG.times,ispc_nomod,EEG.times,squeeze(mean(data4power(1,:,:),3)));
legend({'amplitude modulation';'no amp mod'})
set(h(2),'ylim',[12 36])
set(h(1),'ylim',[0 .4])
xlabel('Time (ms)'), ylabel('ICPS')
title('ICPS')

figure
subplot(121)
plot(power(1,:),ispc_nomod,'.')
axis square
xlabel('Power'), ylabel('ICPS')
title('non-modulated power')

subplot(122)
plot(power(2,:),ispc_mod,'.')
axis square
xlabel('Power'), ylabel('ICPS')
title('modulated power')

%% figure 26.8

% select channels
channel1 = 'fz';
channel2 = 'o1';

% wavelet and FFT parameters
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

chanidx    = zeros(1,2); % always initialize!
chanidx(1) = find(strcmpi(channel1,{EEG.chanlocs.labels}));
chanidx(2) = find(strcmpi(channel2,{EEG.chanlocs.labels}));

% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);


% initialize
spectcoher = zeros(length(freqs2use),length(times2save));

for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % compute power and cross-spectral power
    spec1 = mean(sig1.*conj(sig1),2);
    spec2 = mean(sig2.*conj(sig2),2);
    specX = abs(mean(sig1.*conj(sig2),2)).^2;
    
    % alternative notation for the same procedure, using the Euler-like expression: Me^ik
    %spec1 = mean(abs(sig1).^2,2);
    %spec2 = mean(abs(sig2).^2,2);
    %specX = abs(mean( abs(sig1).*abs(sig2) .* exp(1i*(angle(sig1)-angle(sig2))) ,2)).^2;
    
    % compute spectral coherence, using only requested time points
    spectcoher(fi,:) = specX(times2saveidx)./(spec1(times2saveidx).*spec2(times2saveidx));
    
    % yet another equivalent notation, just FYI
    %spec1 = sum(sig1.*conj(sig1),2);
    %spec2 = sum(sig2.*conj(sig2),2);
    %specX = sum(sig1.*conj(sig2),2);
    %spectcoher(fi,:) = abs(specX(times2saveidx)./sqrt(spec1(times2saveidx).*spec2(times2saveidx))).^2;
    
    
    % imaginary coherence
    %spec1 = sum(sig1.*conj(sig1),2);
    %spec2 = sum(sig2.*conj(sig2),2);
    %specX = sum(sig1.*conj(sig2),2);
    % spectcoher(fi,:) = abs(imag(specX(times2saveidx)./sqrt(spec1(times2saveidx).*spec2(times2saveidx))));
end


figure
subplot(121)
contourf(times2save,freqs2use,spectcoher,20,'linecolor','none') % 
set(gca,'clim',[0 .2],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[times2save(1) times2save(end)])
title('"Raw" spectral coherence')

subplot(122)
contourf(times2save,freqs2use,spectcoher-repmat(mean(spectcoher(:,baselineidx(1):baselineidx(2)),2),1,size(spectcoher,2)),20,'linecolor','none') % 
set(gca,'clim',[-.1 .1],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)),'xlim',[times2save(1) times2save(end)])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Baseline-subtracted spectral coherence')

%% Figure 26.9

% number of "trials"
n = 100;

figure

subplot(221)
phases = rand(n,1)*pi;
polar([phases; phases],repmat([0 1]',n,1),'k');
pli  = abs(mean(sign(imag(exp(1i*phases)))));
ispc = abs(mean(exp(1i*phases)));
title([ 'PLI=' num2str(pli) ', ISPC=' num2str(ispc) ])

subplot(222)
phases = phases-pi/2;
polar([phases; phases],repmat([0 1]',n,1),'k');
pli  = abs(mean(sign(imag(exp(1i*phases)))));
ispc = abs(mean(exp(1i*phases)));
title([ 'PLI=' num2str(pli) ', ISPC=' num2str(ispc) ])


subplot(223)
phases = rand(n,1)/2+pi/3+.25;
polar([phases; phases],repmat([0 1]',n,1),'k');
pli  = abs(mean(sign(imag(exp(1i*phases)))));
ispc = abs(mean(exp(1i*phases)));
title([ 'PLI=' num2str(pli) ', ISPC=' num2str(ispc) ])

subplot(224)
phases = phases-pi/2;
polar([phases; phases],repmat([0 1]',n,1),'k');
pli  = abs(mean(sign(imag(exp(1i*phases)))));
ispc = abs(mean(exp(1i*phases)));
title([ 'PLI=' num2str(pli) ', ISPC=' num2str(ispc) ])


%% Figure 26.10

% select channels
channel1 = 'fz';
channel2 = 'o1';

% specify some time-frequency parameters
freqs2use  = logspace(log10(4),log10(30),15); % 4-30 Hz in 15 steps
times2save = -400:10:800;
timewindow = linspace(1.5,3,length(freqs2use)); % number of cycles on either end of the center point (1.5 means a total of 3 cycles))
baselinetm = [-400 -200];

% wavelet and FFT parameters
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% time in indices
times2saveidx = dsearchn(EEG.times',times2save');
baselineidxF  = dsearchn(EEG.times',baselinetm');  % for the full temporal resolution data (thanks to Daniel Roberts for finding/reporting this bug here!)
baselineidx   = dsearchn(times2save',baselinetm'); % for the temporally downsampled data

chanidx = zeros(1,2); % always initialize!
chanidx(1) = find(strcmpi(channel1,{EEG.chanlocs.labels}));
chanidx(2) = find(strcmpi(channel2,{EEG.chanlocs.labels}));

% data FFTs
data_fft1 = fft(reshape(EEG.data(chanidx(1),:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);

% initialize
ispc    = zeros(length(freqs2use),EEG.pnts);
pli     = zeros(length(freqs2use),EEG.pnts);
wpli    = zeros(length(freqs2use),EEG.pnts);
dwpli   = zeros(length(freqs2use),EEG.pnts);
dwpli_t = zeros(length(freqs2use),length(times2save));
ispc_t  = zeros(length(freqs2use),length(times2save));

for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % cross-spectral density
    cdd = sig1 .* conj(sig2);
    
    % ISPC
    ispc(fi,:) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to ispc(fi,:) = abs(mean(exp(1i*(angle(sig1)-angle(sig2))),2));
    
    
    % take imaginary part of signal only
    cdi = imag(cdd);
    
    % phase-lag index
    pli(fi,:)  = abs(mean(sign(imag(cdd)),2));
    
    % weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
    wpli(fi,:) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);
    
    % debiased weighted phase-lag index (shortcut, as implemented in fieldtrip)
    imagsum      = sum(cdi,2);
    imagsumW     = sum(abs(cdi),2);
    debiasfactor = sum(cdi.^2,2);
    dwpli(fi,:)  = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
    
    % compute time window in indices for this frequency
    time_window_idx = round((1000/freqs2use(fi))*timewindow(fi)/(1000/EEG.srate));

    for ti=1:length(times2save)
        imagsum        = sum(cdi(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:),1);
        imagsumW       = sum(abs(cdi(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:)),1);
        debiasfactor   = sum(cdi(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:).^2,1);
        dwpli_t(fi,ti) = mean((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor));

        % compute phase synchronization
        phasesynch     = abs(mean(exp(1i*angle(cdd(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:))),1));
        ispc_t(fi,ti)  = mean(phasesynch);
    end
end

% baseline subtraction from all measures
ispc    = bsxfun(@minus,ispc,mean(ispc(:,baselineidxF(1):baselineidxF(2)),2)); % not plotted in the book, but you can plot it for comparison with PLI
ispc_t  = bsxfun(@minus,ispc_t,mean(ispc_t(:,baselineidx(1):baselineidx(2)),2));
pli     = bsxfun(@minus,pli,mean(pli(:,baselineidxF(1):baselineidxF(2)),2));
dwpli   = bsxfun(@minus,dwpli,mean(dwpli(:,baselineidxF(1):baselineidxF(2)),2));
dwpli_t = bsxfun(@minus,dwpli_t,mean(dwpli_t(:,baselineidx(1):baselineidx(2)),2));

figure
subplot(221)
contourf(times2save,freqs2use,pli(:,times2saveidx),40,'linecolor','none')
set(gca,'clim',[-.3 .3],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('PLI over trials')

subplot(222)
contourf(times2save,freqs2use,dwpli(:,times2saveidx),40,'linecolor','none')
set(gca,'clim',[-.2 .2],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('dWPLI over trials')

subplot(223)
contourf(times2save,freqs2use,ispc_t,40,'linecolor','none')
set(gca,'clim',[-.1 .1],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('ICPS over time')

subplot(224)
contourf(times2save,freqs2use,dwpli_t,40,'linecolor','none')
set(gca,'clim',[-.1 .1],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
title('dWPLI over time')

%% Figure 26.11

trial2plot = 10; % any trial between 1 and 99 (book uses trial 10)
center_freq = 4.6; % Hz (book uses 4.6)


% create wavelet and take FFT
s = 4.5/(2*pi*center_freq);
wavelet_fft = fft( exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
% phase angles from channel 1 via convolution
convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
sig1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
% phase angles from channel 2 via convolution
convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
sig2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
% cross-spectral density
xsd  = sig1 .* conj(sig2);
xsdi = imag(xsd);

dwpli = zeros(size(EEG.times));
ispc  = zeros(size(EEG.times));

[junk,animate_start] = min(abs(EEG.times-0));
[junk,animate_stop]  = min(abs(EEG.times-1000));

time_window_idx = round(100*timewindow(1)/(1000/EEG.srate));

figure
subplot(121)
hpol = polar(repmat(angle(xsd(animate_start:animate_start+time_window_idx-1,trial2plot)),1,2)',[zeros(time_window_idx,1) ones(time_window_idx,1)]','k-o');
subplot(122)
hplo2 = plot(EEG.times,0,'r');
hold on
hplo1 = plot(EEG.times,0,'b');

for idx=animate_start:animate_stop
    
    % update angles
    for i=1:length(hpol)
        set(hpol(i),'XData',[0; cos(angle(xsd(idx+i,trial2plot)))],'YData',[0; sin(angle(xsd(idx+i,trial2plot)))]);
    end
    title([ num2str(round(EEG.times(idx))) '-' num2str(round(EEG.times(idx+i))) ' ms' ])

    % compute ICPS and dwPLI
    ispc(idx) = abs(mean(exp(1i*angle(xsd(idx:idx+i,trial2plot))),1));
    
    imagsum        = sum(xsdi(idx:idx+i,trial2plot),1);
    imagsumW       = sum(abs(xsdi(idx:idx+i,trial2plot)),1);
    debiasfactor   = sum(xsdi(idx:idx+i,trial2plot).^2,1);
    dwpli(idx) = mean((imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor));
    
    
    set(hplo1,'XData',EEG.times(1:idx),'YData',ispc(1:idx));
    set(hplo2,'XData',EEG.times(1:idx),'YData',dwpli(1:idx));
    set(gca,'xlim',EEG.times([animate_start animate_stop]),'ylim',[-.1 1.1])
    axis square
    
    pause(0.01)
end

%% Figure 26.12

% This figure is generated in the code below.

%% Figure 26.13

% generate inline functions (small functions that you define without being
% saved as general Matlab functions)
vtest  = inline('n.*icpcmag*cos(val).*sqrt(2./n)');
gvtest = inline('n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)))');

% Since initially writing this code, Matlab decided to make the inline
% function obsolete in future versions. The following two lines produce
% the identical 'anonymous functions' as the previous two lines.
%vtest  = @(icpcmag,n,val) n.*icpcmag*cos(val).*sqrt(2./n);
%gvtest = @(icpcmag,n,val) n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)));


% figure
clf
subplot(221)
n = 2:100;
plot(n,1-normcdf(vtest(.3,n,pi/10)))
hold on
plot(n,1-normcdf(gvtest(.3,n,pi/10)),'m')
legend({'v-test';'gv-test'})
xlabel('Number of points'), ylabel('P-value')
set(gca,'ylim',[0 .6])
title('angle = pi/10')

subplot(222)
plot(n,1-normcdf(vtest(.3,n,pi/3)))
hold on
plot(n,1-normcdf(gvtest(.3,n,pi/3)),'m')
legend({'v-test';'gv-test'})
xlabel('Number of points'), ylabel('P-value')
set(gca,'ylim',[0 .6])
title('angle = pi/3')

subplot(223)
x=linspace(-pi,pi,50);
n=15;
polar(x,1-normcdf(vtest(.3,n,x-0)))
hold on
polar(x,1-normcdf(gvtest(.3,n,x-0)),'m')
title([ 'N=' num2str(n) ])

subplot(224)
n=600;
polar(x,1-normcdf(vtest(.3,n,x-0)))
hold on
polar(x,1-normcdf(gvtest(.3,n,x-0)),'m')
set(gca,'xtick',round((-pi:pi/4:pi)*100)/100)
title([ 'N=' num2str(n) ])

% number of simulated data points
numUsims = 10000;

u = zeros(2,numUsims);

for i=1:numUsims
    
    % make some noise
    fake_phase_data = rand(2,EEG.pnts)*pi*2-pi;
    
    % compute ispc
    ispc_mag = abs  (mean(exp(1i*(diff(fake_phase_data,1)))));
    ispc_phs = angle(mean(exp(1i*(diff(fake_phase_data,1)))));
    
    % compute statistics
    u(1,i) = vtest (ispc_mag,EEG.pnts,ispc_phs-0);
    u(2,i) = gvtest(ispc_mag,EEG.pnts,ispc_phs-0);
end


% This figure is also figure 26.12 but with no log-scaling
figure
for i=1:2
    subplot(1,2,i)
    
    [y,x]=hist(u(i,:),100);
    h=bar(x,log10(y),'histc');
    set(h,'linestyle','none')
    
    title([ num2str(100*sum((1-normcdf(u(i,:)))<.05)/length(u)) '% false positive' ])
end


nrange    = 10:300; % here n is number of datapoints
ispcrange = .05:.01:.7;
pvalmat   = zeros(2,length(nrange),length(ispcrange));

for ni=1:length(nrange)
    for mi=1:length(ispcrange)
        
        n = nrange(ni);
        
        pvalmat(1,ni,mi) = 1-normcdf( vtest(ispcrange(mi),nrange(ni),pi/5));
        pvalmat(2,ni,mi) = 1-normcdf(gvtest(ispcrange(mi),nrange(ni),pi/5));
    end
end

figure
subplot(121)
imagesc(ispcrange,nrange,squeeze(pvalmat(1,:,:))), axis xy, 
set(gca,'clim',[0 .5])
xlabel('ICPS strength'), ylabel('N')
title('v-test')

subplot(122)
imagesc(ispcrange,nrange,squeeze(pvalmat(2,:,:))), axis xy, 
title('gv-test')
xlabel('ICPS strength'), ylabel('N')
set(gca,'clim',[0 .5])

%% end.
