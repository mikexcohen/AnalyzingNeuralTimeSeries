%% Analyzing Neural Time Series Data
% Matlab code for Chapter 19
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 


%% figure 19.1

% define angles
a = [0.2 2*pi-.2];

figure
% plot unit vectors defined by those angles
polar([0 0; a],[0 0; 1 1])
hold on
% plot a unit vector with the average angle
polar([0 mean(a)],[0 1],'r')
% plot the average vector
polar([0 angle(mean(exp(1i*a)))],[0 abs(mean(exp(1i*a)))],'m')

%% figure 19.2

% load sample scalp EEG data
load sampleEEGdata

% center frequency
centerfreq = 12; % in Hz
chan2plot  = 'pz';
times2plot = [ 200 800 ]; % in ms, from stimulus onset

% definte convolution parameters
n_wavelet     = EEG.pnts;
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));

% create wavelet. Note that the time vector to create the wavelet uses the
% same number of points as there are in the EEG data. If the EEG data has
% an even number of points (which is the case here--640), the wavelet will
% not have an exact center point. Though not technically incorrect, using
% a wavelet with an even number of points should be avoided when possible.
time    = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;

% get FFT of data
eegfft = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);

% convolution
eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
eegconv = eegconv(1:n_convolution);
% Cut the edges off the result of convolution. Because the wavelet here was
% created with the same number of time points as the EEG data, EEG.pnts
% also corresponds to the number of time points in the wavelet.
eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);

% plot
figure
for subploti=1:2
    subplot(2,2,subploti)
    [junk,idx]=min(abs(EEG.times-times2plot(subploti)));
    polar([zeros(1,EEG.trials); angle(eegconv(idx,:))],[zeros(1,EEG.trials); ones(1,EEG.trials)])
    title([ 'ITPC at ' num2str(times2plot(subploti)) ' ms = ' num2str(round(1000*abs(mean(exp(1i*angle(eegconv(idx,:))))))/1000) ])
    
    subplot(2,2,subploti+2)
    hist(angle(eegconv(idx,:)),20)
    set(gca,'xlim',[-pi pi]*1.1,'xtick',-pi:pi/2:pi,'xticklabel',{'-pi';'-pi/2';'0';'pi';'pi/2'},'ylim',[0 20])
end

%% Figure 19.3

vectors{1} = [0 pi/3];
vectors{2} = [0 pi/2];
vectors{3} = [0 2*pi/3];
vectors{4} = [0 pi*.9];

figure
for i=1:length(vectors)
    subplot(1,length(vectors),i)
    
    % plot individual unit vectors
    polar([0 vectors{i}(1)],[0 1],'k')
    hold on
    polar([0 vectors{i}(2)],[0 1],'k')
    
    % plot mean vector (ignore the math for now, 
    % this will be discussed later in the chapter)
    meanvect = mean(exp(1i*vectors{i}));
    polar([0 angle(meanvect)],[0 abs(meanvect)],'r')
    title([ 'Mean vector length: ' num2str(abs(meanvect)) ])
end

%% Figure 19.4

% get FFT of data
chan2plot  = 'pz';
centerfreq = 12; % in Hz
eegfft = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);

% ITPC at one frequency band
wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;
% convolution
eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
eegconv = eegconv(1:n_convolution);
eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);

figure
plot(EEG.times,abs(mean(exp(1i*angle(eegconv)),2)))
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)')
ylabel('ITPC')


% TF plot of ITPC
frequencies = logspace(log10(4),log10(40),20);
s = logspace(log10(3),log10(10),length(frequencies))./(2*pi*frequencies);
itpc = zeros(length(frequencies),EEG.pnts);

for fi=1:length(frequencies)
    % create wavelet
    wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)))/frequencies(fi);
    
    % convolution
    eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);
    
    % extract ITPC
    itpc(fi,:) = abs(mean(exp(1i*angle(eegconv)),2));
end

figure
contourf(EEG.times,frequencies,itpc,40,'linecolor','none')
set(gca,'clim',[0 .6],'xlim',[-200 1000])
xlabel('Time (ms)'), ylabel('Frequencies (Hz)')

%% figure 19.5

n_trials    = 500;
itpcByNFake = zeros(1,n_trials);

for n=1:n_trials
    for iteri=1:50
        itpcByNFake(n) = itpcByNFake(n) + abs(mean(exp(1i*(rand(1,n)*2*pi-pi))));
    end
end
itpcByNFake = itpcByNFake./iteri;

% Z and P (you will learn more about these formulae in chapter 34)
itpcByNFakeZ = (1:n_trials).*(itpcByNFake.^2);
itpcByNFakeP = exp(sqrt(1+4.*(1:n_trials)+4*( ((1:n_trials).^2) - ((1:n_trials).*itpcByNFake).^2))-(1+2*(1:n_trials))); % used later
itpcFakeCrit = sqrt( -log(.01)./(1:n_trials) );

figure
plot(1:n_trials,itpcByNFake)
hold on
plot(1:n_trials,itpcFakeCrit,'r')
set(gca,'ylim',[0 1])
legend({'ITPC_r_a_w';'ITPC_C_r_i_t'})

%% figure 19.6

% center frequency
centerfreq = 6; % Hz
chan2plot  = 'fcz';

% definte convolution parameters
n_wavelet     = EEG.pnts;
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));

% create wavelet
time    = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;

% get FFT of data
eegfft  = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);

% convolution
eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
eegconv = eegconv(1:n_convolution);
eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);

figure

% compute ITPC as function of N
itpcByN = zeros(1,EEG.trials);
for n=1:EEG.trials
    % multiple iterations to select different sets of random trials
    for iteri=1:50
        trials2use  = randsample(EEG.trials,n);
        itpcByN(n)  = itpcByN(n) + mean(abs(mean(exp(1i*angle(eegconv(282:372,trials2use))),2)),1);
    end
end
itpcByN = itpcByN./iteri;

% Z and P
itpcByNZ = (1:EEG.trials).*(itpcByN.^2);
itpcByNP = exp(sqrt(1+4.*(1:EEG.trials)+4*( ((1:EEG.trials).^2) - ((1:EEG.trials).*itpcByN).^2))-(1+2*(1:EEG.trials)));
itpcCrit = sqrt( -log(.01)./(1:EEG.trials) ); % .01 is the p-value in this case

subplot(211)
plot(itpcByN)
hold on
plot(itpcCrit,'r')
xlabel('Number of trials in analysis')
ylabel('ITPC')
set(gca,'ylim',[.2 1])


% (Note that this figure will look different than that in the book
% because trials are randomly selected.)
subplot(212)
plot(EEG.times,abs(mean(exp(1i*angle(eegconv)),2)))
hold on
randomtrials2plot = randperm(EEG.trials); % could also use randsamp if you have the stats toolbox
plot(EEG.times,abs(mean(exp(1i*angle(eegconv(:,randomtrials2plot(1:20)))),2)),':')
xlabel('Time (ms)'), ylabel('ITPC')
legend({'99 trials';'20 trials'})

%% Figure 34.8

% This cell concerns statistical evaluation of ITPC values. It will be
% discussed in depth in chapter 34, but the code is presented here because
% it relies on calculations from the previous two figures. 

% p-values under assumption of von Mises distribution
approx_pval_fake = exp(-itpcByNFakeZ);
approx_pval_real = exp(-itpcByNZ);

ncutoff = 10;

figure

subplot(121)
plot(itpcByNFakeP(ncutoff+1:end),approx_pval_fake(ncutoff+1:end),'mo','markersize',8)
hold on
plot(itpcByNFakeP(1:ncutoff),approx_pval_fake(1:ncutoff),'k.')
plot([0 1],[0 1],'k')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.25:1,'ytick',0:.25:1)
axis square
ylabel('approximate p')
xlabel('exact p')
title('P-values from fake data')

subplot(122)
plot(itpcByNP(ncutoff+1:end),approx_pval_real(ncutoff+1:end),'mo','markersize',8)
hold on
plot(itpcByNP(1:ncutoff),approx_pval_real(1:ncutoff),'k.')
plot([0.00001 1],[0.00001 1],'k')
set(gca,'xlim',[.0001 1],'ylim',[.0001 1],'xscale','log','yscale','log')
axis square
ylabel('approximate p')
xlabel('exact p')
title('P-values from real data')

%% figure 19.7
% (This figure takes a while to generate. You could also reduce the number
% of iterations to speed it up.)

% Data for this cell are from figure 19.6. If you want to plot 
% results from a different channel, change the electrode in 19.6.

% center frequencies
frequencies = 1:40;

niterations = 50;

% initialize
itpcByNandF = zeros(length(frequencies),EEG.trials);

for fi=1:length(frequencies)
    
    centerfreq = frequencies(fi);
    
    wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;
    
    % convolution
    eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);
    
    for n=1:EEG.trials
        % multiple iterations to select different random sets of trials
        for iteri=1:niterations
            trials2use  = randsample(EEG.trials,n);
            itpcByNandF(fi,n)  = itpcByNandF(fi,n) + mean(abs(mean(exp(1i*angle(eegconv(282:372,trials2use))),2)),1);
        end
    end
end

figure
contourf(1:EEG.trials,frequencies,itpcByNandF./iteri,40,'linecolor','none')
colormap gray
set(gca,'clim',[.1 .5])
xlabel('Trials')
ylabel('Frequency (Hz)')

%% Figure 19.8

figure, set(gcf,'name','Rayleigh''s Z')
subplot(121)
plot(itpcByN)
hold on
plot(itpcByNFake(1:99),'r')
xlabel('Trial count'), ylabel('ITPC')
axis square

subplot(122)
plot(itpcByNZ)
hold on
plot(itpcByNFakeZ(1:99),'r')
xlabel('Trial count'), ylabel('ITPC_Z')
legend({'real data';'fake data'})
axis square

%% Figure 19.9

% pick electrode and frequencies
chan2plot   = 'fcz';
frequencies = 1:40;

% definte convolution parameters
time          = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
n_wavelet     = EEG.pnts;
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));

plotlegends={'without jitter';'with jitter'};

figure

baselinetime = [ -300 -100 ];
baseidx=dsearchn(EEG.times',baselinetime(1)):dsearchn(EEG.times',baselinetime(2));

for simuli=1:2
    
    % add time jitter (or not)
    tempdat = squeeze(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:));
    for ti=1:size(tempdat,2)
        timejitter = ceil(rand*10)*(simuli-1); % (when simuli==1, timejitter==0)
        tempdat(:,ti) = tempdat([ timejitter+1:end 1:timejitter ],ti);
    end
    
    % get FFT of data
    eegfft = fft(reshape(tempdat,1,[]),n_conv_pow2);
    
    % initialize
    itpc = zeros(length(frequencies),EEG.pnts);
    powr = zeros(length(frequencies),EEG.pnts);
    
    for fi=1:length(frequencies)
        
        centerfreq = frequencies(fi);
        
        wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;
        
        % convolution
        eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
        eegconv = eegconv(1:n_convolution);
        eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);
        
        % compute and store itpc and power
        itpc(fi,:) = abs(mean(exp(1i*angle(eegconv)),2));
        powr(fi,:) = mean(abs(eegconv).^2,2);
        powr(fi,:) = 10*log10(powr(fi,:)./mean(powr(fi,baseidx),2));
    end
    
    subplot(2,2,simuli)
    contourf(EEG.times,frequencies,itpc,40,'linecolor','none')
    set(gca,'clim',[.1 .5],'xlim',[-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title([ 'ITPC ' plotlegends{simuli} ])
    
    subplot(2,2,simuli+2)
    contourf(EEG.times,frequencies,powr,40,'linecolor','none')
    set(gca,'clim',[-3 3],'xlim',[-200 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title([ 'DB-power ' plotlegends{simuli} ])
end

%% Figure 19.10

% initialize
data4test  = zeros(2,length(time),EEG.trials);
sensor2use = 'p7';

% amplitude modulation (modulate power by 1 Hz sine wave)
time    = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
amp_mod = (sin(2*pi*1.*time)+2)-1;

for triali=1:EEG.trials
    % each trial is a random channel and trial
    trialdata = EEG.data(strcmpi(sensor2use,{EEG.chanlocs.labels}),:,triali);
    
    % Uncomment the next line of code for band-pass filtered data.
    % This uses the eegfilt function, which is part of the eeglab toolbox.
    trialdata = eegfilt(double(trialdata),EEG.srate,10,20);
    
    data4test(1,:,triali) = trialdata.*amp_mod;
    data4test(2,:,triali) = trialdata;
end

% compute ITPC
itpc_mod   = abs(mean(exp(1i*angle(hilbert(data4test(1,:,:)))),3));
itpc_nomod = abs(mean(exp(1i*angle(hilbert(data4test(2,:,:)))),3));


% plot!
figure
subplot(311)
plot(EEG.times,amp_mod)
set(gca,'ylim',[-.2 2.2])
title('Amplitude modulator')

subplot(312)
plot(EEG.times,data4test(1,:,10))
hold on
plot(EEG.times,data4test(2,:,10),'r')
axis tight
title('Example trials')

subplot(313)
plot(EEG.times,itpc_mod)
hold on
plot(EEG.times,itpc_nomod,'r')
legend({'amplitude modulation';'no amp mod'})
set(gca,'ylim',[0 1])
xlabel('Time (ms)'), ylabel('ITPC')
title('ITPC')

figure
plot(amp_mod,itpc_mod,'.')
xlabel('Power modulation'), ylabel('ITPC')

%% Figure 19.11
% (This figure is populated with randomly generated data, 
%  and so will look different from the book figure.)

randvects = rand(50,1)*2*pi-pi;
vectormod = (randvects + randn(size(randvects))).^2;
vectormod2 = vectormod-min(vectormod)+1; % make sure no negative values


figure

% ITPC
subplot(221)
polar([zeros(size(randvects)) randvects]',[zeros(size(randvects)) ones(size(randvects))]','k')
title([ 'ITPC = ' num2str(abs(mean(exp(1i*randvects)))) ])

% wITPC
subplot(222)
polar([zeros(size(randvects)) randvects]',[zeros(size(randvects)) vectormod]','k')

witpc = abs(mean(vectormod.*exp(1i*randvects)));
perm_witpc = zeros(1,1000);

for i=1:1000
    perm_witpc(i) = abs(mean(vectormod(randperm(length(vectormod))).*exp(1i*randvects)));
end

witpc_z = (witpc-mean(perm_witpc))./std(perm_witpc);
title([ '_wITPC_z = ' num2str(witpc_z) ])

% example of one permutation
subplot(223)
polar([zeros(size(randvects)) randvects]',[zeros(size(randvects)) vectormod(randperm(length(vectormod)))]','k')
title('One null hypothesis iteration')

% histogram of null-hypothesis WITPC
subplot(224)
[y,x]=hist(perm_witpc,50);
h=bar(x,y,'histc');
set(h,'linestyle','none')
hold on
plot([witpc witpc],get(gca,'ylim')/2,'m')
title('Histogram of null hypothesis wITPC')

%% Figure 19.12

centerfreq  = 6;
channel2use = 'po7';
times2save  = -200:50:1200;

% initialize matrix to store RTs
rts = zeros(size(EEG.epoch));

for ei=1:EEG.trials
    
    % find which event is time=0, and take the latency of the event thereafter.
    time0event = find(cell2mat(EEG.epoch(ei).eventlatency)==0);
    
    % use try-catch in case of no response
    try
        rts(ei) = EEG.epoch(ei).eventlatency{time0event+1};
    catch me;
        rts{ei} = NaN;
    end
end


% definte convolution parameters
time          = -EEG.pnts/EEG.srate/2:1/EEG.srate:EEG.pnts/EEG.srate/2-1/EEG.srate;
n_wavelet     = EEG.pnts;
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));

% get FFT of data and wavelet
eegfft  = fft(reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,n_data),n_conv_pow2);
wavefft = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2))),n_conv_pow2);

% convolution
eegconv = ifft(wavefft.*eegfft);
eegconv = eegconv(1:n_convolution);
eegconv = reshape(eegconv(floor((EEG.pnts-1)/2):end-1-ceil((EEG.pnts-1)/2)),EEG.pnts,EEG.trials);

phase_angles = angle(eegconv);

% initialize
itpc    = zeros(size(times2save));
witpc   = zeros(size(times2save));
witpc_z = zeros(size(times2save));

for ti=1:length(times2save)
    
    % find index for this time point
    [junk,timeidx] = min(abs(EEG.times-times2save(ti)));
    
    % ITPC is unmodulated phase clustering
    itpc(ti) = abs(mean(exp(1i*phase_angles(timeidx,:))));
    
    % wITPC is rts modulating the length of phase angles
    witpc(ti) = abs(mean(rts.*exp(1i*phase_angles(timeidx,:))));
    
    % permutation testing
    perm_witpc = zeros(1,1000);
    for i=1:1000
        perm_witpc(i) = abs(mean(rts(randperm(EEG.trials)).*exp(1i*phase_angles(timeidx,:))));
    end
    
    witpc_z(ti) = (witpc(ti)-mean(perm_witpc))./std(perm_witpc);
end

figure

% plot ITPC
subplot(311)
h=plotyy(EEG.times,abs(mean(exp(1i*phase_angles),2)),times2save,witpc);
set(h,'xlim',[times2save(1) times2save(end)]);
title([ 'ITPC and wITPC at ' channel2use ])
legend({'ITPC';'wITPC'})

% plot wITPCz
subplot(312)
plot(times2save,witpc_z)
hold on
plot(get(gca,'xlim'),[0 0],'k')
set(gca,'xlim',[times2save(1) times2save(end)])
xlabel('Time (ms)'), ylabel('wITPCz')
title([ 'wITPCz at ' channel2use ])


subplot(325)
plot(itpc,witpc,'.')
xlabel('ITPC'), ylabel('wITPC')
axis square

subplot(326)
plot(itpc,witpc_z,'.')
xlabel('ITPC'), ylabel('wITPCz')
axis square

%% end.
