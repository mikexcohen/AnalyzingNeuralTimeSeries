%% Analyzing Neural Time Series Data
% Matlab code for Chapter 30
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 30.1

load sampleEEGdata

channel2plot = 'o1';

% wavelet parameters
freq2plot = 25;

% other wavelet parameters
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% find sensor index
sensoridx = find(strcmpi(channel2plot,{EEG.chanlocs.labels}));

% FFT of data
fft_EEG = fft(reshape(EEG.data(sensoridx,:,:),1,EEG.pnts*EEG.trials),n_convolution);

% create wavelet and get its FFT
wavelet = exp(2*1i*pi*freq2plot.*time) .* exp(-time.^2./(2*(4/(2*pi*freq2plot))^2));
fft_wavelet = fft(wavelet,n_convolution);

convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result = reshape(convolution_result,EEG.pnts,EEG.trials);

figure
subplot(211)
plot(EEG.times,real(convolution_result(:,1))) % filtered signal from the first trial
xlabel('Time (ms)'), ylabel('Filtered signal amplitude')
set(gca,'xlim',[-200 1000])

subplot(212)
plot(EEG.times,abs(convolution_result(:,1)).^2) % power from the first trial
xlabel('Time (ms)'), ylabel('Power')
set(gca,'xlim',[-200 1000])

%% figure 30.2

load accumbens_eeg.mat
srate = 1000;

% wavelet parameters
freq2plot = 70;

% other wavelet parameters
time = -1:1/srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = length(eeg);
n_convolution = n_wavelet+n_data-1;

% FFT of data
fft_EEG = fft(eeg,n_convolution);

% create wavelet and get its FFT
wavelet = exp(2*1i*pi*freq2plot.*time) .* exp(-time.^2./(2*(4/(2*pi*freq2plot))^2));
fft_wavelet = fft(wavelet,n_convolution);

convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);

eegtime = (0:length(eeg)-1)/srate;

figure
subplot(311)
plot(eegtime,eeg)
xlabel('Time (ms)'), ylabel('broadband signal amplitude')
set(gca,'xlim',[0 1])

subplot(312)
plot(eegtime,real(convolution_result)) % filtered signal from the first trial
xlabel('Time (ms)'), ylabel('Filtered signal amplitude')
set(gca,'xlim',[0 1])

subplot(313)
plot(eegtime,abs(convolution_result).^2) % power from the first trial
xlabel('Time (ms)'), ylabel('Power')
set(gca,'xlim',[0 1])

%% Figure 30.3


% we will first test for cross-frequency coupling between two specific frequency bands
freq4phase = 10; % in Hz
freq4power = 70; 

% wavelet and FFT parameters
srate = 1000;
time = -1:1/srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = length(eeg);
n_convolution = n_wavelet+n_data-1;
fft_data = fft(eeg,n_convolution);

% wavelet for phase and its FFT
wavelet4phase = exp(2*1i*pi*freq4phase.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4phase))^2));
fft_wavelet4phase = fft(wavelet4phase,n_convolution);

% wavelet for power and its FFT
wavelet4power = exp(2*1i*pi*freq4power.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4power))^2));
fft_wavelet4power = fft(wavelet4power,n_convolution);

% get phase values
convolution_result_fft = ifft(fft_wavelet4phase.*fft_data,n_convolution);
phase = angle(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));

% get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
convolution_result_fft = ifft(fft_wavelet4power.*fft_data,n_convolution);
pwr = abs(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size)).^2;

% plot power and phase
figure
subplot(131)
plot(eegtime,phase)
hold on
plot(eegtime,(pwr-mean(pwr))/std(pwr),'r')
legend({'10 Hz phase';'70 Hz power'})
set(gca,'xlim',[0 1])
axis square

% plot power as a function of phase in polar space
subplot(132)
polar(phase,pwr,'.')

% plot histogram of power over phase
n_hist_bins = 30;

phase_edges=linspace(min(phase),max(phase),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);

for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(pwr(phase>phase_edges(i) & phase<phase_edges(i+1)));
end

subplot(133)
bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)
axis square

%% Figure 30.4

phase_bias = phase;
power_bias = pwr;

phase_bias(phase<-pi/2) = [];
power_bias(phase<-pi/2) = [];


figure
% plot power as a function of phase in polar space
subplot(131)
polar(phase,pwr,'.')
title([ 'PAC = ' num2str(round(abs(mean(pwr.*exp(1i*phase))))) ])

subplot(132)
polar(phase,pwr*10,'.')
title([ 'PAC = ' num2str(round(abs(mean(pwr*10.*exp(1i*phase))))) ])

subplot(133)
polar(phase_bias,power_bias,'.')
title([ 'PAC = ' num2str(round(abs(mean(power_bias.*exp(1i*phase_bias))))) ])

%% Figure 30.5

% observed cross-frequency-coupling (note the similarity to Euler's formula)
obsPAC = abs(mean(pwr.*exp(1i*phase)));
obsPAC_bias = abs(mean(power_bias.*exp(1i*phase_bias)));

num_iter = 1000;

permutedPAC = zeros(2,num_iter);

for i=1:num_iter
    
    % select random time point
    random_timepoint = randsample(round(length(eeg)*.8),1)+round(length(eeg)*.1);
    random_timepoint_bias = randsample(round(length(power_bias)*.8),1)+round(length(power_bias)*.1);
    
    % shuffle power
    timeshiftedpwr      = [ pwr(random_timepoint:end) pwr(1:random_timepoint-1) ];
    timeshiftedpwr_bias = [ power_bias(random_timepoint_bias:end) power_bias(1:random_timepoint_bias-1) ];
    
    % compute PAC
    permutedPAC(1,i) = abs(mean(timeshiftedpwr.*exp(1i*phase)));
    permutedPAC(2,i) = abs(mean(timeshiftedpwr_bias.*exp(1i*phase_bias)));
end

% compute PACz
pacz(1) = (obsPAC-mean(permutedPAC(1,:)))/std(permutedPAC(1,:));
pacz(2) = (obsPAC_bias-mean(permutedPAC(2,:)))/std(permutedPAC(2,:));

figure
subplot(221)
hist(permutedPAC(1,:),50);
hold on
plot([obsPAC obsPAC],get(gca,'ylim')/2,'m','linew',3)
legend({'Permuted values';'Observed value'})
xlabel('Modulation strength'), ylabel('Number of observations')
title([ 'PAC_z = ' num2str(pacz(1)) ])

subplot(222)
hist(permutedPAC(2,:),50)
hold on
plot([obsPAC_bias obsPAC_bias],get(gca,'ylim')/2,'m','linew',3)
legend({'Permuted values';'Observed value'})
xlabel('Modulation strength'), ylabel('Number of observations')
title([ 'PAC_z = ' num2str(pacz(2)) ])

% plot histogram of power over phase
n_hist_bins = 30;
phase_edges=linspace(min(phase),max(phase),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);
for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(pwr(phase>phase_edges(i) & phase<phase_edges(i+1)));
end

subplot(223)
h=bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(h,'linestyle','none'); % turn off black lines around histogram bars
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)


% plot histogram of power over phase
n_hist_bins = 30;
phase_edges=linspace(min(phase_bias),max(phase_bias),n_hist_bins+1);
amp_by_phases=zeros(1,n_hist_bins);
for i=1:n_hist_bins-1
    amp_by_phases(i) = mean(power_bias(phase_bias>phase_edges(i) & phase_bias<phase_edges(i+1)));
end

subplot(224)
h=bar(phase_edges(1:end-1),amp_by_phases,'histc');
set(h,'linestyle','none'); % turn off black lines around histogram bars
set(gca,'xlim',[phase_edges(1) phase_edges(end)])
xlabel([ 'Phase at ' num2str(freq4phase) ' Hz (rad.)' ])
ylabel([ 'Power at ' num2str(freq4power) ' Hz' ])
set(gca,'xlim',[-3.5 3.5],'xtick',-pi:pi/2:pi)

%% Figure 30.6

permutedPAC = zeros(3,num_iter);

for i=1:num_iter
    
    % Permutation method 1: select random time point
    random_timepoint = randsample(round(length(eeg)*.8),1)+round(length(eeg)*.1);
    timeshiftedpwr   = [ pwr(random_timepoint:end) pwr(1:random_timepoint-1) ];
    permutedPAC(1,i) = abs(mean(timeshiftedpwr.*exp(1i*phase)));
    
    % Permutation method 2: totally randomize power time series
    permutedPAC(2,i) = abs(mean(pwr(randperm(length(pwr))).*exp(1i*phase)));
    
    % Permutation method 3: FFT-based power time series randomization
    f = fft(pwr); % compute FFT
    A = abs(f);   % extract amplitudes
    zphs=cos(angle(f))+1i*sin(angle(f)); % extract phases
    powernew=real(ifft(A.*zphs(randperm(length(zphs))))); % recombine using randomized phases (note: use original phases to prove that this method reconstructs the original signal)
    powernew=powernew-min(powernew);
    
    permutedPAC(3,i) = abs(mean(powernew.*exp(1i*phase)));
end

% compute PACz and plot
figure
for i=1:3
    subplot(2,3,i)
    
    % plot example power time series
    switch i
        case 1
            plot(eegtime,timeshiftedpwr)
            title('H_0: Time-shifted')
        case 2
            plot(eegtime,pwr(randperm(length(pwr))))
            title('H_0: randomized')
        case 3
            plot(eegtime,powernew)
            title('H_0: FFT-derived randomization')
    end
    set(gca,'xlim',[0 eegtime(end)],'ylim',[min(pwr) max(pwr)])
    
    % plot null-hypothesis distribution
    subplot(2,3,i+3)
    pacz = (obsPAC-mean(permutedPAC(i,:)))/std(permutedPAC(i,:));
    [y,x]=hist(permutedPAC(i,:),50);
    h=bar(x,y,'histc');
    set(h,'linestyle','none');
    hold on
    plot([obsPAC obsPAC],get(gca,'ylim')/2,'m','linew',3)
    legend({'Permuted values';'Observed value'})
    xlabel('Modulation strength'), ylabel('Number of observations')
    title([ 'PAC_z = ' num2str(pacz) ])
end

%% Figure 30.7

times2plot = -200:100:1200;
freq4phase = 10; % Hz
freq4power = 25; 

cfc_numcycles  = 3;   % number of cycles at phase-frequency

pacz = zeros(size(times2plot));
itpc = zeros(size(times2plot));

% convert cfc times to indices
cfc_time_window     = cfc_numcycles*(1000/freq4phase);
cfc_time_window_idx = round(cfc_time_window/(1000/EEG.srate));

% other wavelet parameters
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% FFT of scalp EEG data
fft_EEG = fft(reshape(EEG.data(sensoridx,:,:),1,EEG.pnts*EEG.trials),n_convolution);

for timei=1:length(times2plot)
    
    cfc_centertime_idx  = dsearchn(EEG.times',times2plot(timei));
    
    % convolution for lower frequency phase
    wavelet            = exp(2*1i*pi*freq4phase.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4phase))^2));
    fft_wavelet        = fft(wavelet,n_convolution);
    convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
    convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    lower_freq_phase   = reshape(convolution_result,EEG.pnts,EEG.trials);
    
    % convolution for upper frequency power
    wavelet            = exp(2*1i*pi*freq4power.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4power))^2));
    fft_wavelet        = fft(wavelet,n_convolution);
    convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
    convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    upper_freq_power   = reshape(convolution_result,EEG.pnts,EEG.trials);
    
    
    
    % extract temporally localized power and phase from task data (not vectorized this time)
    power_ts = abs(upper_freq_power(cfc_centertime_idx-round(cfc_time_window_idx/2):cfc_centertime_idx+round(cfc_time_window_idx/2),:)).^2;
    phase_ts = angle(lower_freq_phase(cfc_centertime_idx-round(cfc_time_window_idx/2):cfc_centertime_idx+round(cfc_time_window_idx/2),:));
    
    % compute observed PAC
    obsPAC = abs(mean( power_ts(:).*exp(1i*phase_ts(:)) ));
    % compute lower frequency ITPC
    itpc(timei) = mean(abs(mean(exp(1i*phase_ts),2)));
    
    num_iter = 1000;
    permutedPAC = zeros(1,num_iter);
    for i=1:num_iter
        
        % in contrast to the previous code, this time-shifts the power time series only within trials. Results are similar using either method.
        random_timepoint = randsample(round(cfc_time_window_idx*.8),EEG.trials,1)+round(cfc_time_window_idx*.1);
        for triali=1:EEG.trials
            power_ts(:,triali) = power_ts([random_timepoint(triali):end 1:random_timepoint(triali)-1],triali);
        end
        
        permutedPAC(i) = abs(mean( power_ts(:).*exp(1i*phase_ts(:)) ));
    end
    
    pacz(timei) = (obsPAC-mean(permutedPAC))/std(permutedPAC);
end


figure
subplot(211)
plot(times2plot,pacz,'-o','markerface','w')
set(gca,'xlim',get(gca,'xlim').*[1.15 1.05]) % open the x-limits a bit

% this next line computes the Z-value threshold at p=0.05, correcting for multiple comparisons across time points (this is a bit conservative because of temporal autocorrelation)
% if you don't have the matlab stats toolbox, use a zval of 2.7131 (p<0.05 correcting for 15 time points/comparisons)
zval = norminv(1-(.05/length(times2plot)));

hold on
plot(get(gca,'xlim'),[zval zval],'k:')
plot(get(gca,'xlim'),[0 0],'k')
xlabel('Time (ms)'), ylabel('PAC_z')

title([ 'PAC_z at electrode ' channel2plot ' between ' num2str(freq4power) ' Hz power and ' num2str(freq4phase) ' Hz phase' ])

% Also plot ITPC for comparison
subplot(212)
plot(times2plot,itpc,'-o','markerface','w')
set(gca,'xlim',get(gca,'xlim').*[1.15 1.05]) % open the x-limits a bit
title([ 'ITPC at electrode ' channel2plot ' at ' num2str(freq4phase) ' Hz' ])

%% Figure 30.8a (takes a while to run)

phase_freqs    = 2:20; % Hz
cfc_centertime = 300; % ms post-stimulus

pacz = zeros(size(phase_freqs));
cfc_centertime_idx  = dsearchn(EEG.times',cfc_centertime);

for fi=1:length(phase_freqs)
    
    % convolution for lower frequency phase
    wavelet     = exp(2*1i*pi*phase_freqs(fi).*time) .* exp(-time.^2./(2*(4/(2*pi*phase_freqs(fi)))^2));
    fft_wavelet = fft(wavelet,n_convolution);
    convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
    convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    lower_freq_phase   = reshape(convolution_result,EEG.pnts,EEG.trials);
    
    % extract temporally localized power and phase from task data (not vectorized this time)
    power_ts = abs(upper_freq_power(cfc_centertime_idx-round(cfc_time_window_idx/2):cfc_centertime_idx+round(cfc_time_window_idx/2),:)).^2;
    phase_ts = angle(lower_freq_phase(cfc_centertime_idx-round(cfc_time_window_idx/2):cfc_centertime_idx+round(cfc_time_window_idx/2),:));
    
    % compute observed PAC
    obsPAC = abs(mean( reshape(power_ts,1,[]).*exp(1i*reshape(phase_ts,1,[])) ));
    
    num_iter = 2000;
    permutedPAC = zeros(1,num_iter);
    for i=1:num_iter
        
        % in contrast to the previous code, this time-shifts the power time series only within trials. Results are similar using either method.
        random_timepoint = randsample(round(cfc_time_window_idx*.8),EEG.trials,1)+round(cfc_time_window_idx*.1);
        for triali=1:EEG.trials
            power_ts(:,triali) = power_ts([random_timepoint(triali):end 1:random_timepoint(triali)-1],triali);
        end
        
        permutedPAC(i) = abs(mean( reshape(power_ts,1,[]).*exp(1i*reshape(phase_ts,1,[])) ));
    end
    
    pacz(fi) = (obsPAC-mean(permutedPAC))/std(permutedPAC);
end

figure
plot(phase_freqs,pacz,'-o')
set(gca,'xlim',get(gca,'xlim').*[.5 1.05]) % open the x-limits a bit
xlabel('Lower frequency for phase (Hz)'), ylabel('PAC_z')

[junk,max_phase_freq] = max(pacz);
title([ 'Best lower frequency phase coupled with ' num2str(freq4power) ' Hz is ' num2str(phase_freqs(max_phase_freq)) ' Hz' ]);

%% Figure 30.8b

power_freqs = 20:5:EEG.srate/2; % Hz

pacz = zeros(size(power_freqs));

% convolution for lower frequency phase
wavelet     = exp(2*1i*pi*freq4phase.*time) .* exp(-time.^2./(2*(4/(2*pi*freq4phase))^2));
fft_wavelet = fft(wavelet,n_convolution);
convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
lower_freq_phase   = reshape(convolution_result,EEG.pnts,EEG.trials);

for fi=1:length(power_freqs)
    
    % convolution for upper frequency power
    wavelet     = exp(2*1i*pi*power_freqs(fi).*time) .* exp(-time.^2./(2*(4/(2*pi*power_freqs(fi)))^2));
    fft_wavelet = fft(wavelet,n_convolution);
    convolution_result = ifft(fft_wavelet.*fft_EEG,n_convolution);
    convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    upper_freq_power   = reshape(convolution_result,EEG.pnts,EEG.trials);
    
    % extract temporally localized power and phase from task data (not vectorized this time)
    power_ts = abs(upper_freq_power(cfc_centertime_idx-round(cfc_time_window_idx/2):cfc_centertime_idx+round(cfc_time_window_idx/2),:)).^2;
    phase_ts = angle(lower_freq_phase(cfc_centertime_idx-round(cfc_time_window_idx/2):cfc_centertime_idx+round(cfc_time_window_idx/2),:));
    
    % compute observed PAC
    obsPAC = abs(mean( reshape(power_ts,1,[]).*exp(1i*reshape(phase_ts,1,[])) ));
    
    num_iter = 2000;
    permutedPAC = zeros(1,num_iter);
    for i=1:num_iter
        
        % in contrast to the previous code, this time-shifts the power time series only within trials. Results are similar using either method.
        random_timepoint = randsample(round(cfc_time_window_idx*.8),EEG.trials,1)+round(cfc_time_window_idx*.1);
        for triali=1:EEG.trials
            power_ts(:,triali) = power_ts([random_timepoint(triali):end 1:random_timepoint(triali)-1],triali);
        end
        
        permutedPAC(i) = abs(mean( reshape(power_ts,1,[]).*exp(1i*reshape(phase_ts,1,[])) ));
    end
    
    pacz(fi) = (obsPAC-mean(permutedPAC))/std(permutedPAC);
end

figure
plot(power_freqs,pacz,'-o')
set(gca,'xlim',[power_freqs(1)-3 power_freqs(end)+3])
xlabel('Upper frequency for power (Hz)'), ylabel('PAC_z')

[junk,max_power_freq] = max(pacz);
title([ 'Best upper frequency power coupled with ' num2str(freq4phase) ' Hz is ' num2str(power_freqs(max_power_freq)) ' Hz' ]);

%% Figure 30.9

% You have to figure this one out on your own!

%% Figure 30.10

% This figure is created by combining the code for figure 30.8. You need
% two loops, one for lower-frequency phase and one for upper-frequency
% power. Compute PAC at each phase-power pair, and then make an image of
% the resulting (z-scored) PAC values. 

%% Figure 30.11

freqs4phase =  1:20;
freqs4power = 25:EEG.srate/2;

powcycles_per_phscycles = zeros(length(freqs4power),length(freqs4phase));

for phsi=1:length(freqs4phase)
    for powi=1:length(freqs4power)
        % number of power cycles per phase cycles, scaled by sampling rate
        powcycles_per_phscycles(powi,phsi) = ( freqs4power(powi)/freqs4phase(phsi) ) / ( 1000/EEG.srate );
    end
end

figure
imagesc(freqs4phase,freqs4power,log10(powcycles_per_phscycles))
set(gca,'clim',log10([1 20]),'ydir','n')
colorbar
colormap gray

%% Figure 30.12

% wavelet parameters
upperfreq = 70;
lowerfreq = 12;

% other wavelet parameters
time          = -1:1/srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = length(eeg);
n_convolution = n_wavelet+n_data-1;

% FFT of data
fft_EEG = fft(eeg,n_convolution);

% convolution for lower frequency phase (with 4 cycles)
waveletL = exp(2*1i*pi*lowerfreq.*time) .* exp(-time.^2./(2*(4/(2*pi*lowerfreq))^2));
fft_waveletL = fft(waveletL,n_convolution);
convolution_result = ifft(fft_waveletL.*fft_EEG,n_convolution);
lowerfreq_phase = angle(convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size));

% convolution for upper frequency (with 4 cycles)
waveletH = exp(2*1i*pi*upperfreq.*time) .* exp(-time.^2./(2*(4/(2*pi*upperfreq))^2));
fft_waveletH = fft(waveletH,n_convolution);
convolution_result = ifft(fft_waveletH.*fft_EEG,n_convolution);
upperfreq_amp = abs(convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size));

% filter the upper frequency power in the lower frequency range (in this example, filter at 12 Hz)
% then take the angle from the result of convolution
% (you could achieve the same result by band-pass filtering around 12 Hz and taking the Hilbert transform)
convolution_result  = ifft(fft_waveletL.*fft(upperfreq_amp,n_convolution),n_convolution);
upperfreq_amp_phase = angle(convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size));

figure

% plot accumbens alpha phase
subplot(411)
plot(eegtime,lowerfreq_phase)
set(gca,'xlim',[0 2])
title([ num2str(lowerfreq) ' Hz phase' ])

subplot(412)
plot(eegtime,upperfreq_amp)
set(gca,'xlim',[0 2])
title([ num2str(upperfreq) ' Hz power' ])

subplot(413)
plot(eegtime,real(convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size)))
set(gca,'xlim',[0 2])
title([ num2str(upperfreq) ' Hz power filtered at ' num2str(lowerfreq) ' Hz' ])

subplot(414)
plot(eegtime,upperfreq_amp_phase)
hold on
plot(eegtime,lowerfreq_phase,'r')
legend({'upper';'lower'})
title([ 'Phase of ' num2str(lowerfreq) ' Hz component in ' num2str(upperfreq) ' Hz power' ])
set(gca,'xlim',[0 2])

% compute synchronization
phasephase_synch = abs(mean(exp(1i*( lowerfreq_phase-upperfreq_amp_phase ))));

disp([ 'Phase-phase coupling between ' num2str(lowerfreq) ' Hz and ' num2str(upperfreq) ' Hz is ' num2str(phasephase_synch) ])

%% end.
