%% Analyzing Neural Time Series Data
% Matlab code for Chapter 20
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 20.1

load sampleEEGdata

channel2plot = 'o1';

% wavelet parameters
min_freq = 2;
max_freq = 30;
num_frex = 20;

% baseline time window
baseline_time = [ -400 -100 ];

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution(1:2) = n_wavelet+n_data-1;
n_convolution(3)   = n_wavelet+EEG.pnts-1; % ERP is only one trial-length

% find sensor index
sensoridx = find(strcmpi(channel2plot,{EEG.chanlocs.labels}));

% compute ERP
erp = squeeze(mean(EEG.data(sensoridx,:,:),3));

% compute induced power by subtracting ERP from each trial
induced_EEG = squeeze(EEG.data(sensoridx,:,:)) - repmat(erp',1,EEG.trials);

% FFT of data
fft_EEG{1} = fft(reshape(EEG.data(sensoridx,:,:),1,EEG.pnts*EEG.trials),n_convolution(1)); % total 
fft_EEG{2} = fft(reshape(induced_EEG,1,EEG.pnts*EEG.trials),n_convolution(2)); % induced
fft_EEG{3} = fft(erp,n_convolution(3)); % evoked; note that it doesn't matter that the FFT is longer than the time series

% convert baseline from ms to indices
[junk,baseidx(1)] = min(abs(EEG.times-baseline_time(1)));
[junk,baseidx(2)] = min(abs(EEG.times-baseline_time(2)));

% initialize output time-frequency data
tf = zeros(4,length(frequencies),EEG.pnts);

for fi=1:length(frequencies)
    
    % create wavelet
    wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(4/(2*pi*frequencies(fi)))^2))/frequencies(fi);
    
    % run convolution for each of total, induced, and evoked
    for i=1:3
        
        % take FFT of data
        fft_wavelet = fft(wavelet,n_convolution(i));
        
        % convolution...
        convolution_result_fft = ifft(fft_wavelet.*fft_EEG{i},n_convolution(i));
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        % reshaping and trial averaging is done only on all trials
        if i<3
            convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
            
            % compute power
            tf(i,fi,:) = mean(abs(convolution_result_fft).^2,2);
        else
            % with only one trial-length, just compute power with no averaging
            tf(i,fi,:) = abs(convolution_result_fft).^2;
        end
        
        % db correct power
        tf(i,fi,:) = 10*log10( squeeze(tf(i,fi,:)) ./ mean(tf(i,fi,baseidx(1):baseidx(2)),3) );
        
        % inter-trial phase consistency on total EEG
        if i==1
            tf(4,fi,:) = abs(mean(exp(1i*angle(convolution_result_fft)),2));
        end
    end % end loop around total, evoked, induced
end % end frequency loop



analysis_labels = {'Total';'Non-phase-locked';'ERP power';'ITPC'};

% color limits
clims = [ -3 3; -3 3; -12 12; 0 .6 ];

% scale ERP for plotting
erpt = (erp-min(erp))./max(erp-min(erp));
erpt = erpt*(frequencies(end)-frequencies(1))+frequencies(1);

figure
for i=1:4
    
    subplot(2,3,i)
    contourf(EEG.times,frequencies,squeeze(tf(i,:,:)),40,'linecolor','none')
    
    set(gca,'clim',clims(i,:),'xlim',[-400 1000],'xtick',-200:200:800)
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})

    % plot ERP on top
    hold on
    plot(EEG.times,erpt,'k')
end

subplot(235)
contourf(EEG.times,frequencies,squeeze(tf(1,:,:)-tf(2,:,:)),40,'linecolor','none')
set(gca,'clim',clims(1,:),'xlim',[-400 1000],'xtick',-200:200:800)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Phase-locked')

% plot ERP on top
hold on
plot(EEG.times,erpt,'k')

%%

figure
subplot(221)
plot(reshape(tf(3,:,:),1,[]),reshape(tf(4,:,:),1,[]),'.')
xlabel('ERP power'), ylabel('ITPC')

subplot(222)
plot(reshape(tf(1,:,:)-tf(2,:,:),1,[]),reshape(tf(4,:,:),1,[]),'.')
xlabel('Phase-locked power'), ylabel('ITPC')
set(gca,'xlim',[-.3 4])

subplot(223)
plot(reshape(tf(1,:,:),1,[]),reshape(tf(2,:,:),1,[]),'.')
xlabel('Total power'), ylabel('Non-phase-locked power')
set(gca,'xlim',[-2 6],'ylim',[-2 6])

subplot(224)
plot(reshape(tf(1,:,:)-tf(2,:,:),1,[]),reshape(tf(3,:,:),1,[]),'.')
xlabel('Phase-locked power'), ylabel('ERP power')
set(gca,'xlim',[-.3 4])


%% end.
