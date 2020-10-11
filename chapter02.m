%% Analyzing Neural Time Series Data
% Matlab code for Chapter 2
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code.

%% Figure 2.1

load sampleEEGdata.mat

nTrials = 6; % modifiable

figure
data = zeros(nTrials,EEG.pnts);

% wavelet... more on this in Chapters 13 and 14
wavetime   = -1:1/EEG.srate:1;
n_conv     = length(wavetime)+EEG.pnts-1;
waveletfft = fft(exp(2*1i*pi*10.*wavetime) .* exp(-wavetime.^2./(2*(5/(2*pi*10))^2))/10,n_conv);
data10hz   = zeros(nTrials,EEG.pnts);

for triali=1:nTrials
    
    % create single trial "ERP"
    data(triali,:) = .5*sin(2*pi*6.*EEG.times/1000 + 2*pi*triali/nTrials-pi) + randn(1,EEG.pnts)/6;
    % add non-phase-locked stimulus potential (note distributed phases in sine wave)
    data(triali,260:360) = data(triali,260:360) + sin(2*pi*10.*EEG.times(260:360)/1000 + 2*pi*triali/nTrials-pi) + randn(1,101)/5;
    
    
    % plot data from this trial
    subplot(nTrials,3,(triali-1)*3+1)
    plot(EEG.times,data(triali,:))
    set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])
    
    % plot ERP from trial 1 to current
    subplot(nTrials,3,(triali-1)*3+2)
    plot(EEG.times,mean(data(1:triali,:),1))
    set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])
    
    % convolve with 10 Hz wavelet (more on convolution in a few chapters...)
    convolution_result_fft = ifft(waveletfft.*fft(data(triali,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft = convolution_result_fft(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz(triali,:) = abs(convolution_result_fft).^2;
    
    % plot 10 Hz power
    subplot(nTrials,3,(triali-1)*3+3)
    plot(EEG.times,mean(data10hz(1:triali,:),1))
    set(gca,'xlim',[-250 850],'ylim',[-.1 .8])
end

%% Figure 2.2
% (This code involves performing convolution with a complex Morlet wavelet,
% which you will learn about in Chapters 10-13.)

srate=1000;

time=(0:1/srate:10);

DCoffset=-.5;

% create multi-frequency signal
a = sin(2*pi*10.*time);             % part 1 of signal (high frequency)
b = .1*sin(2*pi*.3*time)+DCoffset;  % part 2 of signal (low frequency)

data = a.*b; % combined signal
data = data + (2*sin(2*pi*3*time) .* sin(2*pi*.07*time)*.1+DCoffset);

% morlet wavelet convolution (more on this in later chapters)
num_frex = 40;
min_freq =  2;
max_freq = 20;

Ldata  = length(data);
Ltapr  = length(data);
Lconv1 = Ldata+Ltapr-1;
Lconv  = pow2(nextpow2(Lconv1));

frex=logspace(log10(min_freq),log10(max_freq),num_frex);

% initialize
tf=zeros(num_frex,length(data));
datspctra = fft(data,Lconv);

s=4./(2*pi.*frex);
t=-((length(data)-1)/2)/srate:1/srate:((length(data)-2)/2)/srate+1/srate;

for fi=1:length(frex)
    
    wavelet=exp(2*1i*pi*frex(fi).*t).*exp(-t.^2./(2*s(fi)^2));
    
    m = ifft(datspctra.*fft(wavelet,Lconv),Lconv);
    m = m(1:Lconv1);
    m = m(floor((Ltapr-1)/2):end-1-ceil((Ltapr-1)/2));
    
    tf(fi,:) = abs(m).^2;
end


figure

subplot(221)
plot(a)
set(gca,'xlim',[1 8]*1000,'ylim',[-1 1],'xtick',0:1000:10000,'xticklabel',0:10);
title('10 Hz signal, DC=0')

subplot(222)
plot(b)
set(gca,'xlim',[1 8]*1000,'ylim',[-1 1],'xtick',0:1000:10000,'xticklabel',0:10);
title([ '.3 Hz signal, DC=' num2str(DCoffset) ])

subplot(223)
plot(data)
set(gca,'xlim',[1 8]*1000,'ylim',[-1 1],'xtick',0:1000:10000,'xticklabel',0:10);
title('Time-domain signal')

subplot(224)
imagesc(1:length(data),[],tf);
set(gca,'xlim',[1 8]*1000,'ydir','normal','ytick',1:8:num_frex,'yticklabel',round(frex(1:8:end)),'xtick',0:1000:10000,'xticklabel',0:10)
title('Time-frequency representation')


%% Figure 2.3

chan2plot = 'pz'; % you can pick any electrode (type {EEG.chanlocs.labels} for all electrodes)

% compute ERP (time-domain trial average from selected electrode)
erp = squeeze(mean(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),3));


% low-pass filter data (uses signal processing toolbox; you'll learn about
% how filtering works in chapter 14. If you don't have the signal processing 
% toolbox, just comment out the next eight lines; they are not necessary.)
nyquist       = EEG.srate/2;
filter_cutoff = 40;  % Hz
trans_width   = 0.1; % transition width, in fraction of 1

ffrequencies  = [ 0 filter_cutoff filter_cutoff*(1+trans_width) nyquist ]/nyquist;
idealresponse = [ 1 1 0 0 ];
filterweights = firls(100,ffrequencies,idealresponse);
filtered_erp  = filtfilt(filterweights,1,double(erp));


figure
% plot ERP
plot(EEG.times,filtered_erp,'k.-')

% now down-sample and plot
times2plot = dsearchn(EEG.times',(-200:40:1000)');
hold on
plot(EEG.times(1:5:end),filtered_erp(1:5:end),'mo-')

set(gca,'xlim',[-200 1000],'ydir','r')
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
title([ 'ERP from electrode ' chan2plot ])
legend({'256 Hz';'50 Hz'})

%%
