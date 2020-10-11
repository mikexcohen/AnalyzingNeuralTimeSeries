%% Analyzing Neural Time Series Data
% Matlab code for Chapter 5
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% load sample data

load sampleEEGdata.mat

%% time-locked and non-phase-locked

nTrials = 4;

figure
data = zeros(nTrials,EEG.pnts);

% wavelet... more on this in Chapters 13 and 14
wavetime   = -1:1/EEG.srate:1;
n_conv     = length(wavetime)+EEG.pnts-1;
waveletfft = fft(exp(2*1i*pi*10.*wavetime) .* exp(-wavetime.^2./(2*(5/(2*pi*10))^2))/10,n_conv);
data10hz   = zeros(nTrials,EEG.pnts);

for triali=1:nTrials
    
    % create single trial ERP as sine wave plus noise
    data(triali,:) = .15*sin(2*pi*6.*EEG.times/1000 + 2*pi*triali/nTrials-pi) + randn(1,EEG.pnts)/6;
    % add non-phase-locked stimulus potential
    data(triali,260:360) = data(triali,260:360) + sin(2*pi*10.*EEG.times(260:360)/1000 + 2*pi*triali/nTrials-pi) + randn(1,101)/5;
    
    % convolve with 10Hz wavelet (more on convolution in a few chapters...)
    convolution_result_fft = ifft(waveletfft.*fft(data(triali,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft = convolution_result_fft(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz(triali,:) = abs(convolution_result_fft)*2;
    
    
    % plot single trials
    subplot(nTrials,3,(triali-1)*3+1)
    plot(EEG.times,data(triali,:))
    hold on
    plot(EEG.times,data10hz(triali,:),'r')
    set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])
end

subplot(nTrials,3,(triali-1)*3+2)
plot(EEG.times,mean(data,1))
hold on
plot(EEG.times,mean(data10hz,1),'r')
set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])

%% time-locked and phase-locked

nTrials = 4;

figure
data = zeros(nTrials,EEG.pnts);

% wavelet... more on this in Chapters 13 and 14
wavetime   = -1:1/EEG.srate:1;
n_conv     = length(wavetime)+EEG.pnts-1;
waveletfft = fft(exp(2*1i*pi*10.*wavetime) .* exp(-wavetime.^2./(2*(5/(2*pi*10))^2))/10,n_conv);
data10hz   = zeros(nTrials,EEG.pnts);

for triali=1:nTrials
    
    % create single trial ERP
    data(triali,:) = .15*sin(2*pi*6.*EEG.times/1000 + 2*pi*triali/nTrials-pi) + randn(1,EEG.pnts)/6;
    % add non-phase-locked stimulus potential
    data(triali,260:360) = data(triali,260:360) + sin(2*pi*10.*EEG.times(260:360)/1000) + randn(1,101)/5;
    
    % convolve with 10Hz wavelet (more on convolution in a few chapters...)
    convolution_result_fft = ifft(waveletfft.*fft(data(triali,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft = convolution_result_fft(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz(triali,:) = abs(convolution_result_fft)*2;
    
    
    % plot single trials
    subplot(nTrials,3,(triali-1)*3+1)
    plot(EEG.times,data(triali,:))
    hold on
    plot(EEG.times,data10hz(triali,:),'r')
    set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])
end

subplot(nTrials,3,(triali-1)*3+2)
plot(EEG.times,mean(data,1))
hold on
plot(EEG.times,mean(data10hz,1),'r')
set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])

%% non-time-locked and phase-locked

nTrials = 4;

figure
data = zeros(nTrials,EEG.pnts);

% wavelet... more on this in Chapters 13 and 14
wavetime   = -1:1/EEG.srate:1;
n_conv     = length(wavetime)+EEG.pnts-1;
waveletfft = fft(exp(2*1i*pi*10.*wavetime) .* exp(-wavetime.^2./(2*(5/(2*pi*10))^2))/10,n_conv);
data10hz   = zeros(nTrials,EEG.pnts);

for triali=1:nTrials
    
    % create single trial ERP
    data(triali,:) = .15*sin(2*pi*6.*EEG.times/1000 + 2*pi*triali/nTrials-pi) + randn(1,EEG.pnts)/6;
    % add non-phase-locked stimulus potential
    eventtime = randperm(80)+240;
    eventtime = eventtime(1):eventtime(1)+80;
    data(triali,eventtime) = data(triali,eventtime) + sin(2*pi*10.*EEG.times(eventtime)/1000 + 2*pi*triali/nTrials-pi) + randn(1,length(eventtime))/5;
    
    % convolve with 10Hz wavelet (more on convolution in a few chapters...)
    convolution_result_fft = ifft(waveletfft.*fft(data(triali,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft = convolution_result_fft(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz(triali,:) = abs(convolution_result_fft)*2;
    
    
    % plot single trials
    subplot(nTrials,3,(triali-1)*3+1)
    plot(EEG.times,data(triali,:))
    hold on
    plot(EEG.times,data10hz(triali,:),'r')
    set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])
end

subplot(nTrials,3,(triali-1)*3+2)
plot(EEG.times,mean(data,1))
hold on
plot(EEG.times,mean(data10hz,1),'r')
set(gca,'xlim',[-250 850],'ylim',[-2.2 2.2])

%% end.
