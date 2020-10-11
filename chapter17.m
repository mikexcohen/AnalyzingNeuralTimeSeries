%% Analyzing Neural Time Series Data
% Matlab code for Chapter 17
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% load data

load sampleEEGdata

%% figure 17.1

num_frex = 6;
frex = logspace(log10(4),log10(30),num_frex);
s    = 6./(2*pi*frex);
time  = -1:1/EEG.srate:1;

% note that you don't need the wavelet itself, you need the FFT of the wavelet
mwaves = zeros(num_frex,length(time));
swaves = zeros(num_frex,length(time));
for fi = 1:num_frex
    % Create Morlet wavelets and s-transforms. Note that scaling factors
    % are omitted to focus on the shape of the function.
    mwaves(fi,:) = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)));
    swaves(fi,:) = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2.*frex(fi)^2/2);
end


time      = -1:1/EEG.srate:1;
n_wavelet = length(time);
n_data    = EEG.pnts*EEG.trials;
n_conv    = n_wavelet+n_data-1;
half_wave = (length(time)-1)/2;

eegfft    = fft(reshape(EEG.data(64,:,:),1,EEG.pnts*EEG.trials),n_conv);

% convolution
eegconv   = ifft(fft(mwaves(5,:),n_conv).*eegfft);
eegconv   = eegconv(half_wave+1:end-half_wave);

% reshape to time X trials
eegpower  = log10(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2);

figure

plot(EEG.times,eegpower(:,50))
set(gca,'xlim',[-400 1200])
threshold = prctile(eegpower(:),95);
hold on
plot(get(gca,'xlim'),[threshold threshold],'k')
title([ 'power at electrode ' EEG.chanlocs(64).labels ' at ' num2str(frex(5)) ' Hz' ])

%% figure 17.2

figure
for i=1:num_frex
    subplot(2,3,i)
    plot(time,real(mwaves(i,:)))
    hold on
    plot(time,real(swaves(i,:)),'r.')
    title([ num2str(frex(i)) ' Hz' ])
end

%% end.
