%% Analyzing Neural Time Series Data
% Matlab code for Chapter 12
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 12.1

load sampleEEGdata

time = -1:1/EEG.srate:1;
f = 4; % frequency of sine wave in Hz

% create sine wave (actually, a cosine wave, for reasons that will become
% clear in Chapter 13)
sine_wave = cos(2*pi*f.*time);

% make a Gaussian
s=4/(2*pi*f);
gaussian_win = exp(-time.^2./(2*s^2));

figure
plot(time,sine_wave.*gaussian_win)

%% Figure 12.2

figure
subplot(511)
plot(squeeze(EEG.data(47,:,1)))
axis tight

subplot(512)
sine_wave = cos(2*pi*12.*time);
plot(sine_wave)
axis tight

subplot(513)
% boxcar envelope
boxcar = zeros(size(sine_wave));
midpoint = (length(time)-1)/2;
boxcar(midpoint-round(EEG.srate/12/5):midpoint+round(EEG.srate/12/1.25)) = 1;
plot(sine_wave.*boxcar)
axis tight

subplot(514)
% boxcar envelope
boxcar = zeros(size(sine_wave));
midpoint = (length(time)-1)/2;
boxcar(midpoint-50:midpoint+50) = 1;
plot(sine_wave.*boxcar)
axis tight

subplot(515)
% redefine gaussian for new sine wave
s=1.5/(2*pi*12); gaussian_win = exp(-time.^2./(2*s^2));
plot(time,sine_wave.*gaussian_win)
axis tight

%% Figure 12.3

srate = 500; % sampling rate in Hz
f = 10; % frequency of the sine wave in Hz
time = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate

sine_wave = exp(2*pi*1i*f.*time); % complex wavelet

% make a Gaussian
s=6/(2*pi*f);
gaussian_win = exp(-time.^2./(2*s^2));

% and together they make a wavelet!
wavelet = sine_wave .* gaussian_win;

% make plots containing each component
figure
subplot(311)
plot(time,real(sine_wave)) % only plot the real component... we�ll learn about this later
title('Sine wave')

subplot(312)
plot(time,gaussian_win) % plots the Gaussian window
title('Gaussian window')
 
subplot(313)
plot(time,real(wavelet)); % plots the wavelet
title('My first wavelet!')
xlabel('Time (ms)')

%% Figure 12.4

num_wavelets      =  80;  % number of frequency bands
lowest_frequency  =   2;  % in Hz
highest_frequency = 100;  % in Hz

frequencies=linspace(lowest_frequency,highest_frequency,num_wavelets);
% note: the "linspace" function creates linearly spaced numbers between the first and second 
% inputs, with the number of steps corresponding to the third input. 
figure, plot(frequencies,'-*')
xlabel('Frequency order')
ylabel('Frequency in Hz')


% initialize wavelet family
wavelet_family = zeros(num_wavelets,length(time));
 
% Loop through frequencies and make a family of wavelets.
for fi=1:num_wavelets
    
    % create a sine wave at this frequency
    sinewave = exp(2*1i*pi*frequencies(fi).*time); % the "1i" makes it a complex wavelet
    
    % create a Gaussian window
    gaus_win = exp(-time.^2./(2*(6/(2*pi*frequencies(fi)))^2));
    
    % create wavelet via element-by-element multiplication of the sinewave and gaussian window
    wavelet_family(fi,:) = sinewave.*gaus_win;
    
    % note that you can also do this on one line:
    wavelet_family(fi,:) = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(6/(2*pi*frequencies(fi)))^2));
end

% Plot a few wavelets
figure
subplot(2,1,1)
plot(time,real(wavelet_family(1:round(rand*30):end,:))') 
title('A few wavelets...')
 
% Note that in the subplot command you don't need commas if you have fewer than 10 
% rows/cols and if you are not using any variables.
subplot(212)
plot(time,real(wavelet_family(30,:)))
hold on
plot(time,imag(wavelet_family(30,:)),':')
title('Real and imaginary parts of one wavelet')
legend({'real';'imaginary'})


% finally, image the wavelet family.
figure
imagesc(time,frequencies,real(wavelet_family))
axis xy % equivalent to "set(gca,'ydir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% bonus feature

% Running this cell will generate a movie that shows how a line -- in this
% case, a wavelet -- goes from a line that goes up and down, to a colored
% 'flat' line in a 2-D image. This is the principle that underlies viewing
% time-frequency plots. 

figure, set(gcf,'color','k')
surf(repmat(real(wavelet),2,1))
shading interp
axis off
axis([0 length(time) -20 21 -1 1])
view([ -4 4 ])
 
for i=4:2:90
    view([ -4 i ])
    pause(.1)
end
 
rotate3d

%% Figure 12.5

% EEG data from one trial (electrode FCz)
eegdata = squeeze(EEG.data(47,:,10));

% create wavelet
time = -1:1/EEG.srate:1;
f = 6; % frequency of sine wave in Hz
sine_wave = exp(1i*2*pi*f.*time);
s = 4.5/(2*pi*f); 
gaussian_win = exp(-time.^2./(2*s^2));
wavelet = sine_wave .* gaussian_win;
% half of the wavelet size, useful for chopping off edges after convolution.
halfwaveletsize = ceil(length(wavelet)/2);

% convolve with data
% compute Gaussian
n_conv = length(wavelet) + EEG.pnts - 1;

fft_w = fft(wavelet,n_conv);
fft_e = fft(eegdata,n_conv);
ift   = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
wavelet_conv_data = real(ift(halfwaveletsize:end-halfwaveletsize+1));

% create filter and apply to data 
% (more on how to interpret this code in a few chapters!)
nyquist       = EEG.srate/2;
transition_width = 0.2; % percent
filter_low    = 4; % Hz
filter_high   = 8; % Hz
ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(round(3*(EEG.srate/filter_low)),ffrequencies,idealresponse);
eeg_4to8      = filtfilt(filterweights,1,double(eegdata));



% now plot all the pieces
figure

plot(EEG.times,eegdata)
hold on
plot(EEG.times,wavelet_conv_data,'r','linew',2)
plot(EEG.times,eeg_4to8,'m','linew',2)
set(gca,'xlim',[-200 1200],'ydir','r')
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
legend({'Raw data';'wavelet convolved';'band-pass filtered'},0)

%% Figure 12.6

% make a theta-band-centered wavelet
time   = -1:1/EEG.srate:1;
n_conv = EEG.pnts + length(time) - 1;
n2p1   = floor(n_conv/2)+1; % n2p1 = n/2+1

f = 6;
s = 6/(2*pi*f);
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2));
halfwaveletsize = ceil(length(wavelet)/2);


eegdata = squeeze(EEG.data(47,:,10));

figure

subplot(311)
plot(EEG.times,eegdata)
set(gca,'xlim',[-500 1200])

subplot(323)
fft_w = fft(wavelet,n_conv);
hz    = linspace(0,EEG.srate/2,n2p1);
plot(hz,abs(fft_w(1:n2p1))./max(abs(fft_w(1:n2p1))),'k')
hold on

fft_e = fft(eegdata,n_conv);
hz    = linspace(0,EEG.srate/2,n2p1);
plot(hz,abs(fft_e(1:n2p1))./max(abs(fft_e(1:n2p1))),'r')
set(gca,'xlim',[0 40],'ylim',[0 1.05])
title('individual power spectra')

subplot(324)
plot(hz,abs(fft_e(1:n2p1)).*abs(fft_w(1:n2p1)))
set(gca,'xlim',[0 40])

subplot(313)
plot(EEG.times,eegdata)
hold on
ift = ifft(fft_e.*fft_w,n_conv)*sqrt(s)/10;
plot(EEG.times,real(ift(halfwaveletsize:end-halfwaveletsize+1)),'r')
set(gca,'xlim',[-500 1200])

%% Figure 12.7

% create 10 Hz wavelet (kernel)
time = -EEG.pnts/EEG.srate/2 : 1/EEG.srate : EEG.pnts/EEG.srate/2-1/EEG.srate;
f    = 10; % frequency of sine wave in Hz
s    = 4/(2*pi*f);
wavelet = cos(2*pi*f.*time) .* exp(-time.^2./(2*s^2));

% signal is one sine cycle
timeS  = 0:1/EEG.srate:(1/f); % one cycle is 1/f
signal = sin(2*pi*f.*timeS);

% now zero-pad signal
signal = [ zeros(1,EEG.pnts/2-length(timeS)/2) signal zeros(1,EEG.pnts/2-length(timeS)/2) ];

figure

% plot waves
subplot(321)
plot(wavelet)
set(gca,'xlim',[200 length(time)-200])

subplot(323)
plot(signal)
set(gca,'xlim',[200 length(time)-200])

subplot(325)
plot(conv(wavelet,signal,'same'))
set(gca,'xlim',[200 length(time)-200],'ylim',[-12 12])


% now plot dot products at selected phase lags
subplot(322)
plot(wavelet(round(100/f)-2:end),'r')
hold on
plot(signal)
set(gca,'xlim',[200 length(time)-200])
title([ 'dot product: ' num2str( fix(sum(wavelet(round(100/f)-2:end).*signal(1:end-round(100/f)+3))) ) ])

subplot(324)
plot(wavelet(round(2.3*(100/f)-2):end),'r')
hold on
plot(signal)
set(gca,'xlim',[200 length(time)-200])
title([ 'dot product: ' num2str( fix(sum(wavelet(round(2.3*(100/f)-2):end).*signal(1:end-round(2.3*(100/f)-3))) )) ])

subplot(326)
plot(wavelet,'r')
hold on
plot(signal)
set(gca,'xlim',[200 length(time)-200])
title([ 'dot product: ' num2str( fix(sum(wavelet.*signal)) ) ])

%% end
