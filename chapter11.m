%% Analyzing Neural Time Series Data
% Matlab code for Chapter 11
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 11.1

srate = 1000; % sampling rate of 1 kHz
time  = -1:1/srate:1; 
freq  = 10; % in Hz
amp   = 2; % amplitude, or height of the sine wave

sine_wave = amp.*sin(2*pi*freq.*time); % Note that you need the .* for point-wise vector multiplication.

figure
plot(time,sine_wave)
set(gca,'ylim',[-3 3]) % this adjusts the y-axis limits for visibility
title('My first sine wave!')

%% Figure 11.2

% define a sampling rate
srate = 500;
 
% list some frequencies
frex = [ 3   10   5   15   35 ];

% list some random amplitudes... make sure there are 
% the same number of amplitudes as there are frequencies!
amplit = [ 5   15   10   5   7 ];

% phases... list some random numbers between -pi and pi
phases = [  pi/7  pi/8  pi  pi/2  -pi/4 ];

% define time...
time=-1:1/srate:1;


% now we loop through frequencies and create sine waves
sine_waves = zeros(length(frex),length(time)); % remember: always initialize!
for fi=1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*frex(fi).*time + phases(fi));
end

% now plot each wave separately
figure
for fi=1:length(frex)
    subplot(length(frex),1,fi)
    plot(sine_waves(fi,:),'linew',2)
    axis([ 0 length(time) -max(amplit) max(amplit) ])
end

% now plot the result
figure
plot(sum(sine_waves))
axis tight
title('sum of sine waves')

%% Figure 11.3

figure
set(gcf,'Name','Sum of sine waves plus random noise.')
plot(sum(sine_waves+5*randn(size(sine_waves))))
axis([ 0 1020 -40 50 ]) % this sets the x-axis (first two numbers) and y-axis (last two numbers) limits
title('sum of sine waves plus white noise')

%% Figure 11.4

time=-1:1/srate:1;

% create three sine waves
s1 = sin(2*pi*3*time);
s2 = 0.5*sin(2*pi*8*time);
s3 = s1+s2;

% plot the sine waves
figure
for i=1:3
    subplot(2,3,i)
    
    % plot sine waves, using the eval command (evaluate the string)
    eval([ 'plot(time,s' num2str(i) ')' ]);
    set(gca,'ylim',[-1.6 1.6],'ytick',-1.5:.5:1.5)
    
    % plot power
    subplot(2,3,i+3)
    f  = eval([ 'fft(s' num2str(i) ')/length(time)' ]);
    hz = linspace(0,srate/2,floor(length(time)/2)+1);
    bar(hz,abs(f(1:length(hz))*2))
    set(gca,'xlim',[0 11],'xtick',0:10,'ylim',[0 1.2])
end

%% Figure 11.5

N       = 10;         % length of sequence
data    = randn(1,N); % random numbers
srate   = 200;        % sampling rate in Hz
nyquist = srate/2;    % Nyquist frequency -- the highest frequency you can measure in the data


% initialize Fourier output matrix
fourier = zeros(size(data)); 

% These are the actual frequencies in Hz that will be returned by the
% Fourier transform. The number of unique frequencies we can measure is
% exactly 1/2 of the number of data points in the time series (plus DC). 
frequencies = linspace(0,nyquist,N/2+1);
time = ((1:N)-1)/N;

% Fourier transform is dot-product between sine wave and data at each frequency
for fi=1:N
    sine_wave   = exp(-1i*2*pi*(fi-1).*time);
    fourier(fi) = sum(sine_wave.*data);
end
fourier=fourier/N;


figure
subplot(221)
plot(data,'-o')
set(gca,'xlim',[0 N+1])
title('Time domain representation of the data')

subplot(222)
plot3(frequencies,angle(fourier(1:N/2+1)),abs(fourier(1:N/2+1)).^2,'-o','linew',3)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase')
zlabel('power')
title('3-D representation of the Fourier transform')
view([20 20])

subplot(223)
bar(frequencies,abs(fourier(1:N/2+1)).^2)
set(gca,'xlim',[-5 105])
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum derived from discrete Fourier transform')

subplot(224)
bar(frequencies,angle(fourier(1:N/2+1)))
set(gca,'xlim',[-5 105])
xlabel('Frequency (Hz)')
ylabel('Phase angle')
set(gca,'ytick',-pi:pi/2:pi)
title('Phase spectrum derived from discrete Fourier transform')

%% Figure 11.6

% Compute sine waves and sum to recover the original time series
reconstructed_data = zeros(size(data));
for fi=1:N
    % scale sine wave by fourier coefficient
    sine_wave = fourier(fi)*exp(1i*2*pi*(fi-1).*time);
    % sum sine waves together (take only real part)
    reconstructed_data = reconstructed_data + real(sine_wave);
end

figure
plot(data,'-o')
hold on
plot(reconstructed_data,'r-*')

legend({'original data';'inverse Fourier transform data'})

%% Figure 11.7

fft_data = fft(data)/N;

figure
subplot(131)
plot(frequencies,abs(fourier(1:N/2+1)).^2,'*-')
hold on
plot(frequencies,abs(fft_data(1:N/2+1)).^2,'ro-','markersize',8)
% make plot look nice
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum derived from discrete Fourier transform and from FFT')
axis square
legend({'time-domain Fourier';'FFT'})

subplot(132)
plot(frequencies,angle(fourier(1:N/2+1)),'*-')
hold on
plot(frequencies,angle(fft_data(1:N/2+1)),'ro-','markersize',8)
% make plot look nice
xlabel('Frequency (Hz)')
ylabel('Phase')
set(gca,'ytick',-pi:pi/2:pi)
title('Phase spectrum derived from discrete Fourier transform and from FFT')
axis square

subplot(133)
plot(reconstructed_data,'*-')
hold on
plot(ifft(fft(data)),'ro-','markersize',8)
% make plot look nice
xlabel('Time')
ylabel('Amplitude')
title('Manual inverse Fourier transform and ifft')
axis square

%% Figure 11.9

% list some frequencies
frex = [ 3 10 5 7 ];

% list some random amplitudes
amplit = [ 5 15 10 5 ];

% phases... 
phases = [  pi/7  pi/8  pi  pi/2 ];

% create a time series of sequenced sine waves
srate = 500;
time = -1:1/srate:1;
stationary  = zeros(1,length(time)*length(frex));
nonstationary = zeros(1,length(time)*length(frex));

for fi=1:length(frex)
    
    % compute sine wave
    temp_sine_wave = amplit(fi) * sin(2*pi*frex(fi).*time + phases(fi));
    
    % enter into stationary time series
    stationary = stationary + repmat(temp_sine_wave,1,length(frex));
    
    % optional change of amplitude over time
    temp_sine_wave = temp_sine_wave.*(time+1);
    
    % determine start and stop indices for insertion of sine wave
    start_idx = (fi-1)*length(time)+1;
    stop_idx  = (fi-1)*length(time)+length(time);
    
    % enter into non-stationary time series
    nonstationary(start_idx:stop_idx) = temp_sine_wave;
end

figure

% plot stationary signal
subplot(221)
plot(stationary,'r')
set(gca,'xlim',[1 length(stationary)])
title('stationary signal')

% plot non-stationary signal
subplot(222)
plot(nonstationary)
set(gca,'xlim',[1 length(nonstationary)])
title('non-stationary signal')

% perform FFT and plot
frequencies       = linspace(0,srate/2,length(nonstationary)/2+1);
fft_nonstationary = fft(nonstationary)/length(nonstationary);
fft_stationary    = fft(stationary)/length(stationary);
subplot(212)
plot(frequencies,abs(fft_stationary(1:length(frequencies)))*2,'r')
hold on
plot(frequencies,abs(fft_nonstationary(1:length(frequencies)))*2)
set(gca,'xlim',[0 max(frex)*2])
legend({'Power stationary';'Power non-stationary'})

%% Figure 11.10

% these figures produce the unassembled components of figure 10

load sampleEEGdata
eegdat4convol = squeeze(EEG.data(47,:,1));

% create Gaussian (you'll learn more about this formula in the next chapter)
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;

figure
subplot(211)
plot(eegdat4convol)


subplot(212)
plot(gaussian)


figure
subplot(211)
plot(conv(eegdat4convol,gaussian,'same'))

subplot(212)
plot(abs(fft(conv(eegdat4convol,gaussian,'same'))))

figure
subplot(211)
plot(abs(fft(eegdat4convol)))

subplot(212)
plot(abs(fft(gaussian)))

%% Figure 11.11

srate = 1000;
time  = -.5:1/srate:.5-1/srate;
f     = 20;
fg    = [15 5];
s     = sin(2*pi*f*time);

for i=1:2
    
    % compute Gaussian
    g = exp((-time.^2)/(2*(4/(2*pi*fg(i))^2)))/fg(i);
    
    
    figure
    
    subplot(411)
    plot(time,s)
    title('Sine wave (signal)')
    set(gca,'ylim',[-1.1 1.1])
    
    subplot(412)
    plot(time,g)
    title('Gaussian (kernel)')
    
    subplot(413)
    plot(time,conv(s,g,'same'))
    set(gca,'ylim',[-1.1 1.1])
    title('result of convolution')
    
    subplot(427)
    fft_s = abs(fft(s));
    fft_s = fft_s(1:floor(length(fft_s)/2)+1)./max(fft_s(1:floor(length(fft_s)/2)+1));
    bar(0:500,fft_s,'r')
    hold on
    
    fft_g = abs(fft(g));
    fft_g = fft_g(1:floor(length(fft_g)/2)+1)./max(fft_g(1:floor(length(fft_g)/2)+1));
    plot(0:500,fft_g)
    set(gca,'xlim',[0 40],'ylim',[0 1.05])
    title('individual power spectra')
    
    subplot(428)
    bar(0:500,fft_g.*fft_s)
    set(gca,'xlim',[0 40],'ylim',[0 .035])
    title('multiplied power spectra')
end

% inset scaling: axis([15 25 -.01 .11])

%% Figure 11.12

% create Gaussian
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;

figure

% plot EEG data
subplot(411)
plot(EEG.times,eegdat4convol)

% plot Gaussian
subplot(412)
plot(time,gaussian)

% plot result of convolution
subplot(413)
plot(EEG.times,eegdat4convol,'r')
hold on
plot(EEG.times,conv(eegdat4convol,gaussian,'same'))

subplot(427)
nfft = length(eegdat4convol);
fft_s = abs(fft(eegdat4convol,nfft));
fft_s = fft_s(1:floor(nfft/2)+1);
f = linspace(0,EEG.srate/2,floor(nfft/2)+1);
plot(f,fft_s./max(fft_s),'r')
hold on

fft_g = abs(fft(gaussian,nfft));
fft_g = fft_g(1:floor(nfft/2)+1);
plot(f,fft_g./max(fft_g))

set(gca,'xlim',[0 60])

subplot(428)
plot(f,fft_s.*fft_g)
set(gca,'xlim',[0 60])

%% end.
