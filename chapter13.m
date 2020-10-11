%% Analyzing Neural Time Series Data
% Matlab code for Chapter 13
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 13.1

% parameters...
srate = 500; % sampling rate in Hz
f     = 10; % frequency of wavelet in Hz
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate
s     = 6/(2*pi*f);

% and together they make a wavelet
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2)); 

figure
subplot(221)
% show the projection onto the real axis
plot3(time,real(wavelet),imag(wavelet),'m')
xlabel('Time (ms)'), ylabel('real axis')
view(0,90)
title('Projection onto real and time axes')

% show the projection onto the imaginary axis
subplot(222)
plot3(time,real(wavelet),imag(wavelet),'g')
xlabel('Time (ms)'), ylabel('imaginary axis')
view(0,0)
title('Projection onto imaginary and time axes') 
 
% plot projection onto real and imaginary axes
subplot(223)
plot3(time,real(wavelet),imag(wavelet),'k')
ylabel('real axis'), zlabel('imag axis')
view(90,0)
title('Projection onto imaginary and time axes')

% plot real and imaginary projections simultaneously
subplot(224)
plot(time,real(wavelet),'b')
hold on
plot(time,imag(wavelet),'b:')
legend({'real part';'imaginary part'})

%% Figure 13.2

% now show in 3d
figure
plot3(time,real(wavelet),imag(wavelet),'k')
xlabel('Time (ms)'), ylabel('real amplitude'), zlabel('imag amplitude')
title('3-D projection (click and spin with mouse)')
axis equal
rotate3d

%% movie

frequency = 6;         % frequency of the sine wave
srate = 500;           % note: should be the same as the data
time  = -.5:1/srate:.5; % vector of time
 
% make wavelet
wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*(4/(2*pi*frequency))^2));

% make a movie to compare cartesian and polar representation of wavelet
figure, set(gcf,'color',[.6 .2 .7])

timeskip = 1; % if you have a slow computer, set this to, e.g., 5


% setup top row of data (real and imaginary in cartesian plot)
subplot(211)
cplotR = plot(time(1),real(wavelet(1)));
hold on
cplotI = plot(time(1),imag(wavelet(1)),':');
set(gca,'ylim',[-1 1],'xlim',[time(1) time(end)])
xlabel('Time (s)'), ylabel('Amplitude')
title('Cartesian representation')

% setup bottom row of data (polar representation)
subplot(212)
pplot = plot(real(wavelet(1)),imag(wavelet(1)));
set(gca,'ylim',[-1 1],'xlim',[-1 1])
title('Polar representation')
xlabel('Real axis'), ylabel('Imaginary axis')
axis square
 
% loop through time and update data
for ti=1:timeskip:length(time)
    
    % update real part of cartesian plot
    set(cplotR,'XData',time(1:ti),'YData',real(wavelet(1:ti)))
    
    % update imaginary part of cartesian plot
    set(cplotI,'XData',time(1:ti),'YData',imag(wavelet(1:ti)))
    
    % update polar plot
    set(pplot,'XData',real(wavelet(1:ti)),'YData',imag(wavelet(1:ti)))
    drawnow
end

%% Figure 13.4

% Euler's formula: exp(1i*k) gives you a vector on a unit circle with angle k
time = -pi:.25:pi;

figure
subplot(221)
plot(cos(time)+1i*sin(time))
axis square
title('cos\theta + isin\theta')
 
subplot(222)
plot(exp(1i*time))
axis square
title('e^i^\theta')
 
subplot(223)
plot(cos(time)+1i*sin(time),'bo-','markersize',8)
hold on
plot(exp(1i*time),'r.-')
axis square
title('Both')
 
subplot(224)
plot(cos(time)+1i*sin(time))
hold on
plot( (-.35+cos(time)/5) + 1i*(.35+sin(time)/5) ,'m','linew',12) % left eye
plot( (+.35+cos(time)/5) + 1i*(.35+sin(time)/5) ,'r.','markersize',3) % right eye
smile=-pi:.5:0;
plot( (cos(smile)/3) + 1i*(-.35+sin(smile)/5) ,'go','markersize',9) % mouth
xlabel('Real axis')
ylabel('Imaginary axis')
axis square
title('I need a better hobby.')

%% Figure 13.5

% redefine time
time = -.5:1/srate:.5; % vector of time


figure

% plot real and imaginary parts of wavelet
plot(time,real(wavelet),'linew',2)
hold on
plot(time,imag(wavelet),':','linew',2)

% plot cosine and sine
plot(time,cos(2*pi*frequency.*time),'m','linew',2)
plot(time,sin(2*pi*frequency.*time),'m:','linew',2)

% plot gaussian window
gaus_win = exp(-time.^2./(2*s^2));
plot(time,gaus_win,'k')
set(gca,'ylim',[-1.2 1.2])
xlabel('Time (s)')
legend({'real part of wavelet';'imaginary part of wavelet';'cosine';'sine';'Gaussian'})

%% Figure 13.6

load sampleEEGdata % note you don't need the ".mat" (unless the filename contains a period in it)

% create 10 Hz wavelet (kernel)
time = -EEG.pnts/EEG.srate/2 : 1/EEG.srate : EEG.pnts/EEG.srate/2-1/EEG.srate;
f    = 10; % frequency of sine wave in Hz
s    = 4/(2*pi*f);
wavelet = exp(1i*2*pi*f.*time) .* exp(-time.^2./(2*s^2));

% signal is one sine cycle
timeS  = 0:1/EEG.srate:(1/f); % one cycle is 1/f
signal = sin(2*pi*f.*timeS);

% now zero-pad signal
signal = [ zeros(1,EEG.pnts/2-length(timeS)/2) signal zeros(1,EEG.pnts/2-length(timeS)/2) ];

figure

% plot waves
subplot(321)
plot(real(wavelet))
set(gca,'xlim',[200 length(time)-200])

subplot(323)
plot(signal)
set(gca,'xlim',[200 length(time)-200])

subplot(325)
plot(real(conv(wavelet,signal,'same')))
set(gca,'xlim',[200 length(time)-200],'ylim',[-12 12])


% now plot dot products at selected phase lags
subplot(322)
polar(0,12,'-k'), hold on
% compute dot product
dp = sum(wavelet(round(100/f)-2:end).*signal(1:end-round(100/f)+3));
% plot in polar space
polar([angle(dp) angle(dp)],[0 abs(dp)],'k')

subplot(324)
polar(0,12,'-k'), hold on
% compute dot product
dp = sum(wavelet(round(2.3*(100/f)-2):end).*signal(1:end-round(2.3*(100/f)-3)));
% plot in polar space
polar([angle(dp) angle(dp)],[0 abs(dp)],'k')


subplot(326)
polar(0,12,'-k'), hold on
% compute dot product
dp = sum(wavelet.*signal);
% plot in polar space
polar([angle(dp) angle(dp)],[0 abs(dp)],'k')

%% Figure 13.8

% create wavelet
frequency = 6; % in Hz, as usual
time = -1:1/EEG.srate:1;
s    = (4/(2*pi*frequency))^2; % note that s is squared here rather than in the next line...
wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*s)/frequency);

% FFT parameters
n_wavelet            = length(wavelet);
n_data               = EEG.pnts;
n_convolution        = n_wavelet+n_data-1;
half_of_wavelet_size = (length(wavelet)-1)/2;

% FFT of wavelet and EEG data
fft_wavelet = fft(wavelet,n_convolution);
fft_data    = fft(squeeze(EEG.data(47,:,1)),n_convolution); % FCz, trial 1

convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(s);

% cut off edges
convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);

% plot for comparison
figure
subplot(311)
plot(EEG.times,real(convolution_result_fft))
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
title([ 'Projection onto real axis is filtered signal at ' num2str(frequency) ' Hz.' ])

subplot(312)
plot(EEG.times,abs(convolution_result_fft).^2)
xlabel('Time (ms)'), ylabel('Power (\muV^2)')
title([ 'Magnitude of projection vector squared is power at ' num2str(frequency) ' Hz.' ])

subplot(313)
plot(EEG.times,angle(convolution_result_fft))
xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
title([ 'Angle of vector is phase angle time series at ' num2str(frequency) ' Hz.' ])

%% Figure 13.9

figure
plot3(EEG.times,real(convolution_result_fft),imag(convolution_result_fft))
xlabel('Time (ms)'), ylabel('real'), zlabel('imaginary')
grid on
rotate3d

figure
plot3(EEG.times,abs(convolution_result_fft),angle(convolution_result_fft))
title('Click and drag to view phase and amplitude')
xlabel('Time (ms)'), ylabel('amplitude'), zlabel('phase (rad.)')
rotate3d

figure
plot3(EEG.times,angle(convolution_result_fft),abs(convolution_result_fft))
hold on
plot3(EEG.times,angle(convolution_result_fft),real(convolution_result_fft),'r')
title('Click and drag to view phase and amplitude')
xlabel('Time (ms)'), zlabel('amplitude'), ylabel('phase (rad.)')
rotate3d

%% Figure 13.10

srate = 500; % sampling rate in Hz
time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate

% create a 9 Hz wavelet
f = 9; % frequency of wavelet in Hz
s = 6/(2*pi*f);
wavelet9 = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2)); 

% create a 10 Hz wavelet
f = 10; % frequency of wavelet in Hz
s = 6/(2*pi*f);
wavelet10 = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2)); 


figure
subplot(211)
plot(time,real(wavelet9))
hold on
plot(time,real(wavelet10),'r')
xlabel('Time (ms)'), ylabel('Amplitude')

subplot(212)
hz = linspace(0,srate/2,floor(length(time)/2)+1);
fft9  = fft(wavelet9);
fft10 = fft(wavelet10);

plot(hz,abs(fft9(1:length(hz))))
hold on
plot(hz,abs(fft10(1:length(hz))),'r')
set(gca,'xlim',[0 25])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
legend({'9 Hz wavelet';'10 Hz wavelet'})

%% figure 13.11

% definitions, selections...
chan2use = 'fcz';

min_freq =  2;
max_freq = 80;
num_frex = 30;

% define wavelet parameters
time = -1:1/EEG.srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials

baseidx = dsearchn(EEG.times',[-500 -200]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

figure
subplot(121)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title('Logarithmic frequency scaling')

subplot(122)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-200 1000])
title('Linear frequency scaling')

%% IMPORTANT TANGENT HERE ON Y-AXIS SCALING USING IMAGESC!!!

figure
subplot(221)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title('Logarithmic frequency scaling')
ylabel('Frequency (Hz)')

subplot(222)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-200 1000])
title('Linear frequency scaling')

subplot(223)
imagesc(EEG.times,frex,eegpower)
set(gca,'clim',[-3 3],'xlim',[-200 1000],'ydir','norm')
title('WRONG Y-AXIS LABELS!!!!')
ylabel('Frequency (Hz)'), xlabel('Time (ms)')

subplot(224)
imagesc(EEG.times,[],eegpower)
set(gca,'clim',[-3 3],'xlim',[-200 1000],'ydir','norm')
set(gca,'ytick',1:6:num_frex,'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title('CORRECT Y-AXIS LABELS!!!!')
xlabel('Time (ms)')

%% Figure 13.12

frequency = 6;         % frequency of the sine wave
srate = 500;           % note: should be the same as the data
time  = -.5:1/srate:.5; % vector of time
% make wavelet
wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*(4/(2*pi*frequency))^2));

figure
subplot(211)
plot(time,real(wavelet))
title('Good.')

% now make a wavelet that it too short
tooLowFrequency = 2;
wavelet = exp(2*1i*pi*tooLowFrequency.*time) .* exp(-time.^2./(2*(4/(2*pi*tooLowFrequency))^2));

subplot(212)
plot(time,real(wavelet))
xlabel('Time')
title('Not good.')

%% Figure 13.13

frequency = 10;
time      = -.5:1/srate:.5;
numcycles = [ 3 7 ];

wavecolors = 'br';

figure
for i=1:length(numcycles)
    % make wavelet
    wavelet = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*(numcycles(i)/(2*pi*frequency))^2));
    
    subplot(2,2,i)
    plot(time,real(wavelet),wavecolors(i))
    xlabel('Time')
    title([ 'Wavelet at ' num2str(frequency) ' Hz with ' num2str(numcycles(i)) ' cycles' ])
    
    subplot(2,1,2)
    hold on
    fft_wav = 2*abs(fft(wavelet));
    hz_wav  = linspace(0,srate/2,round(length(wavelet)/2)+1);
    plot(hz_wav,fft_wav(1:length(hz_wav)),wavecolors(i))
end

subplot(212)
xlabel('Frequency (Hz)')
set(gca,'xlim',[0 50])
legend({[ num2str(numcycles(1)) ' cycles' ];[ num2str(numcycles(2)) ' cycles' ]})

%% Figure 13.14

% To generate this figure, go up to the code for figure 13.11 and uncomment
% the lines that define the width of the Gaussians for the Morlet wavelets:

% s =  3./(2*pi*frex);
% s = 10./(2*pi*frex);

%% Figure 13.15

frex  = logspace(log10(2),log10(80),30);
srate = 500;
time  = -2:1/srate:2;
N     = length(time);
hz    = linspace(0,srate/2,floor(N/2)+1);
fwhm  = zeros(3,length(frex));

for numcyclesi = 1:3
    
    switch numcyclesi
        case 1
            numcycles=repmat(3,1,length(frex));
        case 2
            numcycles=repmat(10,1,length(frex));
        case 3
            numcycles=logspace(log10(3),log10(10),length(frex));
    end
    
    for fi=1:length(frex)
        
        % make wavelet
        wavelet = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(numcycles(fi)/(2*pi*frex(fi)))^2));
        
        % take FFT of wavelet
        fwave = fft(wavelet);
        fwave = abs(fwave(1:length(hz)))*2;
        
        % normalize power to [0 1]
        fwave = fwave-min(fwave);
        fwave = fwave/max(fwave);
        
        % find left and right 1/2
        [~,peakx]  = max(fwave); % if matlab crashes, replace "~" with "junk"
        [~,left5]  = min(abs(fwave(1:peakx)-.5));
        [~,right5] = min(abs(fwave(peakx:end)-.5));
        right5 = right5+peakx-1;
        
        fwhm(numcyclesi,fi) = hz(right5)-hz(left5);
        
        % plot one example of a wavelet's power spectrum and fwhm
        if fi==ceil(length(frex)/2) && numcyclesi==3
            figure
            
            % plot power spectrum
            plot(hz,fwave,'.-')
            hold on
            
            % plot fwhm
            plot(hz(left5),fwave(left5),'ro','markersize',10)
            plot(hz(right5),fwave(right5),'ro','markersize',10)
            % and draw lines to frequencies
            plot([hz(left5) hz(left5)],[0 fwave(left5)],'r')
            plot([hz(right5) hz(right5)],[0 fwave(right5)],'r')            
            
            set(gca,'xlim',[0 30])
            xlabel('Frequency (Hz)')
            ylabel('Normalized power')
        end
    end
end

figure
plot(frex,fwhm,'.-')
xlabel('Frequency (Hz)')
ylabel('FWHM (Hz)')
legend({'3';'10';'3-10'})

figure
plot(frex,fwhm,'.-')
xlabel('Frequency (Hz)')
ylabel('FWHM (Hz)')
legend({'3';'10';'3-10'})
set(gca,'xlim',[frex(1)*.8 frex(end)*1.2],'ylim',[min(fwhm(:))*.8 max(fwhm(:))*1.2],'xscale','log','xtick',round(logspace(log10(frex(1)),log10(frex(end)),6)),'yscale','log','ytick',round(10*logspace(log10(min(fwhm(:))),log10(max(fwhm(:))),6))/10)

%% end.

