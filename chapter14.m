%% Analyzing Neural Time Series Data
% Matlab code for Chapter 14
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 14.1

% create cosine
time   = 0:.001:1;
cosine = cos(2*pi*5*time);

figure
plot(time,cosine)
hold on
plot(time(1:5:end),real(hilbert(cosine(1:5:end))),'ko') % plot every 5th point because of overlap
plot(time,imag(hilbert(cosine)),'r')
plot(time,angle(hilbert(cosine)),'m')

ylabel('Angle or amplitude')
legend({'cosine';'real part';'imag part';'angle'})

%% the FFT-based hilbert transform

% generate random numbers
n = 21;
randomnumbers = randn(n,1);

% take FFT
f = fft(randomnumbers);
% create a copy that is multiplied by the complex operator
complexf = 1i*f;

% find indices of positive and negative frequencies
posF = 2:floor(n/2)+mod(n,2);
negF = ceil(n/2)+1+~mod(n,2):n;

% rotate Fourier coefficients
% (note 1: this works by computing the iAsin(2pft) component, i.e., the phase quadrature)
% (note 2: positive frequencies are rotated counter-clockwise; negative frequencies are rotated clockwise)
f(posF) = f(posF) + -1i*complexf(posF);
f(negF) = f(negF) +  1i*complexf(negF);
% The next two lines are an alternative and slightly faster method. 
% The book explains why this is equivalent to the previous two lines.
% f(posF) = f(posF)*2;
% f(negF) = f(negF)*0;

% take inverse FFT
hilbertx = ifft(f);

% compare with Matlab function hilbert
hilbertm = hilbert(randomnumbers);

% plot results
figure
subplot(211)
plot(abs(hilbertm))
hold on
plot(abs(hilbertx),'ro')
legend({'Matlab Hilbert function';'"manual" Hilbert'})
title('magnitude of Hilbert transform')

subplot(212)
plot(angle(hilbertm))
hold on
plot(angle(hilbertx),'ro')
legend({'Matlab Hilbert function';'"manual" Hilbert'})
title('phase of Hilbert transform')

%% Figure 14.2

load sampleEEGdata

% first, filter data (filter mechanisms will be explained more below; for now, focus on 
% using the phases from the Hilbert transform to test whether the matrix input was correct)

nyquist = EEG.srate/2;
lower_filter_bound = 4; % Hz
upper_filter_bound = 10; % Hz
transition_width   = 0.2;
filter_order       = round(3*(EEG.srate/lower_filter_bound));

% create the filter shape (this is explained more in the text around figure 14.4)
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

% apply the filter kernal to the data to obtain the band-pass filtered signal
filtered_data = zeros(EEG.nbchan,EEG.pnts);
for chani=1:EEG.nbchan
    filtered_data(chani,:) = filtfilt(filterweights,1,double(EEG.data(chani,:,1)));
end

% apply hilbert transform in correct and incorrect orientations
hilbert_oops = hilbert(filtered_data);
hilbert_yes  = hilbert(filtered_data')'; % time should be in the first dimension. 
% Note that the output of the hilbert transform is transposed to bring us back to an electrode X time matrix.

figure
subplot(221)
plot(EEG.times,angle(hilbert_yes(1,:)'),'b');
title('correct matrix orientation')
xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
set(gca,'xlim',[-1000 1500])

subplot(222)
plot(EEG.times,angle(hilbert_oops(1,:)),'r');
title('incorrect matrix orientation')
xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
set(gca,'xlim',[-1000 1500])

subplot(223)
plot(EEG.times,real(hilbert_yes(1,:)),'b');
title('correct matrix orientation')
xlabel('Time (ms)'), ylabel('Amplitude')
set(gca,'xlim',[-1000 1500])

subplot(224)
plot(EEG.times,real(hilbert_oops(1,:)),'r');
title('incorrect matrix orientation')
xlabel('Time (ms)'), ylabel('Amplitude')
set(gca,'xlim',[-1000 1500])

%% Figure 14.3

center_freq = 20; % in Hz
filter_frequency_spread  = 6; % Hz +/- the center frequency
wavelet_frequency_spread = 4; % number of wavelet cycles

% create wavelet...
time = -1000/EEG.srate/10:1/EEG.srate:1000/EEG.srate/10;
wavelet = zscore(exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(wavelet_frequency_spread/(2*pi*center_freq))^2)));
% ... and compute its power spectrum
fft_wavelet = abs(fft(wavelet));
fft_wavelet = fft_wavelet./max(fft_wavelet); % normalized to 1.0 for visual comparison ease
hz_wavelet  = linspace(0,nyquist,length(time)/2+1);


% construct filter kernel (the mechanics of filter construction will be
% discussed in the text around figure 14.4)
transition_width = 0.2;
ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread) (center_freq-filter_frequency_spread) (center_freq+filter_frequency_spread) (1+transition_width)*(center_freq+filter_frequency_spread) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
filterweights  = zscore(firls(200,ffrequencies,idealresponse));
% also compute weights using fir1
filterweights1 = zscore(fir1(200,[center_freq-filter_frequency_spread center_freq+filter_frequency_spread]./nyquist));
% compute its power spectrum
fft_filtkern  = abs(fft(filterweights));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
fft_filtkern1 = abs(fft(filterweights1));
fft_filtkern1 = fft_filtkern1./max(fft_filtkern1); % normalized to 1.0 for visual comparison ease

hz_filtkern   = linspace(0,nyquist,101); % list of frequencies in Hz corresponding to filter kernel


figure
subplot(211)

% plot wavelet and filter kernel
plot(real(wavelet))
hold on
plot(filterweights,'r')
legend({'wavelet';'filter kernel'})
set(gca,'xlim',[0 200],'ylim',[-5 5])
xlabel('Time')

% plot power spectra
subplot(212)
plot(hz_wavelet,fft_wavelet(1:ceil(length(fft_wavelet)/2)))
hold on
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'r')
legend({'wavelet';'filter kernel'})
set(gca,'ylim',[-.1 1.1])
xlabel('Frequency (Hz)')

%% Figure 14.4

center_freq = 20; % in Hz
filter_frequency_spread_wide = 10; % Hz +/- the center frequency

ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_wide) (center_freq-filter_frequency_spread_wide) (center_freq+filter_frequency_spread_wide) (1+transition_width)*(center_freq+filter_frequency_spread_wide) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
filterweightsW = zscore(firls(200,ffrequencies,idealresponse));

figure
plot(ffrequencies*nyquist,idealresponse,'r')
hold on

fft_filtkern  = abs(fft(filterweightsW));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')

set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
legend({'ideal';'best fit'})

%% Figure 14.5

center_freq = 20; % in Hz
filter_frequency_spread_wide = 10; % Hz +/- the center frequency
filter_frequency_spread_naro =  2; % Hz +/- the center frequency


% construct filter kernels
ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_wide) (center_freq-filter_frequency_spread_wide) (center_freq+filter_frequency_spread_wide) (1+transition_width)*(center_freq+filter_frequency_spread_wide) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
filterweightsW = zscore(firls(200,ffrequencies,idealresponse));

ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_naro) (center_freq-filter_frequency_spread_naro) (center_freq+filter_frequency_spread_naro) (1+transition_width)*(center_freq+filter_frequency_spread_naro) nyquist ]/nyquist;
filterweightsN = zscore(firls(200,ffrequencies,idealresponse));


figure
subplot(211)
fft_filtkern  = abs(fft(filterweightsW));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)))

hold on
fft_filtkern  = abs(fft(filterweightsN));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'r')
set(gca,'ylim',[-.1 1.1])
legend({'10 Hz width';'2 Hz width'})

subplot(223)
plot(filterweightsW)
hold on
plot(filterweightsN,'r')
set(gca,'xlim',[0 200],'ylim',[-4 7])
legend({'10 Hz spread';'2 Hz spread'})
xlabel('Time')

subplot(224)
plot(abs(hilbert(filterweightsW)),'b')
hold on
plot(abs(hilbert(filterweightsN)),'r')
set(gca,'xlim',[0 200],'ylim',[-4 7])
legend({'10 Hz spread';'2 Hz spread'})
xlabel('Time')

%% Figure 14.6

figure
subplot(211)

% plot wavelet and filter kernel
plot(filterweights1)
hold on
plot(filterweights,'r')
legend({'fir1';'firls'})
set(gca,'xlim',[0 200])
xlabel('Time')

% plot power spectra
subplot(212)
plot(hz_filtkern,fft_filtkern1(1:ceil(length(fft_filtkern1)/2)))
hold on
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'r')
legend({'fir1';'firls'})
set(gca,'ylim',[-.1 1.1])
xlabel('Frequency (Hz)')

% This next figure shows that the filter kernel computed by fir1 is the same as
% that produced by firls if you set the transition zone to zero and then
% smooth the filter kernel with a Hamming window. In the plot, the red and magenta
% lines overlap, which is why you don't see the fir1 filter kernel. You can
% subtract them and show that the difference (which is due to re-scaling 
% of the kernel after windowing) is 3-4 orders of magnitude smaller than
% the kernel itself, hence, nearly identical.

freqL = center_freq-filter_frequency_spread;
freqU = center_freq+filter_frequency_spread;

ffrequencies   = [ 0 freqL freqL freqU freqU nyquist ]/nyquist; % transition zone of 0
filterweights  = firls(200,ffrequencies,idealresponse);
filterweights1 = fir1(200,[freqL freqU]./nyquist);

figure
subplot(211)
plot(filterweights)
hold on
plot(filterweights1,'r.')
plot(filterweights.*hamming(length(filterweights))','m')
legend({'firls';'fir1';'firls with hamming window'})

%% Figure 14.7a

center_freq = 60; % in Hz
filter_frequency_spread_wide = 10; % Hz +/- the center frequency

ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_wide) (center_freq-filter_frequency_spread_wide) (center_freq+filter_frequency_spread_wide) (1+transition_width)*(center_freq+filter_frequency_spread_wide) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
filterweightsW = zscore(firls(200,ffrequencies,idealresponse));

figure
plot(ffrequencies*nyquist,idealresponse,'r')
hold on

fft_filtkern  = abs(fft(filterweightsW));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')

set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
legend({'ideal';'best fit'})

freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
title([ 'SSE: ' num2str(sum( (idealresponse-fft_filtkern(freqsidx)).^2 )) ])




center_freq = 60; % in Hz
filter_frequency_spread_wide = 15; % Hz +/- the center frequency

ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_wide) (center_freq-filter_frequency_spread_wide) (center_freq+filter_frequency_spread_wide) (1+transition_width)*(center_freq+filter_frequency_spread_wide) nyquist ]/nyquist;
idealresponse  = [ 0 0 1 1 0 0 ];
filterweightsW = zscore(firls(200,ffrequencies,idealresponse));

figure
plot(ffrequencies*nyquist,idealresponse,'r')
hold on

fft_filtkern  = abs(fft(filterweightsW));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')

set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
legend({'ideal';'best fit'})

freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
title([ 'SSE: ' num2str(sum( (idealresponse-fft_filtkern(freqsidx)).^2 )) ])

%% Figure 14.7b (note: takes a while to run...)

centerfreqs  = linspace(5,80,60);
transwidths  = linspace(0.01,0.2,40);
filterwidths = linspace(0.05,0.3,40);


sse = zeros(length(centerfreqs),length(transwidths));
for centfreqi = 1:length(centerfreqs)
    
    center_freq = centerfreqs(centfreqi);
    
    for transwidi = 1:length(transwidths)
        
        filter_frequency_spread_wide = center_freq*.2;
        transition_width = transwidths(transwidi);
        
        ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_wide) (center_freq-filter_frequency_spread_wide) (center_freq+filter_frequency_spread_wide) (1+transition_width)*(center_freq+filter_frequency_spread_wide) nyquist ]/nyquist;
        filterweightsW = zscore(firls(200,ffrequencies,idealresponse));
        
        fft_filtkern   = abs(fft(filterweightsW));
        fft_filtkern   = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
        freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
        
        sse(centfreqi,transwidi) = sum( (idealresponse-fft_filtkern(freqsidx)).^2 );
    end
end

figure
contourf(transwidths,centerfreqs,sse,40,'linecolor','none')
xlabel('transwidths'), ylabel('center frequencies')
set(gca,'clim',[0 1]), colorbar


sse = zeros(length(centerfreqs),length(transwidths));
for centfreqi = 1:length(centerfreqs)
    
    center_freq = centerfreqs(centfreqi);
    
    for transwidi = 1:length(transwidths)
        
        filter_frequency_spread_wide = center_freq*filterwidths(transwidi);
        transition_width = .2;
        
        ffrequencies   = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread_wide) (center_freq-filter_frequency_spread_wide) (center_freq+filter_frequency_spread_wide) (1+transition_width)*(center_freq+filter_frequency_spread_wide) nyquist ]/nyquist;
        filterweightsW = zscore(firls(200,ffrequencies,idealresponse));
        
        fft_filtkern  = abs(fft(filterweightsW));
        fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
        freqsidx = dsearchn(hz_filtkern',ffrequencies'*nyquist);
        
        sse(centfreqi,transwidi) = sum( (idealresponse-fft_filtkern(freqsidx)).^2 );
    end
end

figure
contourf(filterwidths,centerfreqs,sse,40,'linecolor','none')
xlabel('filterwidths (% center freq.)'), ylabel('center frequencies')
set(gca,'clim',[0 1]), colorbar

%% Figure 14.8

center_freq   = 10;
freqspread    = 4; % Hz +/- the center frequency
transwid      = .10;
ffrequencies  = [ 0 (1-transwid)*(center_freq-freqspread) (center_freq-freqspread) (center_freq+freqspread) (1+transwid)*(center_freq+freqspread) nyquist ]/nyquist;

data2filter   = squeeze(double(EEG.data(47,:,1)));
filterweights = firls(200,ffrequencies,idealresponse); % recompute without z-scoring

filter_result = filtfilt(filterweights,1,data2filter);
convol_result = conv(data2filter,filterweights,'same'); % could also use ifft(fft(data2filter...

figure
plot(EEG.times,filter_result)
hold on
plot(EEG.times,convol_result,'r')
set(gca,'xlim',[-200 1000]) % zoom-in
xlabel('Time (ms)'), ylabel('Amplitude (\muV)');
legend({'filtfilt';'conv'})

%% Figure 14.9

center_freq = 10; % in Hz
nyquist     = EEG.srate/2;

% create short sine wave
time    = -2000/EEG.srate/10:1/EEG.srate:2000/EEG.srate/10;
wavelet = cos(2*pi*center_freq.*time) .* exp(-time.^2./(2*(3/(2*pi*center_freq))^2));


freqspread = 4; % Hz +/- the center frequency
transwid   = .10;

% construct filter kernels
ffrequencies  = [ 0 (1-transwid)*(center_freq-freqspread) (center_freq-freqspread) (center_freq+freqspread) (1+transwid)*(center_freq+freqspread) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(250,ffrequencies,idealresponse);


forward_filt = filter(filterweights,1,wavelet);
reverse_filt = filter(filterweights,1,forward_filt(end:-1:1));
final_filt_result = reverse_filt(end:-1:1); % must reverse time again!

figure
plot(wavelet)
hold on
plot(forward_filt,'r')
plot(reverse_filt,'m')

set(gca,'xlim',[0 length(wavelet)])
legend({'signal';'forward filter';'reverse filter'})

%% Figure 14.10

% 5th-order butterworth filter
[butterB,butterA] = butter(5,[(center_freq-filter_frequency_spread)/nyquist (center_freq+filter_frequency_spread)/nyquist],'bandpass');
butter_filter     = filtfilt(butterB,butterA,data2filter);

figure
subplot(211)

% plot real part of the filtered signal
plot(EEG.times,filter_result)
hold on
plot(EEG.times,butter_filter,'r')
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)'), ylabel('Amplitude (\muV)');
legend({'FIR';'Butterworth'})

% now plot phases
subplot(212)
plot(EEG.times,angle(hilbert(filter_result)))
hold on
plot(EEG.times,angle(hilbert(butter_filter)),'r')
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)'), ylabel('Phase angles (rad.)');
legend({'FIR';'Butterworth'})

%% Figure 14.12

elap_time = [0 0];
num_iter  = 100;


freqspread  =  4; % Hz +/- the center frequency
center_freq = 20;
transwid    = .15;

% construct filter kernels
ffrequencies  = [ 0 (1-transwid)*(center_freq-freqspread) (center_freq-freqspread) (center_freq+freqspread) (1+transwid)*(center_freq+freqspread) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(3*round(EEG.srate/(center_freq-freqspread)),ffrequencies,idealresponse);


for i=1:num_iter
    tic
    data2filter_cat = squeeze(double(reshape(EEG.data(47,:,:),1,EEG.pnts*EEG.trials)));
    filtdat_cat = reshape(filtfilt(filterweights,1,data2filter_cat),EEG.pnts,EEG.trials);
    elap_time(1) = elap_time(1) + toc;
end

for i=1:num_iter
    tic
    data2filter_sep = squeeze(double(EEG.data(47,:,:)));
    filtdat_sep = zeros(size(data2filter_sep));
    for triali=1:EEG.trials
        filtdat_sep(:,triali) = filtfilt(filterweights,1,data2filter_sep(:,triali));
    end
    elap_time(2) = elap_time(2) + toc;
end

elap_time = elap_time/num_iter;

% plot
figure
plot(EEG.times,mean(filtdat_cat,2))
hold on
plot(EEG.times,mean(filtdat_sep,2),'r')
legend({'concatenated';'separated'})

figure
bar(elap_time)
set(gca,'xlim',[0 3],'xticklabel',{'Concatenated';'Separated'})
ylabel('Time (s)')
title([ 'Speed increase of ' num2str(elap_time(2)/elap_time(1)) '!' ])

%% end.
