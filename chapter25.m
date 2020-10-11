%% Analyzing Neural Time Series Data
% Matlab code for Chapter 25
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% setup for figure 25.3

% generate a random signal of 3 seconds
srate = 1000;

randsig1 = randn(1,3*srate); % 3 seconds, with this sampling rate
randsig2 = randn(1,3*srate);

% now filter at 5 Hz
f       = 5; % frequency of wavelet in Hz
time    = -1:1/srate:1; % time for wavelet, from -1 to 1 second in steps of 1/sampling-rate
s       = 6/(2*pi*f); % width of Gaussian
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2)); 

% FFT parameters
n_wavelet            = length(wavelet);
n_data               = length(randsig1);
n_convolution        = n_wavelet+n_data-1;
half_of_wavelet_size = (length(wavelet)-1)/2;

% FFT of wavelet and EEG data
convolution_result_fft = ifft(fft(wavelet,n_convolution).*fft(randsig1,n_convolution),n_convolution)*sqrt(s)/10;
filtsig1  = real(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));
anglesig1 = angle(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));

convolution_result_fft = ifft(fft(wavelet,n_convolution).*fft(randsig2,n_convolution),n_convolution)*sqrt(s)/10;
filtsig2  = real(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));
anglesig2 = angle(convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size));


figure
for i=1:2
    subplot(2,1,i)
    eval([ 'plot(randsig' num2str(i) ')' ]) % "eval" can be useful for increasing control over commands
    hold on
    eval([ 'plot(filtsig' num2str(i) ',''r'')' ])
end

%% Figure 25.3

% initialize output correlation matrix
correlations = zeros(5,round(1000/f));

for i=1:round(1000/f)
    
    % correlation of unfiltered random signal
    temp = corrcoef(randsig1(1:end-i),randsig1(i+1:end));
    correlations(1,i) = temp(1,2);
    
    % correlation of filtered signal
    temp = corrcoef(filtsig1(1:end-i),filtsig1(i+1:end));
    correlations(2,i) = temp(1,2);
    
    % phase clustering
    correlations(3,i) = abs(mean(exp(1i*( angle(anglesig1(1:end-i)-anglesig1(i+1:end))))));
    
    % difference of correlations of filtered signal
    temp = corrcoef(filtsig2(1:end-i),filtsig2(i+1:end));
    correlations(4,i) = temp(1,2) - correlations(2,i);
    
    % difference of phase clusterings
    correlations(5,i) = abs(mean(exp(1i*( angle(anglesig2(1:end-i)-anglesig2(i+1:end)))))) - correlations(3,i);
end

figure
plot(correlations')
xlabel('Lag (ms)')
ylabel('Connectivity strength')
legend({'unfiltered';'power corr';'ISPC';'corr diffs';'ISPC diffs'})

%% end.
