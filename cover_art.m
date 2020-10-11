

% load sample EEG data
load sampleEEGdata

% wavelet parameters
min_freq = 2;
max_freq = 48;
num_frex = 30;

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4;


chan2plot = 'fz';

% define baseline period
baselinetime = [ -500 -200 ]; % in ms


% convert baseline window time to indices
[junk,baselineidx(1)]=min(abs(EEG.times-baselinetime(1)));
[junk,baselineidx(2)]=min(abs(EEG.times-baselinetime(2)));

tf_data = zeros(2,length(frequencies),EEG.pnts);


fft_data = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2));
    fft_wavelet = fft(wavelet,n_conv_pow2);
    fft_wavelet = fft_wavelet./max(fft_wavelet);
    
    % run convolution
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
    convolution_result_fft = convolution_result_fft(1:n_convolution);
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    % put power data into time-frequency matrix
    tf_data(1,fi,:) = mean(abs(convolution_result_fft).^2,2);
    tf_data(2,fi,:) = median(abs(convolution_result_fft).^2,2);
end

baseline_power = squeeze(mean(tf_data(1,:,baselineidx(1):baselineidx(2)),3));
dbconverted = 10*log10( squeeze(bsxfun(@rdivide,tf_data(1,:,:),baseline_power) ));

%%

surf(dbconverted(4:end,155:539))
h=findobj('Type','surface');
%set(h,'CData',double(asdf(4:end,155:539)))
set(gca,'clim',[-3 3],'yscale','log'), shading interp, axis off
set(gcf,'color','k')
cmap=(1+[cos(linspace(0,pi*2,100)); sin(linspace(0,pi*2,100)); cos(linspace(0,pi*2,100))])/2;
colormap(cmap')
