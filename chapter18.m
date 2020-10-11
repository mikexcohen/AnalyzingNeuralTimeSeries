%% Analyzing Neural Time Series Data
% Matlab code for Chapter 18
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 18.1

figure

% 1/f function
c = 1;
x = 1;

plot(c./(1:100).^x)

%% Figure 18.2

% load sample EEG data
load sampleEEGdata

% wavelet parameters
min_freq = 2;
max_freq = 128;
num_frex = 30;

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = EEG.pnts;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; 

% FFT of data (note: this doesn't change on frequency iteration)
fft_data = fft(squeeze(EEG.data(23,:,1)),n_conv_pow2);

% initialize output time-frequency data
tf_data = zeros(length(frequencies),EEG.pnts);

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
    fft_wavelet = fft(wavelet,n_conv_pow2);
    
    % run convolution
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
    convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % put power data into time-frequency matrix
    tf_data(fi,:) = abs(convolution_result_fft).^2;
end

% plot results
ytickskip = 2:4:num_frex; % This will be explained in the text.
figure
subplot(221)
imagesc(EEG.times,[],tf_data)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[0 5000])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Color limit of 0 to 5000')

subplot(222)
imagesc(EEG.times,[],tf_data)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[0 800])
title('Color limit of 0 to 800')

subplot(223)
imagesc(EEG.times,[],tf_data)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[0 25])
title('Color limit of 0 to 25')

subplot(224)
imagesc(EEG.times,[],log10(tf_data))
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[-4 4])
title('Color limit of -4 to 4 (log_1_0 units)')

%% Figure 18.3

% define baseline period
baselinetime = [ -500 -200 ]; % in ms

% convert baseline window time to indices
[~,baselineidx(1)]=min(abs(EEG.times-baselinetime(1)));
[~,baselineidx(2)]=min(abs(EEG.times-baselinetime(2)));

% dB-correct
baseline_power = mean(tf_data(:,baselineidx(1):baselineidx(2)),2);
dbconverted = 10*log10( bsxfun(@rdivide,tf_data,baseline_power));
% FYI: the following lines of code are equivalent to the previous line:
% dbconverted = 10*( bsxfun(@minus,log10(tf_data),log10(baseline_power)));
% dbconverted = 10*log10( tf_data ./ repmat(baseline_power,1,EEG.pnts) );
% dbconverted = 10*( log10(tf_data) - log10(repmat(baseline_power,1,EEG.pnts)) );

figure
contourf(EEG.times,frequencies,dbconverted,40,'linecolor','none')
set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',[-500 1500],'clim',[-12 12])
title('Color limit of -12 to +12 dB')

%% Figure 18.4

time2plot = 300; % in ms

[~,timeidx] = min(abs(EEG.times-time2plot));

% plot frequencies
figure
subplot(211)
plot(frequencies,tf_data(:,timeidx))
title([ 'Power spectrum at ' num2str(EEG.times(timeidx)) ' ms' ])
ylabel('Raw power (\muV^2)')
xlabel('Frequency (Hz)')
 
subplot(212)
plot(frequencies,dbconverted(:,timeidx))
ylabel('Baseline-normalized power (dB)')
xlabel('Frequency (Hz)')

%% Figure 18.5 

% This figure was created by changing the color limits of figure 18.3

%% Figure 18.6

activity = 1:.01:20; % activity
baseline = 10; % baseline

db = 10*log10(activity./ baseline );
pc = 100*( activity-baseline)./ baseline;

figure
plot(db,pc)
xlabel('dB'), ylabel('Percent change')

% find indices where db is closest to -/+2
[~,dbOf2]=min(abs(db-2));
[~,dbOfminus2]=min(abs(db--2));

disp([ 'dB of -2 corresponds to ' num2str(pc(dbOfminus2)) '% change.' ])
disp([ 'dB of +2 corresponds to +' num2str(pc(dbOf2)) '% change.' ])

hold on
axislim=axis;
plot([db(dbOf2) db(dbOf2)],[pc(dbOf2) axislim(3)],'k',[axislim(1) db(dbOf2)],[pc(dbOf2) pc(dbOf2)],'k')
plot([db(dbOfminus2) db(dbOfminus2)],[pc(dbOfminus2) axislim(3)],'k',[axislim(1) db(dbOfminus2)],[pc(dbOfminus2) pc(dbOfminus2)],'k')

% real data: percent change vs. baseline division
figure
baseline_power = mean(tf_data(:,baselineidx(1):baselineidx(2)),2);
pctchange = 100 * (tf_data-repmat(baseline_power,1,EEG.pnts))./ repmat(baseline_power,1,EEG.pnts);
subplot(221)
baselinediv = tf_data ./ repmat(baseline_power,1,EEG.pnts);
plot(dbconverted(1:5:end),baselinediv(1:5:end),'.') % don't need all the datapoints to make a point
xlabel('DB'), ylabel('Baseline division')

% dB vs. baseline division
subplot(222)
plot(pctchange(1:5:end),baselinediv(1:5:end),'.')
xlabel('Percent change'), ylabel('Baseline division')

% Z-transform vs. percent change
subplot(223)
baseline_power = tf_data(:,baselineidx(1):baselineidx(2));
baselineZ = (tf_data-repmat(mean(baseline_power,2),1,size(tf_data,2))) ./ repmat(std(baseline_power,[],2),1,size(tf_data,2));
plot(baselineZ(1:5:end),pctchange(1:5:end),'.')
xlabel('Z-transform'), ylabel('Percent change')

% Z-transform vs. dB
subplot(224)
plot(baselineZ(1:5:end),dbconverted(1:5:end),'.')
xlabel('Z-transform'), ylabel('DB')

%% Figure 18.7

figure

% plot dB-converted power
subplot(221)
imagesc(EEG.times,[],dbconverted)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[-10 10])
title('dB change from baseline')

% plot percent-change
subplot(222)
imagesc(EEG.times,[],pctchange)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[-500 500])
title('Percent change from baseline')

% divide by baseline
subplot(223)
imagesc(EEG.times,[],baselinediv)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[-7.5 7.5])
title('Divide by baseline')

% z-transform
subplot(224)
imagesc(EEG.times,[],baselineZ)
set(gca,'ytick',ytickskip,'yticklabel',round(frequencies(ytickskip)),'ydir','normal','xlim',[-500 1500],'clim',[-3.5 3.5])
title('Z-transform')

%% Figures 18.8

chan2plot = 'fcz'; % p1 for figure 18.11
% define baseline period
baselinetime = [ -500 -200 ]; % in ms


% convert baseline window time to indices
[~,baselineidx(1)] = min(abs(EEG.times-baselinetime(1)));
[~,baselineidx(2)] = min(abs(EEG.times-baselinetime(2)));

tf_data = zeros(2,length(frequencies),EEG.pnts);

n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));

fft_data = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
    fft_wavelet = fft(wavelet,n_conv_pow2);
    
    % run convolution
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
    convolution_result_fft = convolution_result_fft(1:n_convolution);
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    
    if fi==6 % save single-trial data from one frequency band
        convdat2keep = convolution_result_fft;
    end
    
    % put power data into time-frequency matrix
    tf_data(1,fi,:) = mean(abs(convolution_result_fft).^2,2);
    tf_data(2,fi,:) = median(abs(convolution_result_fft).^2,2);
end

% db-correct and plot
labelz = {'mean';'median'};
figure
for i=1:2
    baseline_power = squeeze(mean(tf_data(i,:,baselineidx(1):baselineidx(2)),3));
    dbconverted = 10*log10( squeeze(bsxfun(@rdivide,tf_data(i,:,:),baseline_power) ));
    
    % plot
    subplot(2,2,i)
    contourf(EEG.times,frequencies,dbconverted,40,'linecolor','none')
    set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
    title(labelz{i}), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
end

% plot relationship between mean and median
subplot(223)
db_mean = 10*log10( bsxfun(@rdivide,tf_data(1,:,:),mean(tf_data(1,:,baselineidx(1):baselineidx(2)),3)));
db_medn = 10*log10( bsxfun(@rdivide,tf_data(2,:,:),mean(tf_data(2,:,baselineidx(1):baselineidx(2)),3)));
plot(db_mean(:),db_medn(:),'.')
r=corr(db_mean(:),db_medn(:));
legend([ 'R^2 = ' num2str(r*r) ])
xlabel('dB from Mean')
ylabel('dB from Median')

%% Figure 18.9

% plot all trials, mean, and median

figure
subplot(211)
plot(EEG.times,abs(convdat2keep).^2)
set(gca,'xlim',[-200 1000])
subplot(212)
plot(EEG.times,mean(abs(convdat2keep).^2,2))
hold on
plot(EEG.times,median(abs(convdat2keep).^2,2),'r')
set(gca,'xlim',[-200 1000])

% to make this point even more clear, add an outlier trial
convdat2keep(:,100) = convdat2keep(:,10)*100;
figure
subplot(221)
plot(EEG.times,median(abs(convdat2keep).^2,2))
hold on
plot(EEG.times,median(abs(convdat2keep(:,1:end-1)).^2,2),'r')
set(gca,'xlim',[-200 1000])
legend({'With outlier';'Without outlier'})
title('Median: insensitive to outlier trial')

subplot(222)
plot(EEG.times,mean(abs(convdat2keep).^2,2))
hold on
plot(EEG.times,mean(abs(convdat2keep(:,1:end-1)).^2,2),'r')
set(gca,'xlim',[-200 1000])
legend({'With outlier';'Without outlier'})
title('Mean: sensitive to outlier trial')

%% Figure 18.10

% convenientize power
convdatPower  = abs(convdat2keep).^2;

% single-trial linear baseline correction
convdat2keepB = convdatPower - repmat(mean(convdatPower(baselineidx(1):baselineidx(2),:),1),size(convdatPower,1),1);

% single-trial Z score
convdat2keepZ = (convdatPower - repmat(mean(convdatPower(baselineidx(1):baselineidx(2),:),1),size(convdatPower,1),1)) ./ repmat(std(convdatPower(baselineidx(1):baselineidx(2),:)),size(convdatPower,1),1);

% single-trial log10
convdat2keepL = log10(convdatPower);

figure
subplot(221)
plot(EEG.times,mean(convdat2keepB,2),'r')
hold on
plot(EEG.times,median(convdat2keepB,2),'b')
xlabel('Time (ms)'), ylabel('Power')
set(gca,'xlim',[-200 1000])
title('Linear baseline subtraction'),legend({'mean';'median'})

subplot(222)
plot(EEG.times,mean(convdat2keepZ,2),'r')
hold on
plot(EEG.times,median(convdat2keepZ,2),'b')
xlabel('Time (ms)'), ylabel('Power_Z')
set(gca,'xlim',[-200 1000])
title('Z-transformation'),legend({'mean';'median'})

%% Figure 18.11

% This figure was made by running the code for figure 28.8 but using P1
% instead of FCz. (Fewer frequencies were also plotted.)

%% Figure 18.12

snr_bs = zeros(length(frequencies),EEG.pnts);
snr_tf = zeros(length(frequencies),EEG.pnts);
tf     = zeros(length(frequencies),EEG.pnts);

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
    fft_wavelet = fft(wavelet,n_conv_pow2);
    
    % run convolution
    convolution_result = ifft(fft_wavelet.*fft_data,n_conv_pow2) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result = convolution_result(1:n_convolution);
    convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result = reshape(convolution_result,EEG.pnts,EEG.trials);
    
    % extract SNR in two ways
    snr_tf(fi,:) = mean(abs(convolution_result).^2,2)./std(abs(convolution_result).^2,[],2);
    snr_bs(fi,:) = mean(abs(convolution_result).^2,2)./std(mean(abs(convolution_result(baselineidx(1):baselineidx(2),:)).^2,1),[],2);
    
    % and extract trial-averaged power
    tf(fi,:) = mean(abs(convolution_result).^2,2);
    
end

% plot
figure
subplot(121)
contourf(EEG.times,frequencies,snr_bs,40,'linecolor','none')
set(gca,'clim',[.5 2],'xlim',[-200 1000],'ylim',[frequencies(1) 40])
% colorbar
title('SNR_b_a_s_e_l_i_n_e (mean/std)'), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
axis square


subplot(122)

% In the book, I forgot to re-compute baseline-divided power (variable baselinediv), 
% so the variable was from earlier in the script, which is actually just a single
% trial. Thanks to José Luis Ulloa Fulgeri for catching that bug!
baseline_power = mean(tf(:,baselineidx(1):baselineidx(2)),2);
baselinediv    = tf ./ repmat(baseline_power,1,EEG.pnts);

plot(snr_bs(1:3:end),baselinediv(1:3:end),'.')
xlabel('SNR_b_a_s_e_l_i_n_e'), ylabel('Power (/baseline)')
axis square


figure
subplot(121)
contourf(EEG.times,frequencies,snr_tf,40,'linecolor','none')
set(gca,'clim',[.5 1.25],'xlim',[-200 1000],'ylim',[frequencies(1) 40])
% colorbar
title('SNR_t_f (mean/std)'), ylabel('Frequency (Hz)'), xlabel('Time (ms)')
axis square

subplot(122)
plot(snr_tf(1:3:end),baselinediv(1:3:end),'.')
xlabel('SNR_t_f'), ylabel('Power (/baseline)')
axis square

% time-series of SNR
figure
plot(EEG.times,abs(mean(EEG.data(47,:,:),3) ./ std(EEG.data(47,:,:),[],3)))
title('Time-domain SNR time series')
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)'), ylabel('SNR')

% now compute SNR of peak compared to prestim noise
stimeidx = dsearchn(EEG.times',150);
etimeidx = dsearchn(EEG.times',400);
disp([ 'ERP SNR between 150 and 400 ms at FCz: ' num2str(max(mean(EEG.data(47,stimeidx:etimeidx,:),3)) / std(mean(EEG.data(47,baselineidx(1):baselineidx(2),:),3),[],2)) ])

%% Figure 18.13

iterations=10;
chan2plot = 'P7'; % or FCz
dbcorrect = false;

powerByTrialFreq = zeros(length(frequencies),EEG.trials);

start_time = -200; % in ms
end_time   = 1200;

fft_data = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);
timeidx = dsearchn(EEG.times',[start_time end_time]');

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
    fft_wavelet = fft(wavelet,n_conv_pow2);
    
    % run convolution
    convolution_result = ifft(fft_wavelet.*fft_data,n_conv_pow2) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result = convolution_result(1:n_convolution);
    convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result = abs(reshape(convolution_result,EEG.pnts,EEG.trials)).^2; % reshape and convert to power
    
    % "gold standard" is average of all trials
    if dbcorrect
        template = 10*log10( bsxfun(@rdivide,mean(convolution_result,2),mean(mean(convolution_result(baselineidx(1):baselineidx(2),:),1),2)));
        template = template(timeidx(1):timeidx(2));
    else
        template = mean(convolution_result(timeidx(1):timeidx(2),:),2);
    end
    % normalize template for correlation
    template = bsxfun(@rdivide,bsxfun(@minus,template,mean(template)),std(template))';
    
    for iteri=1:iterations
        for triali=5:EEG.trials % start at 5 trials...
            
            trials2use = randsample(1:EEG.trials,triali);
            % if you don't have the stats toolbox, use the following:
            % trials2use = randperm(EEG.trials); trials2use = trials2use(1:triali);
            
            % compute power time series from the random selection of trials, and then normalization
            if dbcorrect
                tempdat = 10*log10( bsxfun(@rdivide,mean(convolution_result(:,trials2use),2),mean(mean(convolution_result(baselineidx(1):baselineidx(2),trials2use),1),2)));
                tempdat = tempdat(timeidx(1):timeidx(2));
            else
                tempdat = mean(convolution_result(timeidx(1):timeidx(2),trials2use),2);
            end
            tempdat = bsxfun(@rdivide,bsxfun(@minus,tempdat,mean(tempdat)),std(tempdat))';
            
            % compute Pearson correlation. This is a super-fast
            % implementation of a Pearson correlation via least squares
            % fit. You'll learn more about this in Chapter 28. 
            powerByTrialFreq(fi,triali) = powerByTrialFreq(fi,triali) + (tempdat*tempdat')\tempdat*template';
        end
    end
end

powerByTrialFreq = powerByTrialFreq./iterations;


figure
plot(5:EEG.trials,squeeze(powerByTrialFreq(:,5:end)))
xlabel('Number of trials'), ylabel('Power')
set(gca,'ylim',[-.1 1.1])
title('Each line is a frequency band')

figure
contourf(5:EEG.trials,frequencies,squeeze(powerByTrialFreq(:,5:end)),40,'linecolor','none')
set(gca,'clim',[.6 1])
xlabel('Number of trials'), ylabel('Frequency (Hz)')
if dbcorrect, title('DB normalized'), else title('not dB normalized'); end

%% end
