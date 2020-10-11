%% Analyzing Neural Time Series Data
% Matlab code for Chapter 23
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% load sample data

load sampleEEGdata

%% Figure 23.2

% compute ERP
erp = squeeze(mean(EEG.data,3));

% subtract mean and compute covariance
erp = bsxfun(@minus,erp,mean(erp,2));
covar = (erp*erp')./(EEG.pnts-1);

figure
subplot(131)
imagesc(covar)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[-1 5])
title('Covariance of ERP')


% one covariance of all timepoints
subplot(132)
eeg = reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials);
eeg = bsxfun(@minus,eeg,mean(eeg,2));
covar = (eeg*eeg')./(length(eeg)-1);
imagesc(covar)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[20 150])
title('Covariance of single-trial EEG')


% average single-trial covariances
subplot(133)
covar = zeros(EEG.nbchan);
% note that the covariance of each trial is computed separately, then averaged
for i=1:EEG.trials
    eeg = bsxfun(@minus,squeeze(EEG.data(:,:,i)),squeeze(mean(EEG.data(:,:,i),2)));
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./i;
imagesc(covar)
axis square
set(gca,'xticklabel',{EEG.chanlocs(get(gca,'xtick')).labels},'yticklabel',{EEG.chanlocs(get(gca,'ytick')).labels},'clim',[20 150])
title('Average covariance of single-trial EEG')

%% Figure 23.3

% compute covariance of ERP
erp = squeeze(mean(EEG.data,3));
erp = bsxfun(@minus,erp,mean(erp,2));
covar = (erp*erp')./(EEG.pnts-1);

% principle components analysis via eigenvalue decomposition
[pc,eigvals] = eig(covar);

% components are listed in increasing order, and converted here to descending order for convenience
pc      = pc(:,end:-1:1);
eigvals = diag(eigvals);
eigvals = 100*eigvals(end:-1:1)./sum(eigvals); % convert to percent change

for i=1:9 % only first 6 are shown in the real figure
    figure(102)
    subplot(3,3,i)
    topoplot(double(pc(:,i)),EEG.chanlocs,'electrodes','off','plotrad',.53);
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals(i)) ])

    figure(101)
    subplot(9,1,i)
    plot(EEG.times,pc(:,i)'*erp)
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    set(gca,'xlim',[-200 1200])
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals(i)) ])
end

% average single-trial covariances
covar = zeros(EEG.nbchan);
for i=1:EEG.trials
    eeg = bsxfun(@minus,squeeze(EEG.data(:,:,i)),squeeze(mean(EEG.data(:,:,i),2)));
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./i;

% principle components analysis via eigenvalue decomposition
[pc,eigvals] = eig(covar);

% components are listed in increasing order, and converted here to descending order for convenience
pc      = pc(:,end:-1:1);
eigvals = diag(eigvals);
eigvals = 100*eigvals(end:-1:1)./sum(eigvals); % convert to percent change

for i=1:9 % only first 6 are shown in the real figure
    figure(10)
    subplot(3,3,i)
    topoplot(double(pc(:,i)),EEG.chanlocs,'electrodes','off','plotrad',.53);
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals(i)) ])

    figure(11)
    subplot(9,1,i)
    
    % PC time course for each trial, then average together
    pctimes = zeros(1,EEG.pnts);
    for triali=1:EEG.trials
        eeg = bsxfun(@minus,squeeze(EEG.data(:,:,triali)),squeeze(mean(EEG.data(:,:,triali),2)));
        pctimes = pctimes + pc(:,i)'*eeg;
    end
    plot(EEG.times,pctimes./EEG.trials)
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    set(gca,'xlim',[-200 1200])
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals(i)) ])
end

%% Tangent...

% An easily made mistake is to confuse the dimension order of the PC matrix. 
% To be sure you have the correct orientation, plot the first component; 
% it should have a spatially broad ERP-like distribution.

figure

subplot(121)
topoplot(double(pc(:,1)),EEG.chanlocs);
title('Correct orientation!')

subplot(122)
topoplot(double(pc(1,:)),EEG.chanlocs);
title('Incorrect orientation!')

%% Figure 23.4

pcanum = 2; % 1 for panel A; 2 for panel B


% average single-trial covariances
covar = zeros(EEG.nbchan);
for i=1:EEG.trials
    eeg = bsxfun(@minus,squeeze(EEG.data(:,:,i)),squeeze(mean(EEG.data(:,:,i),2)));
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./i;

% principle components analysis via eigenvalue decomposition
[pc,eigvals] = eig(covar);

% components are listed in increasing order, and converted here to descending order for convenience
pc      = pc(:,end:-1:1);
eigvals = diag(eigvals);
eigvals = 100*eigvals(end:-1:1)./sum(eigvals); % convert to percent change

pcadata = zeros(EEG.pnts,EEG.trials);
for i=1:EEG.trials
    pcadata(:,i) = pc(:,pcanum)'*bsxfun(@minus,squeeze(EEG.data(:,:,i)),squeeze(mean(EEG.data(:,:,i),2)));
end



min_freq =  2;
max_freq = 80;
num_frex = 30;

% define wavelet parameters
time = -1:1/EEG.srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(pcadata,1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials

baseidx = dsearchn(EEG.times',[-500 -200]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

figure
subplot(221)
topoplot(double(pc(:,pcanum)),EEG.chanlocs,'plotrad',.53);
subplot(222)
plot(EEG.times,mean(pcadata,2))
set(gca,'xlim',[-200 1200])

subplot(212)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title([ 'TF power from component ' num2str(pcanum) ])

%% Figure 23.5

figure

% recompute PCA
covar = zeros(EEG.nbchan);
for i=1:EEG.trials
    eeg = bsxfun(@minus,squeeze(EEG.data(:,:,i)),squeeze(mean(EEG.data(:,:,i),2)));
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./i;
[pc,eigvals] = eig(covar);
eigvals = diag(eigvals);
eigvals = 100*eigvals(end:-1:1)./sum(eigvals); % convert to percent change



% plot eigenvalues as percent variance accounted for
plot(eigvals,'-o','markerface','w')
set(gca,'ylim',[-1 70],'xlim',[0 65])


% amount of variance expected by chance, computed analytically
hold on
plot([1 EEG.nbchan],repmat(100/EEG.nbchan,1,2),'k')


% amount of variance expected by chance, computed based on random data
nperms = 1000;
permeigvals = zeros(nperms,EEG.nbchan);
for permi=1:nperms
    % random data by randomizing ERP time points/channels
    randdat = reshape(erp(randperm(numel(erp))),size(erp));
    covar = (randdat*randdat')./(EEG.pnts-1);
    [pc,eigvals] = eig(covar);
    eigvals = diag(eigvals);
    permeigvals(permi,:) = 100*eigvals(end:-1:1)./sum(eigvals);
end

plot(mean(permeigvals,1),'r-o','markerface','w')

legend({'% var. accounted for';'chance-level (alg)';'chance-level (perm.test)'})

%% Figure 23.6

whichcomp = 1; % 1 for panel A; 2 for panel B


centertimes = -200:50:1200;
timewindow  = 200; % ms on either side of center times
if whichcomp==1
    maptimes = [ -100 200 500 1000 ]; % times for plotting topomaps, based on visual inspection of PCA-coherence time courses
    clim     = [.08 .15]; % color limits also based on visual inspection
else
    maptimes = [ 0 300 750 1000 ];
    clim     = [-.2 .2]; % color limits also based on visual inspection
end



pcvariance = zeros(size(centertimes));
firstpcas  = zeros(length(centertimes),EEG.nbchan);

timesidx   = dsearchn(EEG.times',centertimes');
timewinidx = round(timewindow/(1000/EEG.srate));
mapsidx    = dsearchn(centertimes',maptimes');


for ti=1:length(centertimes)
    
    % temporally localized covariance
    covar = zeros(EEG.nbchan);
    for i=1:EEG.trials
        eeg = squeeze(EEG.data(:,timesidx(ti)-timewinidx:timesidx(ti)+timewinidx,i));
        eeg = bsxfun(@minus,eeg,mean(eeg,2));
        covar = covar + (eeg*eeg')./(EEG.pnts-1);
    end
    covar = covar./i;
    
    % principle components analysis via eigenvalue decomposition
    [pc,eigvals] = eig(covar);
    pc      = pc(:,end:-1:1);
    eigvals = diag(eigvals);
    eigvals = 100*eigvals(end:-1:1)./sum(eigvals); % convert to percent change
    
    pcvariance(ti)  = eigvals(whichcomp);
    firstpcas(ti,:) = pc(:,whichcomp);
end

figure
plot(centertimes,pcvariance)
xlabel('Time (ms)'), ylabel([ '% variance from PC' num2str(whichcomp) ])


figure
for i=1:length(maptimes)
    subplot(ceil(length(maptimes)/ceil(sqrt(length(maptimes)))),ceil(sqrt(length(maptimes))),i)
    topoplot(firstpcas(mapsidx(i),:),EEG.chanlocs,'plotrad',.53,'maplimits',clim);
    title([ 'PC' num2str(whichcomp) ' from ' num2str(maptimes(i)) ' ms' ]);
end

%% Figure 23.7

% filter data
center_freq = 12; % in Hz
filter_frequency_spread  = 3; % Hz +/- the center frequency
transition_width  = 0.2;

% construct filter kernel
nyquist       = EEG.srate/2;
filter_order  = round(3*(EEG.srate/(center_freq-filter_frequency_spread)));

ffrequencies  = [ 0 (1-transition_width)*(center_freq-filter_frequency_spread) (center_freq-filter_frequency_spread) (center_freq+filter_frequency_spread) (1+transition_width)*(center_freq+filter_frequency_spread) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

filter_result = filtfilt(filterweights,1,double(reshape(EEG.data,EEG.nbchan,EEG.pnts*EEG.trials))')';
filter_result = reshape(filter_result,EEG.nbchan,EEG.pnts,EEG.trials);


% average single-trial covariances
covar = zeros(EEG.nbchan);
for i=1:EEG.trials
    eeg = squeeze(filter_result(:,:,i)) - repmat(squeeze(mean(filter_result(:,:,i),2)),1,EEG.pnts);
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./i;

% principle components analysis via eigenvalue decomposition
[pc,eigvals] = eig(covar);

% components are listed in increasing order, and converted here to descending order for convenience
pc      = pc(:,end:-1:1);
eigvals = diag(eigvals);
eigvals = 100*eigvals(end:-1:1)./sum(eigvals); % convert to percent change

for i=1:9 % only first 6 are shown in the real figure
    figure(10)
    subplot(3,3,i)
    topoplot(pc(:,i),EEG.chanlocs,'electrodes','off','plotrad',.53);
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals(i)) ])

    figure(11)
    subplot(9,1,i)
    plot(EEG.times,pc(:,i)'*squeeze(mean(filter_result,3)))
    hold on
    plot(get(gca,'xlim'),[0 0],'k')
    set(gca,'xlim',[-200 1200])
    title([ 'PC #' num2str(i) ', eigval=' num2str(eigvals(i)) ])
end

%% Figure 23.8

% Note about this code: The legend of this figure (page 304) states that the PCA computed on all trials and then the weights were
% separately applied to the first and last 30 trials. However, in the code below (and thus, in the book figure), the PCA is 
% actually computed separately on the first and last 30 trials. Thanks to Matt Mollison for bringing this mistake to my attention.

figure

% PCA on first 30 trials
covar = zeros(EEG.nbchan);
for i=1:30
    eeg = squeeze(filter_result(:,:,i)) - repmat(squeeze(mean(filter_result(:,:,i),2)),1,EEG.pnts);
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./i;

% principle components analysis via eigenvalue decomposition
[pc,eigvals] = eig(covar);
pc = pc(:,end:-1:1);

for i=1:9 % only first 6 are shown in the real figure
    subplot(3,3,i)
    plot(EEG.times,pc(:,i)'*squeeze(mean(filter_result(:,:,1:30),3)))
    hold on
    set(gca,'xlim',[-200 1200])
    title([ 'PC #' num2str(i) ])
end


% PCA on last 30 trials
covar = zeros(EEG.nbchan);
for i=EEG.trials-29:EEG.trials
    eeg = squeeze(filter_result(:,:,i)) - repmat(squeeze(mean(filter_result(:,:,i),2)),1,EEG.pnts);
    covar = covar + (eeg*eeg')./(EEG.pnts-1);
end
covar = covar./30;

% principle components analysis via eigenvalue decomposition
[pc,eigvals] = eig(covar);
pc = pc(:,end:-1:1);

for i=1:9
    subplot(3,3,i)
    plot(EEG.times,pc(:,i)'*squeeze(mean(filter_result(:,:,end-30:end),3)),'r')
    set(gca,'xlim',[-200 1200])
    title([ 'PC #' num2str(i) ])
end
legend({'first 30 trials';'last 30 trials'})

%% end.
