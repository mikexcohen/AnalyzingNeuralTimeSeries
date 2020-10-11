%% Analyzing Neural Time Series Data
% Matlab code for Chapter 34
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% load sample data

load sampleEEGdata

% note: Most of these figures take a while to generate. Have patience!

%% extract TF power (create data that are used for the rest of this chapter)

% definitions, selections...
chan2use = 'fcz';

min_freq = 3;
max_freq = 30;
num_frex = 20;


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

% note that you don't need the wavelet itself, you need the FFT of the wavelet
wavelets = zeros(num_frex,n_conv_pow2);
for fi = 1:num_frex
    wavelets(fi,:)  = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
end

% get FFT of data
eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts,EEG.trials); % frequencies X time X trials
eegphase = zeros(num_frex,EEG.pnts,EEG.trials); % frequencies X time X trials

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    % convolution
    eegconv = ifft(wavelets(fi,:).*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % reshape to time X trials
    eegpower(fi,:,:) = abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2;
    eegphase(fi,:,:) = exp(1i*angle(reshape(eegconv,EEG.pnts,EEG.trials)));
end

% remove edge artifacts
time_s = dsearchn(EEG.times',-500);
time_e = dsearchn(EEG.times',1200);

eegpower = eegpower(:,time_s:time_e,:);
tftimes  = EEG.times(time_s:time_e);
nTimepoints = numel(tftimes);

%% Figure 34.1

voxel_pval   = 0.01;
cluster_pval = 0.05;

% note: try to use 1000 or more permutations for real data
n_permutes = 1000;

baseidx(1) = dsearchn(tftimes',-500);
baseidx(2) = dsearchn(tftimes',-100);

% compute actual t-test of difference
realbaselines = squeeze(mean(eegpower(:,baseidx(1):baseidx(2),:),2));
realmean      = 10*log10(bsxfun(@rdivide, mean(eegpower,3), mean(realbaselines,2)));

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_frex);
permuted_vals    = zeros(n_permutes,num_frex,numel(tftimes));
max_clust_info   = zeros(n_permutes,1);


for permi=1:n_permutes
    cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1);
    permuted_vals(permi,:,:) = 10*log10(bsxfun(@rdivide,mean(eegpower(:,[cutpoint:end 1:cutpoint-1],:),3),mean(realbaselines,2)) );
    % btw, using bsxfun instead of repmat increases the speed 
    % of this loop by a factor of ~5.
end

zmap = (realmean-squeeze(mean(permuted_vals))) ./ squeeze(std(permuted_vals));
threshmean = realmean;
threshmean(abs(zmap)<norminv(1-voxel_pval))=0;

figure
subplot(221)
contourf(tftimes,frex,realmean,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('power map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

subplot(222)
contourf(tftimes,frex,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

subplot(223)
contourf(tftimes,frex,threshmean,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('Uncorrected power map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')



% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    fakecorrsz = squeeze((permuted_vals(permi,:,:)-mean(permuted_vals,1)) ./ std(permuted_vals,[],1) );
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(fakecorrsz);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    % using cellfun here eliminates the need for a slower loop over cells
end

% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

subplot(224)
contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('Cluster-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%% Figure 34.3

voxel_pval = 0.05;
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
mcc_cluster_pval = 0.05;


% note: try to use 1000 or more permutations for real data
n_permutes = 1000;

real_condition_mapping = [ -ones(1,floor(EEG.trials/2)) ones(1,ceil(EEG.trials/2)) ];

% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(mean(eegpower(:,:,real_condition_mapping==-1),3) - mean(eegpower(:,:,real_condition_mapping==1),3));
tdenom = sqrt( (std(eegpower(:,:,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) + (std(eegpower(:,:,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1) );
real_t = tnum./tdenom;


% initialize null hypothesis matrices
permuted_tvals  = zeros(n_permutes,num_frex,nTimepoints);
max_pixel_pvals = zeros(n_permutes,2);
max_clust_info  = zeros(n_permutes,1);

% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    fake_condition_mapping = sign(randn(EEG.trials,1));
    
    % compute t-map of null hypothesis
    tnum   = squeeze(mean(eegpower(:,:,fake_condition_mapping==-1),3)-mean(eegpower(:,:,fake_condition_mapping==1),3));
    tdenom = sqrt( (std(eegpower(:,:,fake_condition_mapping==-1),0,3).^2)./sum(fake_condition_mapping==-1) + (std(eegpower(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) );
    tmap   = tnum./tdenom;
    
    % save all permuted values
    permuted_tvals(permi,:,:) = tmap;
    
    % save maximum pixel values
    max_pixel_pvals(permi,:) = [ min(tmap(:)) max(tmap(:)) ];
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,EEG.trials-1))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(tmap);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
end

% now compute Z-map
zmap = (real_t-squeeze(mean(permuted_tvals,1)))./squeeze(std(permuted_tvals));

figure
subplot(221)
contourf(tftimes,frex,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('Unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% apply uncorrected threshold
subplot(222)
contourf(tftimes,frex,zmap,40,'linecolor','none')
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false;
zmapthresh=logical(zmapthresh);
hold on
contour(tftimes,frex,zmapthresh,1,'linecolor','k')

axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('Unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')



% apply pixel-level corrected threshold
lower_threshold = prctile(max_pixel_pvals(:,1),    mcc_voxel_pval*100/2);
upper_threshold = prctile(max_pixel_pvals(:,2),100-mcc_voxel_pval*100/2);

zmapthresh = zmap;
zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=0;
subplot(223)
contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('Pixel-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

subplot(224)
contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[-500 1200])
title('Cluster-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%% Figure 34.4

voxel_pval = 0.01;
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
mcc_cluster_pval = 0.05;


% note: try to use 1000 or more permutations for real data
n_permutes = 1000;

rts=zeros(1,EEG.trials);
for ei=1:EEG.trials
    % In this task, the button press always followed the stimulus at
    % time=0. Thus, finding the RT involves finding the latency of the
    % event that occurs after the time=0 event.
    % If you follow a procedure like this in your data, you may need to
    % include special exceptions, e.g., if there was no response or if
    % a non-response marker could have occurred between stimulus and response.
    time0event = find(cell2mat(EEG.epoch(ei).eventlatency)==0);
    rts(ei)    = EEG.epoch(ei).eventlatency{time0event+1};
end

% rank-transform RTs
rtsrank = tiedrank(rts);

% rank-transform power data (must be transformed)
eegpowerreshaped = reshape(eegpower,num_frex*nTimepoints,EEG.trials)';
eegpowerrank = tiedrank(eegpowerreshaped);

% technically, you want to perform a correlation, but a linear least-squares fit provides the same conceptual results 
% while being ~15 times faster, and we don't care here about the actual scale of the data. For completeness,
% the following line shows you how to compute the Spearman correlation coefficient
% realcorrs = 1-6*sum((eegpowerrank-repmat(rtsrank',1,size(eegpowerrank,2))).^2)/(EEG.trials*(EEG.trials^2-1));
realcorrs = (rtsrank*rtsrank')\rtsrank*eegpowerrank;
realcorrs = reshape(realcorrs,num_frex,nTimepoints);

% initialize null hypothesis matrices
permuted_rvals  = zeros(n_permutes,num_frex,nTimepoints);
max_pixel_rvals = zeros(n_permutes,2);
max_clust_info  = zeros(n_permutes,1);

% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    fake_rt_mapping = rtsrank(randperm(EEG.trials));
    
    % compute t-map of null hypothesis
    fakecorrs = (fake_rt_mapping*fake_rt_mapping')\fake_rt_mapping*eegpowerrank;
    
    % reshape to 2D map for cluster-correction
    fakecorrs = reshape(fakecorrs,num_frex,nTimepoints);
    
    % save all permuted values
    permuted_rvals(permi,:,:) = fakecorrs;
    
    % save maximum pixel values
    max_pixel_rvals(permi,:) = [ min(fakecorrs(:)) max(fakecorrs(:)) ];
end

% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    % indices of permutations to include in thresholding at this iteration
    perms2use4distribution = true(1,n_permutes);
    perms2use4distribution(permi) = 0;
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    fakecorrsz = squeeze((permuted_rvals(permi,:,:)-mean(permuted_rvals(perms2use4distribution,:,:),1)) ./ std(permuted_rvals(perms2use4distribution,:,:),[],1) );
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(fakecorrsz);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
end

% now compute Z-map
zmap = (realcorrs-squeeze(mean(permuted_rvals,1)))./squeeze(std(permuted_rvals));

figure
subplot(221)
contourf(tftimes,frex,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-4 4],'xlim',[-500 1200])
title('Unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% apply uncorrected threshold
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
zmapthresh=logical(zmapthresh);
subplot(222)
contourf(tftimes,frex,zmap,40,'linecolor','none')
hold on
contour(tftimes,frex,zmapthresh,1,'linecolor','k')
axis square
set(gca,'clim',[-4 4],'xlim',[-500 1200])
title('Uncorrected thresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% apply pixel-level corrected threshold
lower_threshold = prctile(max_pixel_rvals(:,1),    mcc_voxel_pval*100/2);
upper_threshold = prctile(max_pixel_rvals(:,2),100-mcc_voxel_pval*100/2);

zmapthresh = zmap;
zmapthresh(realcorrs>lower_threshold & realcorrs<upper_threshold)=0;
subplot(223)
contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-4 4],'xlim',[-500 1200])
title('Pixel-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

subplot(224)
contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-4 4],'xlim',[-500 1200])
title('Cluster-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%% Figure 34.5

chan2use = 'o1';
time2use = dsearchn(EEG.times',[0 250]');
freq2use = dsearchn(frex',10);

eegfft = fft(reshape(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
eegconv = ifft(wavelets(freq2use,:).*eegfft);
eegconv = eegconv(1:n_convolution);
eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);

% reshape to time X trials
temp  = abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2;
o1power = zscore(mean(temp(time2use(1):time2use(2),:),1));

% define covariates (RT and trial number)
X = [ zscore(rts') o1power' ]';
eegpowerrank = tiedrank(eegpowerreshaped)';

figure
subplot(311)
imagesc(X)

subplot(212)
imagesc(eegpowerrank)

%% Figure 34.6

voxel_pval = 0.01;
mcc_cluster_pval = 0.05;

% note: try to use 1000 or more permutations for real data
n_permutes = 1000;

realbeta = (X*X')\X*eegpowerrank';
realbeta = reshape(realbeta,[2 num_frex nTimepoints]);

% initialize null hypothesis matrices
permuted_bvals = zeros(n_permutes,2,num_frex,nTimepoints);
max_clust_info = zeros(n_permutes,2);

% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    
    % randomly shuffle trial order
    fakeX = X(:,randperm(EEG.trials));
    
    % compute beta-map of null hypothesis
    fakebeta = (fakeX*fakeX')\fakeX*eegpowerrank';
    
    % reshape to 2D map for cluster-correction
    fakebeta = reshape(fakebeta,[2 num_frex nTimepoints ]);
    
    % save all permuted values
    permuted_bvals(permi,:,:,:) = fakebeta;
end

% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    for testi=1:2
        % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
        fakecorrsz = squeeze((permuted_bvals(permi,testi,:,:)-mean(permuted_bvals(:,testi,:,:),1)) ./ std(permuted_bvals(:,testi,:,:),[],1) );
        fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(fakecorrsz);
        max_clust_info(permi,testi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    end
end


figure
for testi=1:2
    
    % now compute Z-map
    zmap = (squeeze(realbeta(testi,:,:))-squeeze(mean(permuted_bvals(:,testi,:,:),1))) ./ squeeze(std(permuted_bvals(:,testi,:,:),[],1));
    
    subplot(2,3,1+(testi-1)*3)
    contourf(tftimes,frex,zmap,40,'linecolor','none')
    axis square
    set(gca,'clim',[-3 3],'xlim',[-500 1200])
    title('Unthresholded Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    
    % apply uncorrected threshold
    zmapthresh = zmap;
    zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
    subplot(2,3,2+(testi-1)*3)
    contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
    axis square
    set(gca,'clim',[-3 3],'xlim',[-500 1200])
    title('Uncorrected thresholded Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    
    % apply cluster-level corrected threshold
    zmapthresh = zmap;
    % uncorrected pixel-level threshold
    zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(zmapthresh);
    clust_info = cellfun(@numel,clustinfo.PixelIdxList);
    clust_threshold = prctile(max_clust_info(:,testi),100-mcc_cluster_pval*100);
    
    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    
    subplot(2,3,3+(testi-1)*3)
    contourf(tftimes,frex,zmapthresh,40,'linecolor','none')
    axis square
    set(gca,'clim',[-3 3],'xlim',[-500 1200])
    title('Cluster-corrected Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
end

%% Figure 34.7

a = rand(10000,1);
b = rand(10000,1);

figure
clear h

subplot(221)
[y,x] = hist(a,50);
h(1)=bar(x,y,'histc');
set(gca,'xlim',[-.05 1.05])

subplot(222)
[y,x] = hist(b,50);
h(2)=bar(x,y,'histc');
set(gca,'xlim',[-.05 1.05])


subplot(212)
[y,x] = hist(atanh(a-b),50);
h(3)=bar(x,y,'histc');
set(gca,'xlim',[-2 2])
title('ITPC differences')
xlabel('difference value'), ylabel('Count')

set(h,'linestyle','none','facecolor','k')

%% Figure 34.8

% The code to produce this figure is presented in chapter19.m, between the
% cells for figures 19.6 and 19.7. To generate figure 34.8, you will need to run
% the code for figures 19.2-6. 

%% Figure 34.9

voxel_pval   = 0.01;
cluster_pval = 0.05;

% note: try to use 1000 or more permutations for real data
n_permutes = 1000;

% compute actual t-test of difference
realitpc = squeeze(abs(mean(eegphase,3)));

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_frex);
permuted_vals    = zeros(n_permutes,num_frex,EEG.pnts);
max_clust_info   = zeros(n_permutes,1);
eegtemp          = zeros(size(realitpc));

for permi=1:n_permutes
    for triali=1:EEG.trials
        cutpoint = randsample(2:nTimepoints-2,1);
        eegtemp(:,:,triali) = eegphase(:,[cutpoint:end 1:cutpoint-1],triali);
    end
    permuted_vals(permi,:,:) = squeeze(abs(mean(eegtemp,3)));
    
    % note: the following lines produce fairly similar results as the loop above
    % cutpoint = randsample(2:nTimepoints-2,1);
    %permuted_vals(permi,:,:) = squeeze(abs(mean(eegphase(:,[cutpoint:end 1:cutpoint-1],:),3)));
end

zmap = (realitpc-squeeze(mean(permuted_vals))) ./ squeeze(std(permuted_vals));
threshmean = realitpc;
threshmean(abs(zmap)<norminv(1-voxel_pval))=0;

figure
subplot(221)
contourf(EEG.times,frex,realitpc,40,'linecolor','none')
axis square
set(gca,'clim',[0 .5],'xlim',[-200 1000])
title('power map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

subplot(222)
contourf(EEG.times,frex,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[-200 1000])
title('unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

subplot(223)
contourf(EEG.times,frex,threshmean,40,'linecolor','none')
axis square
set(gca,'clim',[0 .5],'xlim',[-200 1000])
title('Uncorrected power map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')



% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    fakecorrsz = squeeze((permuted_vals(permi,:,:)-mean(permuted_vals,1)) ./ std(permuted_vals,[],1) );
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(fakecorrsz);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    % using cellfun here eliminates the need for a slower loop over cells
end

% apply cluster-level corrected threshold
zmapthresh = realitpc;
% uncorrected pixel-level threshold
zmapthresh(abs(zmap)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

subplot(224)
contourf(EEG.times,frex,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[0 .5],'xlim',[-200 1000])
title('Cluster-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%% Figure 33.5

% The code for figures 33.5/6 are presented here and in the next cell. You
% will need first to run the code for figure 34.3 to run this code.

% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(mean(eegpower(4,400,real_condition_mapping==-1),3) - mean(eegpower(4,400,real_condition_mapping==1),3));
tdenom = sqrt( (std(eegpower(4,400,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) + (std(eegpower(4,400,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1) );
real_t = tnum./tdenom;

n_permutes = round(linspace(100,3000,200));
zvals   = zeros(size(n_permutes));

for grandpermi=1:length(n_permutes)
    
    % initialize null hypothesis matrices
    permuted_tvals  = zeros(n_permutes(grandpermi),1);
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:n_permutes(grandpermi)
        fake_condition_mapping = sign(randn(EEG.trials,1));
        % compute t-map of null hypothesis
        tnum   = squeeze(mean(eegpower(4,400,fake_condition_mapping==-1),3)-mean(eegpower(4,400,fake_condition_mapping==1),3));
        tdenom = sqrt( (std(eegpower(4,400,fake_condition_mapping==-1),0,3).^2)./sum(fake_condition_mapping==-1) + (std(eegpower(4,400,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) );
        % save all permuted values
        permuted_tvals(permi) = tnum./tdenom;
    end
    
    zvals(grandpermi) = (real_t-mean(permuted_tvals))/std(permuted_tvals);
    
    % display progress
    if mod(grandpermi,20)==0, disp([ 'Metapermutation #' num2str(grandpermi) ]); end
end
    
figure
subplot(211)
plot(n_permutes,zvals)
xlabel('Number of iterations'), ylabel('Z-value')

subplot(212)
hist(zvals,30)
xlabel('Z-value at different runs of permutation test'), ylabel('Count')

%% Figure 33.6

% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(mean(eegpower(:,:,real_condition_mapping==-1),3) - mean(eegpower(:,:,real_condition_mapping==1),3));
tdenom = sqrt( (std(eegpower(:,:,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) + (std(eegpower(:,:,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1) );
real_t = tnum./tdenom;


n_permutes = round(linspace(100,3000,200));
zvals = zeros(length(n_permutes),size(tnum,1),size(tnum,2));

for grandpermi=1:length(n_permutes)
    
    % initialize null hypothesis matrices
    permuted_tvals  = zeros(n_permutes(grandpermi),size(tnum,1),size(tnum,2));
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:n_permutes(grandpermi)
        fake_condition_mapping = sign(randn(EEG.trials,1));
        % compute t-map of null hypothesis
        tnum   = squeeze(mean(eegpower(:,:,fake_condition_mapping==-1),3)-mean(eegpower(:,:,fake_condition_mapping==1),3));
        tdenom = sqrt( (std(eegpower(:,:,fake_condition_mapping==-1),0,3).^2)./sum(fake_condition_mapping==-1) + (std(eegpower(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) );
        % save all permuted values
        permuted_tvals(permi,:,:) = tnum./tdenom;
    end
    
    zvals(grandpermi,:,:) = (real_t-squeeze(mean(permuted_tvals)))./squeeze(std(permuted_tvals));
    
    % display progress
    if mod(grandpermi,20)==0, disp([ 'Metapermutation # ' num2str(grandpermi) ]); end
end



figure
subplot(211)
plot(n_permutes,squeeze(zvals(:,4,400)))
subplot(212)
hist(squeeze(zvals(:,4,400)),30)


zvalsall = reshape(zvals,length(n_permutes),size(tnum,1)*size(tnum,2));
figure
plot(mean(zvalsall,1),std(zvalsall,[],1),'.')
set(gca,'xlim',[-3.5 3.5],'ylim',[0 .12])
xlabel('Average Z-statistic')
ylabel('Standard deviation of Z-statistics')

figure
clear h

[~,z0]=min(abs(mean(zvalsall,1)));
[x0,y0]=ind2sub(size(tnum),z0);
[yy,xx]=hist(squeeze(zvals(:,x0,y0)),30);
h(1) = bar(xx,yy,'histc');
hold on
plot([0 0],get(gca,'ylim'),'k:')

[~,z2]=min(abs(mean(zvalsall,1)-2));
[x2,y2]=ind2sub(size(tnum),z2);
[yy,xx]=hist(squeeze(zvals(:,x2,y2)),30);
h(2) = bar(xx,yy,'histc');
plot([2 2],get(gca,'ylim'),'k:')

[~,z3]=min(abs(mean(zvalsall,1)--3));
[x3,y3]=ind2sub(size(tnum),z3);
[yy,xx]=hist(squeeze(zvals(:,x3,y3)),30);
hold on
h(3) = bar(xx,yy,'histc');
plot([-3 -3],get(gca,'ylim'),'k:')

set(h,'linestyle','none')
set(gca,'xlim',[-3.5 3.5])
xlabel('Z-value')
ylabel('Count of possible z-values')

%% end.
