%% Analyzing Neural Time Series Data
% Matlab code for Chapter 31
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 31.2

% load data
load sampleEEGdata

% specify some time-frequency parameters
center_freq   = 10; % Hz
time2analyze  = 200; % in ms

% wavelet and FFT parameters
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;

% initialize connectivity output matrix
connectivitymat = zeros(EEG.nbchan,EEG.nbchan);

% time in indices
[junk,tidx] = min(abs(EEG.times-time2analyze));

% create wavelet and take FFT
s = 5/(2*pi*center_freq);
wavelet_fft = fft( exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);

% compute analytic signal for all channels
analyticsignals = zeros(EEG.nbchan,EEG.pnts,EEG.trials);
for chani=1:EEG.nbchan
    
    % FFT of data
    data_fft = fft(reshape(EEG.data(chani,:,:),1,n_data),n_convolution);
    
    % convolution
    convolution_result = ifft(wavelet_fft.*data_fft,n_convolution);
    convolution_result = convolution_result(half_wavelet+1:end-half_wavelet);
    
    analyticsignals(chani,:,:) = reshape(convolution_result,EEG.pnts,EEG.trials);
end

% now compute all-to-all connectivity
for chani=1:EEG.nbchan
    for chanj=chani:EEG.nbchan % note that you don't need to start at 1
        xsd = squeeze(analyticsignals(chani,tidx,:) .* conj(analyticsignals(chanj,tidx,:)));
        
        % connectivity matrix (phase-lag index on upper triangle; ISPC on lower triangle)
        connectivitymat(chani,chanj) = abs(mean(sign(imag(xsd))));
        connectivitymat(chanj,chani) = abs(mean(exp(1i*angle(xsd))));
    end
end

figure
imagesc(connectivitymat)
set(gca,'clim',[0 .7],'xtick',1:8:EEG.nbchan,'xticklabel',{EEG.chanlocs(1:8:end).labels},'ytick',1:8:EEG.nbchan,'yticklabel',{EEG.chanlocs(1:8:end).labels});
axis square
colorbar

%% Figure 31.4

figure

subplot(221)
% get upper part of matrix
temp = nonzeros(triu(connectivitymat));
temp(temp==1)=[]; % clear 1's from the diagonal

% threshold is one std above median connectivity value
pli_thresh = std(temp)+median(temp);

% plot histogram and vertical line at threshold
[y,x] = hist(temp,30);
h = bar(x,y,'histc');
hold on
plot([pli_thresh pli_thresh],get(gca,'ylim'),'m','linew',2)

% make nice
set(h,'linestyle','none')
set(gca,'xtick',0:.2:1,'xlim',[0 1])
xlabel('Phase-lag index'), ylabel('Count')



subplot(222)
% get upper part of matrix
temp = nonzeros(tril(connectivitymat));
temp(temp==1)=[]; % clear 1's on the diagonal

% find 1 std above median connectivity value
ispc_thresh = std(temp)+median(temp);

% plot histogram and vertical line at threshold
[y,x] = hist(temp,30);
h=bar(x,y,'histc');
hold on
plot([ispc_thresh ispc_thresh],get(gca,'ylim'),'m','linew',2)

% make nice
set(h,'linestyle','none')
set(gca,'xtick',0:.2:1,'xlim',[0 1])
xlabel('ISPC'), ylabel('Count')


subplot(223)

% make symmetric phase-lag index connectivity matrix
pli_mat = connectivitymat;
pli_mat(logical(tril(pli_mat))) = 0; % eliminate lower triangle
pli_mat = pli_mat + triu(pli_mat)';  % mirror lower triangle to upper triangle
pli_mat(pli_mat<pli_thresh)=0;
imagesc(pli_mat)
set(gca,'clim',[0 .7],'xtick',1:8:EEG.nbchan,'xticklabel',{EEG.chanlocs(1:8:end).labels},'ytick',1:8:EEG.nbchan,'yticklabel',{EEG.chanlocs(1:8:end).labels});
axis square


subplot(224)

% make symmetric phase-lag index connectivity matrix
ispc_mat = connectivitymat;
ispc_mat(logical(triu(ispc_mat))) = 0; % eliminate lower triangle
ispc_mat = ispc_mat + tril(ispc_mat)';  % mirror lower triangle to upper triangle
ispc_mat(ispc_mat<ispc_thresh)=0;
imagesc(logical(ispc_mat)) % logical converts to 0's and 1's, thus binarizing connectivity matrix
set(gca,'clim',[0 .7],'xtick',1:8:EEG.nbchan,'xticklabel',{EEG.chanlocs(1:8:end).labels},'ytick',1:8:EEG.nbchan,'yticklabel',{EEG.chanlocs(1:8:end).labels});
axis square

%% Figure 31.6

figure
subplot(121)
% note: logical() below converts the matrix to 0's and 1's
topoplot(sum(logical(pli_mat)),EEG.chanlocs,'plotrad',.53,'maplimits',[0 25]);
title('Connectivity degree, phase-lag index')

subplot(122)
topoplot(sum(logical(ispc_mat)),EEG.chanlocs,'plotrad',.53,'maplimits',[0 25]);
title('Connectivity degree, ISPC')

%% Figure 31.7 (analyses are here; next cell plots results)

frex = logspace(log10(3),log10(40),25);
times2save = -300:50:1200;
thresh = zeros(size(frex)); % added threshold vector

% wavelet and FFT parameters
time          = -1:1/EEG.srate:1;
half_wavelet  = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv2       = pow2(nextpow2(n_convolution));

% create wavelet (and take FFT)
wavelets_fft = zeros(length(frex),n_conv2);
s = logspace(log10(4),log10(10),length(frex))./(2*pi.*frex);
for fi=1:length(frex)
    wavelets_fft(fi,:) = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) ,n_conv2);
end

% find time indices
times2saveidx = dsearchn(EEG.times',times2save');

% initialize matrices
alldata    = zeros(EEG.nbchan,length(frex),length(times2save),EEG.trials);
tf_all2all = zeros(EEG.nbchan,EEG.nbchan,length(frex),length(times2save));
tf_degree  = zeros(EEG.nbchan,length(frex),length(times2save));

% first, run convolution for all electrodes and save results
for chani=1:EEG.nbchan
    
    % FFT of activity at this electrode (note that 
    % this is done outside the frequency loop)
    eeg_fft = fft(reshape(EEG.data(chani,:,:),1,[]),n_conv2);
    
    % loop over frequencies
    for fi=1:length(frex)
        
        % analytic signal from target
        conv_res = ifft(wavelets_fft(fi,:).*eeg_fft,n_conv2);
        conv_res = conv_res(1:n_convolution);
        asig     = reshape(conv_res(half_wavelet+1:end-half_wavelet),EEG.pnts,EEG.trials);
        
        % store the required time points
        alldata(chani,fi,:,:) = asig(times2saveidx,:);
    end % end frequency loop
end % end channel loop


% now that we have all the data, compute all-to-all connectivity
for chani=1:EEG.nbchan
    for chanj=chani+1:EEG.nbchan
        
        % compute connectivity
        xsd = squeeze(alldata(chani,:,:,:).*conj(alldata(chanj,:,:,:)));
        
        % connectivity matrix (phase-lag index or ISPC; comment one or the other line)
        % tf_all2all(chani,chanj,:,:) = abs(mean(sign(imag(xsd)),3)); % pli
        tf_all2all(chani,chanj,:,:) = abs(mean(exp(1i*angle(xsd)),3)); % ispc
    end
end

% now that we have a one-to-all connectivity, threshold the 
% connectivity matrix (separate threshold for each frequency
for fi=1:length(frex)
    
    % Note about the threshold: The original code had a bug where the
    % threshold for all analyses and all frequencies was based only on the
    % final frequency. Thanks to Lars Benschop for finding that bug!
    
    tempsynch = nonzeros(tf_all2all(:,:,fi,:));
    thresh(fi) = median(tempsynch) + std(tempsynch);
    
    % isolate, threshold, binarize
    for ti=1:size(tf_all2all,4)
        temp = squeeze(tf_all2all(:,:,fi,ti));
        temp = temp + triu(temp)';      % make symmetric matrix
        temp = temp>thresh(fi);         % threshold and binarize
        tf_degree(:,fi,ti) = sum(temp); % compute degree (sum of suprathreshold connections)
    end
end

%% Figure 31.7 (plotting)

% show topographical maps
freqs2plot = [ 5 9 ]; % Hz
times2plot = [ 100 200 300 ]; % ms
clim = [0 20];

figure
for fi=1:length(freqs2plot)
    for ti=1:length(times2plot)
        
        subplot(length(freqs2plot),length(times2plot),ti+(fi-1)*length(times2plot))
        topoplot(squeeze(tf_degree(:,dsearchn(frex',freqs2plot(fi)),dsearchn(times2save',times2plot(ti)))),EEG.chanlocs,'plotrad',.53,'maplimits',clim,'electrodes','off');
        
        title([ num2str(times2plot(ti)) ' ms, ' num2str(freqs2plot(fi)) ' Hz' ])
        
    end
end

%% Figure 31.8

electrode2plot = 'fcz';'oz';
baselineperiod = [ -300 -100 ];
clim = [-10 10];


% convert baseline period time to idx
baseidx = dsearchn(times2save',baselineperiod');
% subtract baseline
tf_degree_base = tf_degree - repmat(mean(tf_degree(:,:,baseidx(1):baseidx(2)),3),[1 1 size(tf_degree,3)]);


figure
contourf(times2save,frex,squeeze(tf_degree_base(strcmpi(electrode2plot,{EEG.chanlocs.labels}),:,:)),20,'linecolor','none')
set(gca,'clim',clim,'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),6)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'TF connectivity degree at electrode ' electrode2plot ])

%% Figure 31.10

freqs2plot = [ 5 9 ]; % Hz
times2plot = [ 100 200 300 ]; % ms
clim = [.4 .8];

% open two figures ('h' for handle)
figh1 = figure;
figh2 = figure;

for fi=1:length(freqs2plot)
    for ti=1:length(times2plot)
        
        % frequency index
        fidx = dsearchn(frex',freqs2plot(fi));
        
        % extract thresholded connectivity matrix
        connmat = squeeze(tf_all2all(:,:,fidx,dsearchn(times2save',times2plot(ti))));
        connmat = connmat + triu(connmat)';  % make symmetric matrix
        connmat = connmat>thresh(fidx);      % threshold and binarize
        
        
        % initialize
        clustcoef = zeros(size(EEG.nbchan));
        
        for chani=1:EEG.nbchan
            
            % find neighbors (suprathreshold connections)
            neighbors = find(connmat(chani,:));
            n = length(neighbors);
            
            % cluster coefficient not computed for islands
            if n>1
                % "local" network of neighbors
                localnetwork = connmat(neighbors,neighbors);
                % localnetwork is symmetric; remove redundant values by replacing with NaN
                localnetwork = localnetwork + tril(nan(n));
                
                % compute cluster coefficient (neighbor connectivity scaled)
                clustcoef(chani) = 2*nansum(localnetwork(:)) / ((n-1)*n);
            end
        end
        
        % topoplots
        figure(figh1)
        subplot(length(freqs2plot),length(times2plot),ti+(fi-1)*length(times2plot))
        topoplot(clustcoef,EEG.chanlocs,'plotrad',.53,'maplimits',clim,'electrodes','off');
        title([ 'Cluster coefficient: ' num2str(times2plot(ti)) ' ms, ' num2str(freqs2plot(fi)) ' Hz' ])
        
        % relationship between degree and cluster coefficient
        figure(figh2)
        subplot(length(freqs2plot),length(times2plot),ti+(fi-1)*length(times2plot))
        plot(squeeze(tf_degree(:,dsearchn(frex',freqs2plot(fi)),dsearchn(times2save',times2plot(ti)))),clustcoef,'.')
        [r,p] = corr(squeeze(tf_degree(:,dsearchn(frex',freqs2plot(fi)),dsearchn(times2save',times2plot(ti)))),clustcoef','type','s');
        set(gca,'ylim',[-.05 1.05],'xlim',[0 30])
        xlabel('degree'), ylabel('cluster coefficient')
        title([ 'Correlation: r=' num2str(r) ', p=' num2str(p) ])
        
    end
end

%% Figure 31.11

thresholds = linspace(0,3,20);

clustercoefficients = zeros(EEG.nbchan,length(thresholds));

for ti=1:length(thresholds)
    
    % find frequency index
    fidx = dsearchn(frex',freqs2plot(1));
    
    % extract thresholded connectivity matrix
    connmat = squeeze(tf_all2all(:,:,fidx,dsearchn(times2save',times2plot(3))));
    tmpsynch = nonzeros(tf_all2all(:,:,fidx,:));
    tthresh = median(tmpsynch)+thresholds(ti)*std(tmpsynch);
    connmat = connmat + triu(connmat)';  % make symmetric matrix
    connmat = connmat>tthresh;           % threshold and binarize
    
    
    for chani=1:EEG.nbchan
        
        % find neighbors (suprathreshold connections)
        neighbors = find(connmat(chani,:));
        n = length(neighbors);
        
        % cluster coefficient not computed for islands
        if n>1
            % "local" network of neighbors
            localnetwork = connmat(neighbors,neighbors);
            % localnetwork is symmetric; remove redundant values by replacing with NaN
            localnetwork = localnetwork + tril(nan(n));
            
            % compute cluster coefficient (neighbor connectivity scaled)
            clustercoefficients(chani,ti) = 2*nansum(localnetwork(:)) / ((n-1)*n);
        end
    end
end

figure
plot(thresholds,clustercoefficients)
hold on
plot(thresholds,mean(clustercoefficients,1),'k','linew',3)
xlabel('Threshold (number of standard deviations above median)')
ylabel('Clustering coefficient')
set(gca,'ylim',[-.025 1.025])

%% Figures 31.13/14 (this cell takes a long time to run...)

use_real_network = false; % if false, simulate a network (true/false for Figures 13/14)

freq2use = 10;
time2use = 100;

% frequency index
fidx = dsearchn(frex',freq2use);

connmat = squeeze(tf_all2all(:,:,fidx,dsearchn(times2save',time2use)));
connmat = connmat + triu(connmat)';  % make symmetric matrix
connmat = connmat>thresh(fidx);      % threshold and binarize

if ~use_real_network
    n = 1000; % number of nodes
    k = 10;  % neighbor connectivity
    connmat = zeros(n);
    for i=1:n
        % set neighbors to 1 (special cases for start and end of network)
        connmat(i,max(1,i-k/2):min(n,i+k/2)) = 1;
    end
end

% probabilities of rewiring
probs = logspace(log10(0.0001),log10(1),20);

cp = zeros(size(probs));
lp = zeros(size(probs));

% loop over a few networks to make plots look a bit nicer
for neti=1:10
    for probi=1:length(probs)
        
        % rewire
        % find which edges to rewire
        real_edges = find(tril(connmat));
        real_edges = real_edges(randperm(length(real_edges)));
        edges2rewire = real_edges(1:round(probs(probi)*length(real_edges)));
        
        % rewired connectivity matrix
        connmat_rewired = connmat;
        
        % loop through edges and change target
        for ei=1:length(edges2rewire)
            
            % find XY coordinates
            [x,y] = ind2sub(size(connmat),edges2rewire(ei));
            
            % find possible edges to change (cannot already be an edge)
            edges2change = find(~connmat_rewired(x,:));
            
            % rewire
            y2rewire = randsample(edges2change,1);
            connmat_rewired(x,y2rewire) = 1;
            connmat_rewired(y2rewire,x) = 1;
            
            % set original to zero
            connmat_rewired(x,y) = 0;
            connmat_rewired(y,x) = 0;
        end
        
        % mirror matrix
        connmat_rewired = logical(tril(connmat_rewired) + tril(connmat_rewired)');
        
        
        clustcoef_rewired  = zeros(1,EEG.nbchan);
        pathlength_rewired = zeros(1,EEG.nbchan);
        for chani=1:EEG.nbchan
            
            % cluster coefficient
            neighbors = find(connmat_rewired(chani,:));
            n = length(neighbors);
            
            if n>1
                % "local" network of neighbors
                localnetwork = connmat_rewired(neighbors,neighbors);
                % localnetwork is symmetric; remove redundant values by replacing with NaN
                localnetwork = localnetwork + tril(nan(n));
                % compute cluster coefficient (neighbor connectivity scaled)
                clustcoef_rewired(chani) = 2*nansum(localnetwork(:)) / ((n-1)*n);
            end
            
        end
        % average clustering coefficient over channels
        cp(probi) = cp(probi) + nanmean(clustcoef_rewired);
        
        % average path length (remove zeros and Inf's)
        temppathlengths = nonzeros(pathlength(double(connmat_rewired)));
        lp(probi) = lp(probi) + mean(temppathlengths(isfinite(temppathlengths)));
        
        
        % save example networks from select probabilities
        if probi==1
            network1 = connmat_rewired;
        elseif probi==10
            network10 = connmat_rewired;
        elseif probi==length(probs)
            networkend = connmat_rewired;
        end
        
    end
end

cp = cp./neti;
lp = lp./neti;

figure
plot(probs,cp./cp(1),'-o')
hold on
plot(probs,lp./lp(1),'r-o')
xlabel('Probability of rewiring')
legend({'C_r/C';'L_r/L'})
set(gca,'xlim',[probs(1)/1.5 1.2],'ylim',[0 1.02],'xscale','lo')
plot([probs(1) probs(1)],[0 1],'k:')
plot([probs(10) probs(10)],[0 1],'k:')
plot([probs(end) probs(end)],[0 1],'k:')


if use_real_network
    xylims = [1 EEG.nbchan];
else xylims = [1 100];
end

figure
subplot(131)
imagesc(network1)
set(gca,'xlim',xylims,'ylim',xylims), axis square, colormap(1-gray)

subplot(132)
imagesc(network10)
set(gca,'xlim',xylims,'ylim',xylims), axis square

subplot(133)
imagesc(networkend)
set(gca,'xlim',xylims,'ylim',xylims), axis square

%% Figures 31.15/6

freq2use = 6;
time2use = 300;
n_permutations = 1000;

thresholds=linspace(0,2,20);
swn_z = zeros(size(thresholds));

for threshi=1:length(thresholds)
    
    % frequency index
    fidx = dsearchn(frex',freq2use);
    
    connmat = squeeze(tf_all2all(:,:,fidx,dsearchn(times2save',time2use)));
    tmpsynch = nonzeros(tf_all2all(:,:,fidx,:));
    tthresh = median(tmpsynch)+thresholds(threshi)*std(tmpsynch);
    connmat = connmat + triu(connmat)';  % make symmetric matrix
    connmat = connmat>tthresh;           % threshold and binarize
    nconnections = sum(connmat(:))/2;    % number of connections (for random network)
    
    % find locations of lower triangle
    matrix_locs = find(tril(connmat+1-eye(length(connmat))));
    
    swn_permutations  = zeros(1,n_permutations);
    swn_permutations2 = zeros(1,n_permutations);
    
    
    
    % compute real cluster coefficient and path lengths
    clustcoef_real = zeros(1,EEG.nbchan);
    for chani=1:EEG.nbchan
        neighbors = find(connmat(chani,:));
        n = length(neighbors);
        if n>1
            localnetwork = connmat(neighbors,neighbors);
            localnetwork = localnetwork + tril(nan(n));
            clustcoef_real(chani) = 2*nansum(localnetwork(:)) / ((n-1)*n);
        end
    end
    % average clustering coefficient over channels
    clustcoef_real = nanmean(clustcoef_real);
    
    % average path length (remove zeros and Inf's)
    temppathlengths = nonzeros(pathlength(double(connmat)));
    pathlengths_real = mean(temppathlengths(isfinite(temppathlengths)));
    
    
    
    
    % first create 1000 random graphs and compute their CC and PL
    clustercoefficients_random = zeros(1,n_permutations);
    pathlengths_random = zeros(1,n_permutations);
    
    for permi=1:n_permutations
        
        % generate random network...
        connmat_random = zeros(size(connmat));
        edges2fill     = randsample(matrix_locs,nconnections);
        connmat_random(edges2fill) = 1;
        % mirror matrix
        connmat_random = logical(tril(connmat_random) + tril(connmat_random)');
        
        clustcoef_random  = zeros(1,EEG.nbchan);
        for chani=1:EEG.nbchan
            neighbors = find(connmat_random(chani,:));
            n = length(neighbors);
            if n>1
                localnetwork = connmat_random(neighbors,neighbors);
                localnetwork = localnetwork + tril(nan(n));
                clustcoef_random(chani) = 2*nansum(localnetwork(:)) / ((n-1)*n);
            end
        end
        clustercoefficients_random(permi) = nanmean(clustcoef_random);
        
        % average path length (remove zeros and Inf's)
        temppathlengths = nonzeros(pathlength(double(connmat_random)));
        pathlengths_random(permi) = mean(temppathlengths(isfinite(temppathlengths)));
    end
    
    % now compute permuted small-world-network-ness
    for permi=1:n_permutations
        whichnetworks2use = randsample(1:n_permutations,2);
        swn_permutations(permi)  = ( clustercoefficients_random(whichnetworks2use(1))/clustercoefficients_random(whichnetworks2use(2)) ) / ( pathlengths_random(whichnetworks2use(1))/pathlengths_random(whichnetworks2use(2)) );
    end
    
    
    % true swn and its z-value
    swn_real       = ( clustcoef_real/mean(clustercoefficients_random) ) / ( pathlengths_real/mean(pathlengths_random) );
    swn_z(threshi) = (swn_real-mean(swn_permutations))/std(swn_permutations);
    
    if threshi==length(thresholds)/2
        figure
        [y,x] = hist(swn_permutations,50);
        h=bar(x,y,'histc');
        set(h,'linestyle','none')
        hold on
        plot(repmat(swn_real,1,2),get(gca,'ylim')/2,'m','linew',3)
        title([ 'Small-world-networkness: Z=' num2str(swn_z(threshi)) ', p=' num2str(1-normcdf(abs(swn_z(threshi)))) ])
    end
end

figure
plot(thresholds,swn_z,'-o','markerface','w')
xlabel('Threshold (Standard deviations above median)')
ylabel('SWN_z')
set(gca,'xlim',[-.05 2.05])

%% end.
