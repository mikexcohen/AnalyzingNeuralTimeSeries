%% Analyzing Neural Time Series Data
% Matlab code for Chapter 29
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 29.2

% create two signals
time    = 0:.0001:1;
signal1 = sin(2*pi*10*time);
signal2 = rand(size(signal1))*2-1; % uniform random numbers in the same scale as the sine wave

% plot signals
figure
subplot(221)
plot(time,signal1)
set(gca,'xlim',[time(1) time(end)])
subplot(222)
plot(time,signal2,'r')
set(gca,'xlim',[time(1) time(end)])

% bin data (via Matlab function hist)
nbins = 50;
[hdat1,x1] = hist(signal1,nbins);
[hdat2,x2] = hist(signal2,nbins);

% convert histograms to probability values
hdat1 = hdat1./sum(hdat1);
hdat2 = hdat2./sum(hdat2);

% plot histograms
subplot(223)
plot(x1,hdat1)
hold on
plot(x2,hdat2,'r')
legend({'Sine wave';'Random data'})
xlabel('Value bins')
ylabel('Probability')

subplot(224)
plot(sort(signal1))
hold on
plot(sort(signal2),'r')
set(gca,'xlim',[0 length(signal1)])

% The following code using the 'eval' command may seem a bit needlessly
% complex, but it introduces you to this useful Matlab function. The
% main advantage here is to call different variables inside a loop,
% which would otherwise not be possible because it is the variable
% names that are different, not indices into variables. For example,
% notice how the variables hdat1 and hdat2 are called.
for i=1:2
    eval([ 'entro(' num2str(i) ') = -sum(hdat' num2str(i) '.*log2(hdat' num2str(i) '+eps));' ]);
end

% the following code will do this same thing:
% entro(1) = -sum(hdat1.*log2(hdat1+eps));
% entro(2) = -sum(hdat2.*log2(hdat2+eps));

disp([ 'Entropies of sine wave and random noise are ' num2str(entro(1)) ' and ' num2str(entro(2)) '.' ]);

%% Figure 29.3

% range of bin numbers
nbins = 10:2000;

entropyByBinSize=zeros(size(nbins));

for nbini=1:length(nbins)
    
    % bin data, transform to probability, and eliminate zeros
    hdat = hist(signal1,nbins(nbini));
    hdat = hdat./sum(hdat);
    
    % compute entropy
    entropyByBinSize(nbini) = -sum(hdat.*log2(hdat+eps));
end

figure
plot(nbins,entropyByBinSize)
xlabel('Number of bins'), ylabel('Entropy')

%% Figure 29.4

% optimal number of bins for histogram based on a few different guidelines
n = length(signal1);
maxmin_range = max(signal1)-min(signal1);

% note: the function iqr (inter-quartile range) is in the stats toolbox. 
% If you don't have this toolbox, you can write your own similar function
% by sorting the values, finding the values that are 25% and 75% of the
% sorted distribution, and then subtracting the 25% number from the 75% number.
fd_bins      = ceil(maxmin_range/(2.0*iqr(signal1)*n^(-1/3))); % Freedman-Diaconis 
scott_bins   = ceil(maxmin_range/(3.5*std(signal1)*n^(-1/3))); % Scott
sturges_bins = ceil(1+log2(n)); % Sturges

figure
subplot(211)
% plot up to 50 bins
[junk,maxNbins] = min(abs(nbins-50)); % index of nbins that most closely matches 50
plot(nbins(1:maxNbins),entropyByBinSize(1:maxNbins))

hold on
plot([fd_bins fd_bins],get(gca,'ylim'),'m','linew',2)
plot([scott_bins scott_bins],get(gca,'ylim'),'k','linew',2)
plot([sturges_bins sturges_bins],get(gca,'ylim'),'r','linew',2)

legend({'entropy';'Freedman-Diaconis';'Scott';'Sturges'})
xlabel('Number of bins'), ylabel('Entropy')

subplot(223)
[y,x]=hist(signal1,fd_bins);
h=bar(x,y,'histc');
set(h,'linestyle','none');
set(gca,'xlim',[min(signal1) max(signal1)*1.1])
xlabel('Value'), ylabel('Count')
title('Optimal number of bins (FD rule)')

subplot(224)
[y,x]=hist(signal1,2000);
h=bar(x,y,'histc');
set(gca,'xlim',[min(signal1) max(signal1)]*1.05)
xlabel('Value'), ylabel('Count')
title('Too many bins (2000)')

%% Figure 29.5

load sampleEEGdata

% entropy over all sensors
time4entropy = [  100  400 ]; % in ms
base4entropy = [ -400 -100 ]; % in ms
topo_entropy = zeros(size(EEG.chanlocs));

timeidx=zeros(size(time4entropy)); baseidx=zeros(size(base4entropy)); 
for i=1:2
    [junk,timeidx(i)] = min(abs(EEG.times-time4entropy(i)));
    [junk,baseidx(i)] = min(abs(EEG.times-base4entropy(i)));
end

for chani=1:EEG.nbchan
    
    % entropy during task
    tempdat = EEG.data(chani,timeidx(1):timeidx(2),:);
    hdat = hist(tempdat(:),25);
    hdat = hdat./sum(hdat);
    task_entropy =  -sum(hdat.*log2(hdat+eps));
    
    % entropy during pre-stim baseline
    tempdat = EEG.data(chani,baseidx(1):baseidx(2),:);
    hdat = hist(tempdat(:),25);
    hdat = hdat./sum(hdat);
    base_entropy =  -sum(hdat.*log2(hdat+eps));
    
    % compute entropy
    topo_entropy(chani) = task_entropy - base_entropy;
end

figure
% Note: topoplot is a function in the eeglab toolbox
topoplot(topo_entropy,EEG.chanlocs,'maplimits',[-.5 .5],'plotrad',.53,'electrodes','off','numcontour',0);



% entropy over time in one electrode (fcz or po8)
sensor4entropy = 'fcz';
times2save     = -300:50:1200;
timewindow     = 400; % ms

timewindowidx = round(timewindow/(1000/EEG.srate)/2);
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end

electrodeidx = find(strcmpi(sensor4entropy,{EEG.chanlocs.labels}));

timeEntropy = zeros(1,length(times2save));

for timei = 1:length(times2save)
    
    tempdata = EEG.data(electrodeidx,times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,1:30);
    hdat = hist(tempdata(:),25);
    hdat = hdat./sum(hdat);
    
    % compute entropy
    timeEntropy(1,timei) = -sum(hdat.*log2(hdat+eps));
    
    
    tempdata = EEG.data(electrodeidx,times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,end-30:end);
    hdat = hist(tempdata(:),25);
    hdat = hdat./sum(hdat);
    
    % compute entropy
    timeEntropy(2,timei) = -sum(hdat.*log2(hdat+eps));
end

figure
plot(times2save,timeEntropy)
xlabel('Time (ms)'), ylabel('Entropy (bits)')
title([ 'Entropy over time from electrode ' sensor4entropy ])
legend({'First 30 trials';'last 30 trials'})

%% Figure 29.6

% Note about this figure: The panels use the same code but with different
% input signals. Comment out some of the lines below to recreate each panel.

% right panel: random noise
signal1 = rand(size(signal1))*2-1;
signal2 = rand(size(signal1))*2-1;

% center panel: one pure sine wave and one sine wave plus random noise
signal1 = sin(2*pi*10*time);
signal2 = signal1 + randn(size(signal1))/2;

% left panel: one pure sine wave and its inverse
signal1 = sin(2*pi*10*time);
signal2 = -signal1;


% determine the optimal number of bins for each variable
n            = length(signal1);
maxmin_range = max(signal1)-min(signal1);
fd_bins1     = ceil(maxmin_range/(2.0*iqr(signal1)*n^(-1/3))); % Freedman-Diaconis 

n            = length(signal2);
maxmin_range = max(signal2)-min(signal2);
fd_bins2     = ceil(maxmin_range/(2.0*iqr(signal2)*n^(-1/3))); % Freedman-Diaconis 

% and use the average...
fd_bins = ceil((fd_bins1+fd_bins2)/2);


% bin data (using histc this time)
edges = linspace(min(signal1),max(signal1),fd_bins+1);
[nPerBin1,bins1] = histc(signal1,edges);

edges = linspace(min(signal2),max(signal2),fd_bins+1);
[nPerBin2,bins2] = histc(signal2,edges);

% compute joint frequency table
jointprobs = zeros(fd_bins);
for i1=1:fd_bins
    for i2=1:fd_bins
        jointprobs(i1,i2) = sum(bins1==i1 & bins2==i2);
    end
end
jointprobs=jointprobs./sum(jointprobs(:));

figure
subplot(211)
plot(time,signal1)
subplot(212)
plot(time,signal2)

figure
imagesc(jointprobs)
colormap gray
set(gca,'clim',[0 .01],'ydir','normal')
colorbar
xlabel('Signal 2 bin'), ylabel('Signal 1 bin')

%% Figure 29.7

figure
subplot(221)
x = 0:.001:1;
y = x;
plot(x,y,'.')
title([ 'MI=' num2str(mutualinformationx(x,y)) ', r_s=' num2str(corr(x(:),y(:),'type','s')) ])
axis square

subplot(222)
x = 0:.001:1;
y = -x.^3;
plot(x,y,'.')
title([ 'MI=' num2str(mutualinformationx(x,y)) ', r_s=' num2str(corr(x(:),y(:),'type','s')) ])
axis square

subplot(223)
x=cos(0:.01:2*pi);
y=sin(0:.01:2*pi);
plot(x,y,'.')
title([ 'MI=' num2str(mutualinformationx(x,y)) ', r_s=' num2str(corr(x(:),y(:),'type','s')) ])
axis square

subplot(224)
x=[cos(0:.01:2*pi) cos(0:.01:2*pi)+1];
y=[sin(0:.01:2*pi) sin(0:.01:2*pi)-1];
plot(x,y,'.')
title([ 'MI=' num2str(mutualinformationx(x,y)) ', r_s=' num2str(corr(x(:),y(:),'type','s')) ])
axis square

%% Figure 29.8

% Theoretical size of errors as a function of histogram bins and N
% the Matlab 'inline' function is a useful tool to substitute small
% functions (1-2 lines of code) that only need to be run locally.
% However, as mentioned in chapter 26, this function will be removed in
% future versions of Matlab. 
entropy_error = inline('(b-1)./(2.*n.*log(2))');
mutinfo_error = inline('(b-1).^2./(2.*n.*log(2))');

n = 20:300;
nfixed = 15;

figure
subplot(211)
plot(n,entropy_error(nfixed,n))
hold on
plot(n,entropy_error(ceil(1+log2(n)),n),'r')
legend({[ num2str(nfixed) ' bins' ];'Sturges'' rule'})
xlabel('Number of data points'), ylabel('Entropy error (bits)')
set(gca,'xlim',[n(1) n(end)])

subplot(212)
plot(n,mutinfo_error(nfixed,n))
hold on
plot(n,mutinfo_error(ceil(1+log2(n)),n),'r')
legend({[ num2str(nfixed) ' bins' ];'Sturges'' rule'})
xlabel('Number of data points'), ylabel('Mutual information error')
set(gca,'xlim',[n(1) n(end)])

%% Figure 29.9

x = [cos(0:.01:2*pi) cos(0:.01:2*pi)+1];
y = [sin(0:.01:2*pi) sin(0:.01:2*pi)-1];

figure
subplot(221)
plot(x,y,'.')
axis([-3 3 -3 3])

subplot(212)

noiselevels = 0:.01:1;
mi = zeros(size(noiselevels));
for ni=1:length(noiselevels)
    mi(ni) = mutualinformationx(x,y+randn(size(y))*noiselevels(ni),20);
end

plot(noiselevels,mi)
xlabel('Noise level'), ylabel('Mutual information')

subplot(222)
plot(x,y+randn(size(y))*noiselevels(round(ni/2)),'.');
axis([-3 3 -3 3])

%% Figure 29.9d (takes a while to run)

nrange = 300:205:3000;
noiselevels = 0:.01:1;

mi = zeros(length(nrange),length(noiselevels));
b  = zeros(1,length(nrange)); % number of histogram bins

for ni=1:length(nrange)
    
    % define time
    t = linspace(0,2*pi,nrange(ni));
    % define signals
    x = [cos(t) cos(t)+1];
    y = [sin(t) sin(t)-1];
    
    for noi=1:length(noiselevels)
        if noi==1 % keep number of bins constant across noise levels within each number of points
            [mi(ni,noi),~,b(ni)] = mutualinformationx(x,y+randn(size(y))*noiselevels(noi));
        else
            mi(ni,noi) = mutualinformationx(x,y+randn(size(y))*noiselevels(noi),b(ni));
        end
    end
end

% convert to % change from best-case scenario (no noise, large N)
mip=100.*(mi-mi(end,1))./mi(end,1);

figure
contourf(noiselevels,nrange,mip,40,'linecolor','none')
set(gca,'clim',[-100 0])
colorbar
xlabel('Noise level'), ylabel('N (data length)')
title('Percent decrease in MI due to noise')

%% Figure 29.10

electrodes4mi = {'fz';'o1'};
timewindow = 400; % in ms
times2save = -400:100:1200;


% convert ms to indices
timewindowidx = round(timewindow/(1000/EEG.srate)/2);
times2saveidx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2saveidx(i)]=min(abs(EEG.times-times2save(i)));
end

electrodesidx(1) = find(strcmpi(electrodes4mi{1},{EEG.chanlocs.labels}));
electrodesidx(2) = find(strcmpi(electrodes4mi{2},{EEG.chanlocs.labels}));

% initialize outputs
entropy = zeros(3,length(times2save));
mi      = zeros(2,length(times2save));
nbins   = zeros(1,length(times2save));

for timei = 1:length(times2save)
    datax = EEG.data(electrodesidx(1),times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,:);
    datay = EEG.data(electrodesidx(2),times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,:);
    
    [mi(1,timei),entropy(:,timei),nbins(timei)] = mutualinformationx(datax,datay);
    [mi(2,timei),entropy(:,timei)             ] = mutualinformationx(datax,datay,70);
end

figure
set(gcf,'name',[ 'Mutual information between ' electrodes4mi{1} ' and ' electrodes4mi{2} ])

subplot(221)
plot(times2save,mi(1,:))
xlabel('Time (ms)'), ylabel('MI (bits)')
title('Variable bin length')
set(gca,'xlim',[times2save(1)-50 times2save(end)+50],'ylim',[min(mi(:))-.01 max(mi(:))+.01])

subplot(222)
plot(nbins,mi(1,:),'.')
xlabel('bin length'), ylabel('MI (bits)')
title('Bin length vs. MI')
set(gca,'ylim',[min(mi(:))-.01 max(mi(:))+.01])

subplot(223)
plot(times2save,mi(2,:))
xlabel('Time (ms)'), ylabel('MI (bits)')
title('Constant bin length')
set(gca,'xlim',[times2save(1)-50 times2save(end)+50],'ylim',[min(mi(:))-.01 max(mi(:))+.01])

%% Figure 29.11

% (Figure 10 must be generated before running this cell.)
frex = logspace(log10(4),log10(40),20);
baselinetime = [-500 -200];


% baseline from ms to idx
[junk,baseidx(1)] = min(abs(times2save-baselinetime(1)));
[junk,baseidx(2)] = min(abs(times2save-baselinetime(2)));

% specify convolution and wavelet info
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
% FFT of data
fft_EEG1 = fft(reshape(EEG.data(electrodesidx(1),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_EEG2 = fft(reshape(EEG.data(electrodesidx(2),:,:),1,EEG.pnts*EEG.trials),n_convolution);


% initialize outputs
mi   = zeros(2,length(frex),length(times2save));
ispc = zeros(length(frex),length(times2save));
powc = zeros(length(frex),length(times2save));

for fi=1:length(frex)
    % create wavelet and get its FFT
    fft_wavelet = fft(exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(4/(2*pi*frex(fi)))^2)),n_convolution);
    
    % convolution of each electrode with wavelet
    convres     = ifft(fft_wavelet.*fft_EEG1,n_convolution);
    analytic1   = reshape(convres(half_of_wavelet_size+1:end-half_of_wavelet_size),EEG.pnts,EEG.trials);
    convres     = ifft(fft_wavelet.*fft_EEG2,n_convolution);
    analytic2   = reshape(convres(half_of_wavelet_size+1:end-half_of_wavelet_size),EEG.pnts,EEG.trials);

    for timei = 1:length(times2save)
        datax = analytic1(times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,:);
        datay = analytic2(times2saveidx(timei)-timewindowidx:times2saveidx(timei)+timewindowidx,:);
        
        % compute MI
        mi(1,fi,timei) = mutualinformationx(log10(abs(datax).^2),log10(abs(datay).^2),50);
        mi(2,fi,timei) = mutualinformationx(angle(datax),angle(datay),20);
        
        % also compute ISPC-time for comparison
        ispc(fi,timei) = mean(abs(mean(exp(1i*(angle(datay)-angle(datax))),2)),1);
        
        % also compute power correlations-time
        dataxr = tiedrank(abs(datax));
        datayr = tiedrank(abs(datay));
        n = timewindowidx*2+1;
        powc(fi,timei) = mean( 1-6*sum((dataxr-datayr).^2)/(n*(n^2-1)) );
    end
    
    disp([ 'Finished frequency ' num2str(fi) ' out of ' num2str(length(frex)) ]);
end


figure, set(gcf,'name',[ 'Mutual information between ' electrodes4mi{1} ' and ' electrodes4mi{2} ])

for i=1:2
    % plot baseline-subtracted corrected MI
    subplot(2,2,i)
    contourf(times2save,frex,squeeze(mi(i,:,:))-repmat(mean(mi(i,:,baseidx(1):baseidx(2)),3)',1,length(times2save)),40,'linecolor','none')
    set(gca,'clim',[-.075 .075],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),6)))
end

% plot power correlations
subplot(223)
contourf(times2save,frex,powc-repmat(mean(powc(:,baseidx(1):baseidx(2)),2),1,length(times2save)),40,'linecolor','none')
set(gca,'clim',[-.2 .2],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),6)))

% plot phase synchronization
subplot(224)
contourf(times2save,frex,ispc-repmat(mean(ispc(:,baseidx(1):baseidx(2)),2),1,length(times2save)),40,'linecolor','none')
set(gca,'clim',[-.1 .1],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),6)))

%% Figure 29.12

time    = 0:.0001:1;
signal1 = sin(2*pi*10*time);
signal2 = -signal1;

lagz=1:10:1500;
milags=zeros(size(lagz));

for li=1:length(lagz)
    milags(li) = mutualinformationx(signal1,[signal2(lagz(li):end) signal2(1:lagz(li)-1)],15);
end

figure
plot(lagz/(1/mean(diff(time))),milags)
xlabel('Lag (seconds)'), ylabel('mutual information (bits)')


% now on real data (6 Hz power MI from figure 11)

time = -1:1/EEG.srate:1;
fft_wavelet = fft(exp(2*1i*pi*frex(5).*time) .* exp(-time.^2./(2*(4/(2*pi*frex(5)))^2)),n_convolution);

% convolution of each electrode with wavelet
convres = ifft(fft_wavelet.*fft_EEG1,n_convolution);
pow1    = log10(abs(reshape(convres(half_of_wavelet_size+1:end-half_of_wavelet_size),EEG.pnts,EEG.trials)).^2);
convres = ifft(fft_wavelet.*fft_EEG2,n_convolution);
pow2    = log10(abs(reshape(convres(half_of_wavelet_size+1:end-half_of_wavelet_size),EEG.pnts,EEG.trials)).^2);


onecycle = round(1000/frex(5));
onecycleidx = round(onecycle/(1000/EEG.srate));

lagz=-onecycleidx:onecycleidx;
milags=zeros(size(lagz));

for li=1:length(lagz)
    
    if lagz(li)<0
        milags(li) = mutualinformationx(pow1(1:end+lagz(li),:),pow2(-lagz(li)+1:end,:),30); % reverse sign for negative lags
    elseif lagz(li)==0
        milags(li) = mutualinformationx(pow1,pow2,30); % special case for no lag
    elseif lagz(li)>0
        milags(li) = mutualinformationx(pow1(lagz(li)+1:end,:),pow2(1:end-lagz(li),:),30);
    end
end

figure
plot(1000*lagz/(1/mean(diff(time))),milags)
xlabel([ electrodes4mi{1} ' leads ' electrodes4mi{2} ' ... Lag (ms) ... ' electrodes4mi{2} ' leads ' electrodes4mi{1} ]), ylabel('mutual information (bits)')
set(gca,'xlim',[-onecycle onecycle])

%% end.
