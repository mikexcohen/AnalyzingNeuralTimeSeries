%% Analyzing Neural Time Series Data
% Matlab code for Chapter 28
% Mike X Cohen
%
% This code accompanies the book, titled "Analyzing Neural Time Series Data"
% (MIT Press). Using the code without following the book may lead to confusion,
% incorrect data analyses, and misinterpretations of results.
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code.

%% Figure 28.1

figure

% - Stationary-like process - %

x = randn;
for i=2:30
    x(i) = exp(cos(pi*x(i-1))) + randn;
end

subplot(221)
plot(x,'mo-','linewid',1,'markerface','k','markersize',6)
set(gca,'xlim',[0 31])
legend('x_t = e^c^o^s^(^\pi^xt-1^) + randn')
title('Stationary autoregressive process')

x=randn(2,1);
for i=3:30
    x(i) = .2*x(i-1) - .4*x(i-2) + randn;
end

subplot(222)
plot(x,'mo-','linewid',1,'markerface','k','markersize',6)
set(gca,'xlim',[0 31])
legend('x_t = .2x_t_-_1 - .4x_t_-_2 + randn')
title('Stationary autoregressive process')



% - Non-stationary process - %

x=1;
for i=2:30
    x(i) = 1.1*x(i-1) + randn;
end

subplot(223)
plot(x,'kp-','linewid',1,'markerface','g','markersize',10)
set(gca,'xlim',[0 31])
title('Univariate autoregression')
legend('x_t = 1.1*x_t_-_1 + randn')
title('Non-stationary autoregressive process')

x=[1 1.5];
for i=3:30
    x(i) = 1.2*x(i-2) + -.3*x(i-1) + randn;
end

subplot(224)
plot(x,'kp-','linewid',1,'markerface','g','markersize',10)
set(gca,'xlim',[0 31])
legend('x_t = -0.3*x_t_-_1 + 1.2*x_t_-_2 + randn')
title('Non-stationary autoregressive process')

%% Figure 28.2

% - Stationary-like process - %

% define X
x=randn(2,1);
for i=3:30
    x(i) = .2*x(i-1) - .4*x(i-2) + randn;
end

% define y
y=rand(2,1); % random initial conditions
for i=3:length(x)
    y(i) = .25*y(i-1) - .8*x(i-2) + 1.5*x(i-1) + randn;
end

figure
subplot(221)
plot(x,'mo-','linewid',1,'markerface','k','markersize',6)
hold on
plot(y,'go-','linewid',1,'markerface','k','markersize',6)
set(gca,'xlim',[0 31])
title('Stationary bivariate autoregression')
legend({'x_t = .2x_t_-_1 - .4x_t_-_2 + randn','y=0.25*y_t_-_1 - 0.8*x_t_-_2 + 1.5*x_t_-_1 + randn'})



% - Non-stationary process - %

% define X
x=[1 1.5];
for i=3:30
    x(i) = 1.2*x(i-2) + -.3*x(i-1) + randn;
end

% define y
y=rand(2,1); % random initial conditions
for i=3:length(x)
    y(i) = .25*y(i-1) - 1.2*x(i-2) + randn;
end

subplot(222)
plot(x,'mo-','linewid',1,'markerface','k','markersize',6)
hold on
plot(y,'go-','linewid',1,'markerface','k','markersize',6)
set(gca,'xlim',[0 31])
title('Non-stationary bivariate autoregression')
legend({'x_t = -.3x_t_-_1 - 1.2x_t_-_2 + randn','y=0.25*y_t_-_1 - 0.8*x_t_-_2 + 1.5*x_t_-_1 + randn'})

%% Figure 28.3 (note: this cell takes a while to run)

% load sample EEG data
load sampleEEGdata

% define channels for granger prediction
chan1name = 'o1';
chan2name = 'f5';

% Granger prediction parameters
timewin = 200; % in ms
order   =  27; % in ms

% temporal down-sample results (but not data!)
times2save = -400:20:1000; % in ms


% convert parameters to indices
timewin_points = round(timewin/(1000/EEG.srate));
order_points   = round(order/(1000/EEG.srate));

% find the index of those channels
chan1 = find(strcmpi(chan1name,{EEG.chanlocs.labels}));
chan2 = find(strcmpi(chan2name,{EEG.chanlocs.labels}));

% remove ERP from selected electrodes to improve stationarity
eegdata = bsxfun(@minus,EEG.data([chan1 chan2],:,:),mean(EEG.data([chan1 chan2],:,:),3));


% convert requested times to indices
times2saveidx = dsearchn(EEG.times',times2save');

% initialize
[x2y,y2x] = deal(zeros(1,length(times2save))); % the function deal assigns inputs to all outputs
bic = zeros(length(times2save),15); % Bayes info criteria (hard-coded to order=15)

for timei=1:length(times2save)

    % data from all trials in this time window
    tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));

    % detrend and zscore all data
    for triali=1:size(tempdata,3)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));

        % At this point with real data, you should check for stationarity
        % and possibly discard or mark data epochs that are extreme stationary violations.
    end

    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);

    % fit AR models (model estimation from bsmart toolbox)
    [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
    [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
    [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);

    % time-domain causal estimate
    y2x(timei)=log(Ex/E(1,1));
    x2y(timei)=log(Ey/E(2,2));

    % test BIC for optimal model order at each time point
    % (this code is used for the following cell)
    for bici=1:size(bic,2)
        % run model
        [Axy,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
        % compute Bayes Information Criteria
        bic(timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
    end
end

% draw lines
figure
plot(times2save,x2y)
hold on
plot(times2save,y2x,'r')
legend({[ 'GP: ' chan1name ' -> ' chan2name ];[ 'GP: ' chan2name ' -> ' chan1name ]})
title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ])
xlabel('Time (ms)')
ylabel('Granger prediction estimate')

%% Figure 28.4

% see "bici" for-loop above for code to compute BIC

figure

subplot(121)
plot((1:size(bic,2))*(1000/EEG.srate),mean(bic,1),'--.')
xlabel('Order (converted to ms)')
ylabel('Mean BIC over all time points')

[bestbicVal,bestbicIdx]=min(mean(bic,1));
hold on
plot(bestbicIdx*(1000/EEG.srate),bestbicVal,'mo','markersize',15)

title([ 'Optimal order is ' num2str(bestbicIdx) ' (' num2str(bestbicIdx*(1000/EEG.srate)) ' ms)' ])

subplot(122)
[junk,bic_per_timepoint] = min(bic,[],2);
plot(times2save,bic_per_timepoint*(1000/EEG.srate),'--.')
xlabel('Time (ms)')
ylabel('Optimal order (converted to ms)')
title('Optimal order (in ms) at each time point')

%% Figure 28.5

min_freq = 10; % in Hz, using minimum of 10 Hz because of 200-ms window
max_freq = 40;

order_points = 15;

frequencies = logspace(log10(min_freq),log10(max_freq),15);

% initialize
tf_granger=zeros(2,length(frequencies),length(times2save));


for timei=1:length(times2save)

    % data from all trials in this time window
    % (note that the ERP-subtracted data are used)
    tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));

    % detrend and zscore all data
    for triali=1:size(tempdata,3)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));

        % At this point with real data, you might want to check for stationarity
        % and possibly discard or mark data epochs that are non-stationary.
    end

    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);

    % fit AR models
    [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
    [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
    [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);

    % code below is adapted from bsmart toolbox function pwcausal.m
    % corrected covariance
    eyx = E(2,2) - E(1,2)^2/E(1,1);
    exy = E(1,1) - E(2,1)^2/E(2,2);
    N = size(E,1);

    for fi=1:length(frequencies)

        % transfer matrix (note the similarity to Fourier transform)
        H = eye(N);
        for m = 1:order_points
            H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/EEG.srate);
        end

        Hi = inv(H);
        S  = H\E*Hi'/EEG.srate;

        % granger prediction per frequency
        tf_granger(1,fi,timei) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/EEG.srate) );
        tf_granger(2,fi,timei) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/EEG.srate) );
    end
end % end time loop


% - plot - %
figure, set(gcf,'Name',[ 'Granger prediction between electrodes' chan1name ' and ' chan2name '.' ]);
subplot(211)
contourf(times2save,frequencies,squeeze(tf_granger(1,:,:)),40,'linecolor','none')
set(gca,'clim',[0 .025])
colorbar
title([ chan1name '->' chan2name '; ''raw'' units' ])

subplot(212)
contourf(times2save,frequencies,squeeze(tf_granger(2,:,:)),40,'linecolor','none')
set(gca,'clim',[0 .025])
colorbar
title([ chan2name '->' chan1name '; ''raw'' units' ])


%% figure 28.6

% plot of cycles per frequency
figure
for fi=1:length(frequencies)
    subplot(4,4,fi)
    plot((1:order_points)*(1000/EEG.srate),real(exp(-1i*(1:order_points)*2*pi*frequencies(fi)/EEG.srate)))
    set(gca,'xlim',[.5 1.015*order_points]*1000/EEG.srate,'ylim',[-1.1 1.1])
    title([ num2str(round(frequencies(fi))) ' Hz' ])
end
xlabel('Time (ms)')

%% Figure 28.7

figure

electrode2plot = strcmpi('o2',{EEG.chanlocs.labels});
erp = squeeze(mean(EEG.data(electrode2plot,:,:),3));
subplot(211)
plot(EEG.times,squeeze(EEG.data(electrode2plot,:,1)))
hold on
plot(EEG.times,squeeze(EEG.data(electrode2plot,:,1))-erp,'r')
plot(EEG.times,erp,'k')
set(gca,'xlim',[-200 1200],'ydir','r')
legend({'single trial';'single trial - ERP';'ERP'})
title([ 'Data from electrode ' EEG.chanlocs(electrode2plot).labels ]);

subplot(212)
electrode2plot = strcmpi('afz',{EEG.chanlocs.labels});
erp = squeeze(mean(EEG.data(electrode2plot,:,:),3));
plot(EEG.times,squeeze(EEG.data(electrode2plot,:,1)))
hold on
plot(EEG.times,squeeze(EEG.data(electrode2plot,:,1))-erp,'r')
plot(EEG.times,erp,'k')
set(gca,'xlim',[-200 1200],'ydir','r')
legend({'single trial';'single trial - ERP';'ERP'})
title([ 'Data from electrode ' EEG.chanlocs(electrode2plot).labels ]);

%% Figure 28.8

% baseline time window
baseline_period = [ -400 -100 ];

% convert to indices
[junk,baseidx(1)] = min(abs(times2save-baseline_period(1)));
[junk,baseidx(2)] = min(abs(times2save-baseline_period(2)));

% plot as % changes from baseline
figure
plot(times2save,100*(x2y-mean(x2y(baseidx(1):baseidx(2))))/mean(x2y(baseidx(1):baseidx(2))))
hold on
plot(times2save,100*(y2x-mean(y2x(baseidx(1):baseidx(2))))/mean(y2x(baseidx(1):baseidx(2))),'r')
legend({[ 'GC: ' chan1name ' -> ' chan2name ];[ 'GC: ' chan2name ' -> ' chan1name ]})
xlabel('Time (ms)')
ylabel('Granger prediction estimate (% change from baseline)')



% also convert time-frequency domain to percent change

tf_grangerPC = tf_granger;
for i=1:2
    meangranger = mean(tf_grangerPC(i,:,baseidx(1):baseidx(2)),3);
    % fancy bsxfun code for subtracting the mean and dividing by the mean
    tf_grangerPC(i,:,:) = 100*(bsxfun(@rdivide, bsxfun(@minus,tf_grangerPC(i,:,:),meangranger) ,meangranger));
end

figure
subplot(211)
contourf(times2save,frequencies,squeeze(tf_grangerPC(1,:,:)),40,'linecolor','none')
set(gca,'clim',[-100 100])
colorbar
title([ chan1name '->' chan2name '; percent change' ])

subplot(212)
contourf(times2save,frequencies,squeeze(tf_grangerPC(2,:,:)),40,'linecolor','none')
set(gca,'clim',[-150 150])
colorbar
title([ chan2name '->' chan1name '; percent change' ])

%% Figure 28.9

% generate two chi-square distributed random numbers
% (note: requires Matlab statistics toolbox)
d1 = chi2rnd(2,1,1000);
d2 = chi2rnd(2,1,1000);

% get histograms
[y1,x1]=hist(d1,50);
[y2,x2]=hist(d2,50);
[y3,x3]=hist(d1-d2,50);

figure
subplot(221)
plot(x1,y1,'k')
hold on
plot(x2,y2,'r')

subplot(223)
plot(x3,y3)
set(gca,'xlim',[-10 10])



% once more, with new distributions
d1 = chi2rnd(7,1,1000);
d2 = chi2rnd(7,1,1000);

% get histograms
[y1,x1]=hist(d1,50);
[y2,x2]=hist(d2,50);
[y3,x3]=hist(d1-d2,50);

subplot(222)
plot(x1,y1,'k')
hold on
plot(x2,y2,'r')

subplot(224)
plot(x3,y3)
set(gca,'xlim',[-20 20])

%% end.
