%% Analyzing Neural Time Series Data
% Matlab code for Chapter 9
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 9.1a

% load EEG data
load sampleEEGdata.mat

% plot a few trials from one channel...

% specify the label of the channel to plot
which_channel_to_plot = 'fcz';
% and find the index (channel number) of that label
channel_index = strcmpi(which_channel_to_plot,{EEG.chanlocs.labels});

x_axis_limit = [-200 1000]; % in ms

num_trials2plot = 12;


figure
set(gcf,'Name',[ num2str(num_trials2plot) ' random trials from channel ' which_channel_to_plot ],'Number','off')
for i=1:num_trials2plot
    
    % figure out how many subplots we need
    subplot(ceil(num_trials2plot/ceil(sqrt(num_trials2plot))),ceil(sqrt(num_trials2plot)),i)
    
    % pick a random trial (using randsample, which is in the stats toolbox)
    random_trial_to_plot = randsample(EEG.trials,1);
    % if you don't have the stats toolbox, use the following two lines:
    %random_trial_to_plot = randperm(EEG.trials);
    %random_trial_to_plot = random_trial_to_plot(1);
    
    % plot trial and specify x-axis and title
    plot(EEG.times,squeeze(EEG.data(channel_index,:,random_trial_to_plot)));
    set(gca,'xlim',x_axis_limit,'ytick',[])
    title([ 'Trial ' num2str(random_trial_to_plot) ])
end

%% Figure 9.1b

figure
% plot all trials
plot(EEG.times,squeeze(EEG.data(channel_index,:,:)),'y')

hold on
% plot ERP (simply the average time-domain signal)
plot(EEG.times,squeeze(mean(EEG.data(channel_index,:,:),3)),'k','linew',2)
set(gca,'xlim',[-300 1000],'ylim',[-60 60],'ydir','reverse')


% now plot only the ERP
figure
plot(EEG.times,squeeze(mean(EEG.data(channel_index,:,:),3))) % Note the "3" as second input to "mean"; this takes the average of the 3rd dimension.

% plot lines indicating baseline activity and stim onset
hold on
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0],get(gca,'ylim'),'k:')

% add axis labels and title
xlabel('Time (ms)')
ylabel('\muV') % note that matlab interprets "\mu" as the Greek character for micro
title([ 'ERP (average of ' num2str(EEG.trials) ' trials) from electrode ' EEG.chanlocs(channel_index).labels ])

% plot upside down, following ERP convention
set(gca,'ydir','reverse')
% axis ij % this does the same as the previous trial

% below is some advanced but flexible code to change the x-axis label
set(gca,'xlim',[-300 1000])
xticklabel=cellstr(get(gca,'xticklabel'));
xticklabel{str2double(xticklabel)==0}='stim';
set(gca,'xticklabel',xticklabel)

%% Figure 9.2

% pick a channel
chan2plot = 'p7';

% compute ERP
erp = double(squeeze(mean(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),3)));

% low-pass filter data (requires signal processing toolbox) 
% you'll learn about how filtering works, and what this code means, in chapter 14. 
nyquist          = EEG.srate/2;
transition_width = 0.15; % percent

% first, filter from 0-40
filter_cutoff = 40; % Hz
ffrequencies  = [ 0 filter_cutoff filter_cutoff*(1+transition_width) nyquist ]/nyquist;
idealresponse = [ 1 1 0 0 ];
filterweights = firls(100,ffrequencies,idealresponse);
erp_0to40     = filtfilt(filterweights,1,double(erp));

% next, filter from 0-10
filter_cutoff = 10; % Hz
ffrequencies  = [ 0 filter_cutoff filter_cutoff*(1+transition_width) nyquist ]/nyquist;
idealresponse = [ 1 1 0 0 ];
filterweights = firls(100,ffrequencies,idealresponse);
erp_0to10     = filtfilt(filterweights,1,double(erp));

% finally, filter from 5-15
filter_low    =  5; % Hz
filter_high   = 15; % Hz
ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(round(3*(EEG.srate/filter_low)),ffrequencies,idealresponse);
erp_5to15     = filtfilt(filterweights,1,double(erp));



% now plot all filtered ERPs
figure
plot(EEG.times,erp,'k')
hold on
plot(EEG.times,erp_0to40,'b','linew',2)
plot(EEG.times,erp_0to10,'r')
plot(EEG.times,erp_5to15,'m')

set(gca,'xlim',[-200 1200],'ydir','r')
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
title([ 'ERP from electrode ' chan2plot ])
legend({'None';'0-40';'0-10';'5-15'})

%% Figure 9.3

figure

% Butterfly plot
subplot(211)
plot(EEG.times,squeeze(mean(EEG.data,3)))
set(gca,'xlim',[-200 1000],'ydir','reverse')
xlabel('Time (ms)'), ylabel('\muV')
title('ERP from all sensors')

% topographical variance plot
subplot(212)
plot(EEG.times,squeeze(var( mean(EEG.data,3) ))) % note: change var to std for global field power
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)'), ylabel('var(\muV)')
title('Topographical variance')

%% figure 9.4

% topoplot with colored dots vs. interpolated surface
% topoplot is a function in the eeglab toolbox, which you can download for
% free from the internet.

figure
subplot(121)
% To get the following line to work, you need to modify the eeglab function 'topoplot'
% Replace the first line of code (the function definition) with:
% function [handle,Zi,grid,Xi,Yi,x,y,pltchans] = topoplot(Values,loc_file,varargin)
[~,~,~,~,~,xi,yi,pltchans] = topoplot(zeros(64,1),EEG.chanlocs,'electrodes','off','plotrad',.53);

hold on
c = -squeeze(mean(EEG.data(:,300,:),3));
c = c-min(c); % these lines scale the variable c
c = c./max(c);
c = repmat(c,1,3);
for i=1:length(xi)
    hold on
    plot3(yi(i),xi(i),.1,'o','markerfacecolor',c(i,:),'markersize',5,'markeredge','none');
    % here, plot3 is used to plot the electrode on top of the topomap
    % (compare with the plot function) (a trick I learned from topoplot.m)
end

subplot(122)
c = -squeeze(mean(EEG.data(:,300,:),3));
c = c-min(c);
c = c./max(c);
c = repmat(c,1,3);
topoplot(double(c(:,1)),EEG.chanlocs,'plotrad',.53,'electrodes','off','numcontour',0);

colormap gray

%% introduction to topographical plotting

timepoint2plot = 100; % in ms
trial2plot = randsample(EEG.trials,1);

color_limit = 20; % more-or-less arbitrary, but this is a good value

% convert time point from ms to index
[junk,timepointidx] = min(abs(EEG.times-timepoint2plot));


% First step is to get X and Y coordinates of electrodes
% These must be converted to polar coordinates
th=pi/180*[EEG.chanlocs.theta];
[electrode_locs_X,electrode_locs_Y] = pol2cart(th,[EEG.chanlocs.radius]);

% interpolate to get nice surface
interpolation_level = 100; % you can try changing this number, but 100 is probably good
interpX = linspace(min(electrode_locs_X),max(electrode_locs_X),interpolation_level);
interpY = linspace(min(electrode_locs_Y),max(electrode_locs_Y),interpolation_level);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);

% let's look at these matrices
figure
subplot(121)
imagesc(gridX)

subplot(122)
imagesc(gridY)


% now interpolate the data on a 2D grid
interpolated_EEG_data = griddata(electrode_locs_Y,electrode_locs_X,double(squeeze(EEG.data(:,timepointidx,trial2plot))),gridX,gridY);

figure
set(gcf,'number','off','name',[ 'Topographical data from trial ' num2str(trial2plot) ', time=' num2str(round(EEG.times(timepointidx))) ]);
subplot(121)
contourf(interpY,interpX,interpolated_EEG_data,100,'linecolor','none');
axis square
set(gca,'clim',[-color_limit color_limit],'xlim',[min(interpY) max(interpY)]*1.1,'ylim',[min(interpX) max(interpX)]*1.1)
title('Interpolated data in matlab')

subplot(122)
topoplot(double(squeeze(EEG.data(:,timepointidx,trial2plot))),EEG.chanlocs,'maplimits',[-color_limit color_limit]); % eeglab's topoplot function
title('eeglab ''topoplot'' function')

figure
set(gcf,'name','a landscape of cortical electrophysiological dynamics')
surf(interpY,interpX,interpolated_EEG_data);
xlabel('left-right of scalp'), ylabel('anterior-posterior of scalp'), zlabel('\muV')
shading interp, axis tight
set(gca,'clim',[-color_limit color_limit],'xlim',[min(interpY) max(interpY)]*1.1,'ylim',[min(interpX) max(interpX)]*1.1)
rotate3d on, view(0,90)

%% Figure 9.5

% find indices for requested plot times
times2plot = dsearchn(EEG.times',(-100:50:600)');

figure
for i=1:length(times2plot)
    subplot(3,5,i)
    
    % extract EEG data and replace FC4 data with noise
    eegdata2plot = double(squeeze(mean(EEG.data(:,times2plot(i),:),3)));
    eegdata2plot(strcmpi('fc4',{EEG.chanlocs.labels})) = randn*10;

    topoplot(eegdata2plot,EEG.chanlocs,'maplimits',[-8 8]);
    title([ num2str(round(EEG.times(times2plot(i)))) ' ms' ])
end


%% Figure 9.6

useRTs = true; % or false

% get RTs from each trial, to use for sorting trials. In this experiment,
% the RT was always the first event after the stimulus (the time=0 event).
% Normally, you should build in exceptions in case there was no response or
% another event occured between the stimulus and response. This was already
% done for the current dataset. 

rts = zeros(size(EEG.epoch));
for ei=1:length(EEG.epoch)
    
    % first, find the index at which the time=0 event occurs
    time0event = find(cell2mat(EEG.epoch(ei).eventlatency)==0);
    
    % then get reaction time
    rts(ei) = EEG.epoch(ei).eventlatency{time0event+1};
end

% now find the sorting indices for RTs
if useRTs
    [dontneed,rts_idx] = sort(rts);
else
    [dontneed,rts_idx] = sort(squeeze(EEG.data(47,334,:)));
end

% now plot
figure
imagesc(EEG.times,1:EEG.trials,squeeze(EEG.data(47,:,rts_idx))');
set(gca,'clim',[-30 30],'xlim',[-200 1200],'ydir','n')

% also plot the RTs on each trial
if useRTs
    hold on
    plot(rts(rts_idx),1:EEG.trials,'k','linew',3)
end

%% end.
