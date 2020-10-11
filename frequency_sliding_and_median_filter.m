%% This script shows you how to compute "frequency sliding,"
%  time-varying changes in the peak frequency of an oscillation.
%  This script corresponds to a paper:

% This script assumes that your data are in eeglab format. However, the
% eeglab toolbox is not necessary and it should be fairly straight forward
% to adapt this script to whatever format you use for your data.

% Questions? -> mikexcohen@gmail.com

%%

% Because the median filter involves computing the median many many times,
% I recommend installing a free Matlab toolbox called "nth_element"
% http://www.mathworks.com/matlabcentral/fileexchange/29453-nth-element
% The toolbox will speed the analysis but is not necessary

addpath('C:\Users\mxc\Documents\MATLAB\nth\')

%% setup preliminaries

% time points to save (speeds analysis by reducing the number of medians
% computed)
times2save = -300:25:1200;

% define boundaries for frequency bands
freq_bands = {
    [  4   8 ];
    [  8  12 ];
    [ 30  50 ];
    };


num_freq_bands = size(freq_bands,1);

% the next two lines are parameters for the median filter. The idea of the
% median filter is to compute the median several times using different
% sized windows. 'n_order' is the number of times to compute the median, and
% 'orders' is the size of the windows.
n_order = 10;
orders = linspace(10,400,n_order)/2; % recommended: 10 steps between 10 and 400 ms
% orders are cut in half to get n/2 datapoints before and n/2 datapoints
% after each center time point.
% Note that the widths are here specified in ms; they are
% converted to time points later on.

%% load in EEG data and do whatever processing/cleaning/etc is necessary



%% setup more parameters

% setup time indexing
times2saveidx = dsearchn(EEG.times',times2save');

% convert orders from ms to timepoints
orders = round( orders/(1000/EEG.srate) );

%% initialize output matrices

% time-frequency Hz and power
[tfhz,tfpw]  = deal(zeros(EEG.nbchan,num_freq_bands,length(times2save)));

%% loop through frequency bands

for fi=1:num_freq_bands
    
    %% filter data
    
    % apply a band-pass filter with 15% transition zones.
    
    trans_width    = .15;
    idealresponse  = [ 0 0 1 1 0 0 ];
    filtfreqbounds = [ 0 (1-trans_width)*freq_bands{fi,1}(1) freq_bands{fi,1}(1) freq_bands{fi,1}(2) freq_bands{fi,1}(2)*(1+trans_width) EEG.srate/2 ]/(EEG.srate/2);
    filt_order     = round(3*(EEG.srate/freq_bands{fi,1}(1)));
    filterweights  = firls(filt_order,filtfreqbounds,idealresponse);
    
    % this part does the actual filtering
    filterdata = zeros(size(EEG.data));
    for chani=1:EEG.nbchan
        filterdata(chani,:,:) = reshape( filtfilt(filterweights,1,double(reshape(EEG.data(chani,:,:),1,EEG.pnts*EEG.trials))) ,EEG.pnts,EEG.trials);
    end
    
    %% compute pre-filtered frequency sliding and band-specific power
    
    freqslide_prefilt = zeros(EEG.nbchan,EEG.pnts,EEG.trials);
    temppow = zeros(EEG.nbchan,length(times2save),EEG.trials);
    
    % loop over trials
    for triali=1:EEG.trials
        
        % get analytic signal via Hilbert transform
        temphilbert = hilbert(squeeze(filterdata(:,:,triali))')';
        
        % compute frequency sliding (note that this signal may be noisy;
        % median filtering is recommended before interpretation)
        freqslide_prefilt(:,1:end-1,triali) = diff(EEG.srate*unwrap(angle(temphilbert'),[],2)',1,2)/(2*pi);
        
        % and power (only at requested time points)
        temppow(:,:,triali) = abs(temphilbert(:,times2saveidx)).^2;
    end
    
    % power
    tfpw(:,fi,:) = squeeze(mean(temppow,3));
    
    %% apply median filter
    
    % temporary matrix involve in filtering
    phasedmed = zeros(EEG.nbchan,length(orders),length(times2save),EEG.trials);
    
    
    % median filter (only on requested time points to decrease computation time)
    for oi=1:n_order
        for ti=1:length(times2save)
            
            %% use compiled fast_median if available
            for triali=1:EEG.trials
                phasedmed(:,oi,ti,triali) = fast_median(freqslide_prefilt(:,max(times2saveidx(ti)-orders(oi),1):min(times2saveidx(ti)+orders(oi),EEG.pnts-1),triali)');
            end
            
            %% use 'manual median' otherwise
            %temp = sort(freqslide_prefilt(:,max(times2saveidx(ti)-orders(oi),1):min(times2saveidx(ti)+orders(oi),EEG.pnts-1),:),2);
            %phasedmed(:,oi,ti,:) = temp(:,floor(size(temp,2)/2)+1,:);
            
        end
    end
    
    % the final step is to take the mean of medians
    tfhz(:,fi,:) = mean(median(phasedmed,2),4);
    
end % end frequency band loop

%% end
