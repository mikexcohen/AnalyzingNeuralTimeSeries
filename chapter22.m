%% Analyzing Neural Time Series Data
% Matlab code for Chapter 22
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% figure 22.2

% This figure was made by 'stepping-in' to the function laplacian_perrinX
% and then creating topographical maps of the Legendre polynomial.

%% Figure 22.3

% load sample EEG dataset
load sampleEEGdata

% compute inter-electrode distances
interelectrodedist=zeros(EEG.nbchan);
for chani=1:EEG.nbchan
    for chanj=chani+1:EEG.nbchan
        interelectrodedist(chani,chanj) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chanj).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chanj).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chanj).Z)^2);
    end
end

valid_gridpoints = find(interelectrodedist);

% extract XYZ coordinates from EEG structure
X = [EEG.chanlocs.X];
Y = [EEG.chanlocs.Y];
Z = [EEG.chanlocs.Z];



% create G and H matrices
[junk,G,H] = laplacian_perrinX(rand(size(X)),X,Y,Z,[],1e-6);


figure
subplot(221)
imagesc(G)
axis square
title('G')

subplot(222)
imagesc(H)
axis square
title('H')

subplot(223)
plot(interelectrodedist(valid_gridpoints),G(valid_gridpoints),'r.')
hold on
plot(interelectrodedist(valid_gridpoints),H(valid_gridpoints),'m.')
legend({'G';'H'})
set(gca,'ylim',[-.065 .065])
xlabel('Inter-electrode distances (mm)'), ylabel('G or H')
axis square

%% Figure 22.4

% In the book, figure 4 uses chan1 as Pz. The book also mentions that tha
% surface Laplacian will attenuate the impact of EOG artifacts. This can be
% simulated here by setting chan1 to FPz. 

% Note also that you'll need the eeglab function topoplot for this figure.

chan1 = 'pz';
chan2 = 'c4';
chan3 = 'c3';

eucdist1 = zeros(1,64);
eucdist2 = zeros(1,64);
eucdist3 = zeros(1,64);

chan1idx = strcmpi(chan1,{EEG.chanlocs.labels});
chan2idx = strcmpi(chan2,{EEG.chanlocs.labels});
chan3idx = strcmpi(chan3,{EEG.chanlocs.labels});

for chani=1:EEG.nbchan
    eucdist1(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan1idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan1idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan1idx).Z)^2 );
    eucdist2(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan2idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan2idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan2idx).Z)^2 );
    eucdist3(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan3idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan3idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan3idx).Z)^2 );
end

hi_spatfreq  = 2*exp(- (eucdist1.^2)/(2*95^2) ); 
lo_spatfreq  =   exp(- (eucdist2.^2)/(2*50^2) )  +  exp(- (eucdist3.^2)/(2*50^2) );
surf_lap_all = laplacian_perrinX(hi_spatfreq+lo_spatfreq,X,Y,Z);

figure
subplot(221)
topoplot(hi_spatfreq,EEG.chanlocs,'plotrad',.53);
title('Low spatial frequency feature')

subplot(222)
topoplot(lo_spatfreq,EEG.chanlocs,'plotrad',.53);
title('High spatial frequency features')

subplot(223)
topoplot(hi_spatfreq+lo_spatfreq,EEG.chanlocs,'plotrad',.53);
title('Low+high features')

subplot(224)
topoplot(surf_lap_all,EEG.chanlocs,'plotrad',.53);
title('Laplacian of low+high features')

%% another example similar to Figure 4

chan1 = 'cz';
chan2 = 'p5';


eucdist1 = zeros(1,64);
eucdist2 = zeros(1,64);

chan1idx = strcmpi(chan1,{EEG.chanlocs.labels});
chan2idx = strcmpi(chan2,{EEG.chanlocs.labels});

for chani=1:EEG.nbchan
    eucdist1(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan1idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan1idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan1idx).Z)^2 );
    eucdist2(chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan2idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan2idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan2idx).Z)^2 );
end

data2use = exp(- (eucdist1.^2)/(2*65^2) )  +  exp(- (eucdist2.^2)/(2*50^2) );

surf_lap = laplacian_perrinX(data2use,X,Y,Z);

figure
subplot(121)
topoplot(data2use,EEG.chanlocs,'plotrad',.53);
title('Spatially unfiltered')
subplot(122)
topoplot(surf_lap,EEG.chanlocs,'plotrad',.53);
title('surface Laplacian')

%% Figure 22.5

data2use = double(mean(EEG.data(:,321,:),3));

surf_lapN = laplacian_nola(X,Y,Z,data2use,100);
surf_lapP = laplacian_perrinX(data2use,X,Y,Z,[],1e-5); 
% note: try changing the smoothing parameter above (last input argument) to
% see the effects of the smoothing (lambda) parameter. Reasonable values
% are 1e-4 to 1e-6, and the default parameter is 1e-5.

figure
subplot(131)
topoplot(data2use,EEG.chanlocs,'plotrad',.53,'electrodes','off');
title('Raw data')

subplot(132)
topoplot(surf_lapN,EEG.chanlocs,'plotrad',.53,'electrodes','off');
title('Laplacian (Nunez book)')

subplot(133)
topoplot(surf_lapP,EEG.chanlocs,'plotrad',.53,'electrodes','off');
title('Laplacian (Perrin et al)')

disp([ 'Spatial correlation: r=' num2str(corr(surf_lapN,surf_lapP)) ])

%% Figure 22.6

% tic/toc are included in case you want to test the Perrin and New Orleans methods
timetest(1) = tic; lap_data = laplacian_perrinX(EEG.data,X,Y,Z); t(1) = toc;
timetest(2) = tic; lap_data2 = laplacian_nola(X,Y,Z,EEG.data); t(2) = toc;

times2plot = -100:100:800;


figure
for i=1:length(times2plot)
    
    % find time index
    [junk,timeidx] = min(abs(EEG.times-times2plot(i)));
    
    tempdata = double(squeeze(mean(EEG.data(:,timeidx,:),3)));
    
    % plot voltage map (spatially unfiltered)
    subplot(2,length(times2plot),i)
    topoplot(tempdata,EEG.chanlocs,'plotrad',.53,'maplimits',[-10 10],'electrodes','off');
    title([ 'Voltage, ' num2str(times2plot(i)) ' ms' ])
    
    % plot Laplacian map (spatially filtered)
    subplot(2,length(times2plot),i+length(times2plot))
    topoplot(laplacian_perrinX(tempdata,X,Y,Z),EEG.chanlocs,'plotrad',.53,'maplimits',[-40 40],'electrodes','off');
    title([ 'surface Laplacian, ' num2str(times2plot(i)) ' ms' ])
end

%% brief aside: 

% This figure shows that computing the Laplacian of the ERP is the same as
% computing the Laplacian of single trials and then taking the ERP. This is
% not surprising: the ERP is a linear transform of the single trials.

figure
subplot(121)
topoplot(laplacian_perrinX(double(squeeze(mean(EEG.data(:,321,:),3))),X,Y,Z),EEG.chanlocs,'plotrad',.53,'maplimits',[-40 40],'electrodes','off');

subplot(122)
topoplot(squeeze(double(mean(lap_data(:,321,:),3))),EEG.chanlocs,'plotrad',.53,'maplimits',[-40 40],'electrodes','off');

%% Figure 22.7

freq2use =   8; % Hz
time2use = 400; % ms

% FFT parameters
time          = -1:1/EEG.srate:1;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv2       = pow2(nextpow2(n_convolution));

% create wavelet, etc
wavelet_fft = fft(exp(2*1i*pi*freq2use.*time) .* exp(-time.^2./(2*(4/(2*pi*freq2use))^2))/freq2use,n_conv2);
half_of_wavelet_size = (length(time)-1)/2;



% initialize
allphases_pre = zeros(size(EEG.data));
allphases_lap = zeros(size(EEG.data));
ispc_pre = zeros(EEG.nbchan);
ispc_lap = zeros(EEG.nbchan);
timeidx  = dsearchn(EEG.times',time2use');

% get all phases
for chani=1:EEG.nbchan
    
    % first for nonspatially filtered data
    fft_data = fft(reshape(EEG.data(chani,:,:),1,EEG.pnts*EEG.trials),n_conv2);
    conv_res = ifft(wavelet_fft.*fft_data,n_conv2);
    conv_res = conv_res(1:n_convolution);
    conv_res = conv_res(half_of_wavelet_size+1:end-half_of_wavelet_size);
    % collect analytic signal
    allphases_pre(chani,:,:) = reshape(conv_res,EEG.pnts,EEG.trials);

    % then for laplacian filtered data
    fft_data = fft(reshape(lap_data(chani,:,:),1,EEG.pnts*EEG.trials),n_conv2);
    conv_res = ifft(wavelet_fft.*fft_data,n_conv2);
    conv_res = conv_res(1:n_convolution);
    conv_res = conv_res(half_of_wavelet_size+1:end-half_of_wavelet_size);
    % collect analytic signal
    allphases_lap(chani,:,:) = reshape(conv_res,EEG.pnts,EEG.trials);
end

% compute synchrony
for chani=1:EEG.nbchan
    for chanj=chani+1:EEG.nbchan
        
        % cross-spectral density
        cd = squeeze(allphases_pre(chani,timeidx,:).*conj(allphases_pre(chanj,timeidx,:)));
        ispc_pre(chani,chanj) = abs(mean(exp(1i*angle(cd))));
        
        cd = squeeze(allphases_lap(chani,timeidx,:).*conj(allphases_lap(chanj,timeidx,:)));
        ispc_lap(chani,chanj) = abs(mean(exp(1i*angle(cd))));
    end
end

% mirror connectivity matrices
ispc_pre = ispc_pre + ispc_pre' + eye(EEG.nbchan);
ispc_lap = ispc_lap + ispc_lap' + eye(EEG.nbchan);


figure
subplot(121)
plot(interelectrodedist(valid_gridpoints),ispc_pre(valid_gridpoints),'.')
xlabel('Electrode distances (mm)'), ylabel('ISPC')
title([ 'Spatially unfiltered ISPC at ' num2str(freq2use) ' Hz' ])
set(gca,'ylim',[0 1])
r=corr(interelectrodedist(valid_gridpoints),ispc_pre(valid_gridpoints),'type','s');
legend([ 'R^2 = ' num2str(r^2) ])
axis square

subplot(122)
plot(interelectrodedist(valid_gridpoints),ispc_lap(valid_gridpoints),'.')
xlabel('Electrode distances (mm)'), ylabel('ISPC')
title([ 'Laplacian filtered ISPC at ' num2str(freq2use) ' Hz' ])
set(gca,'ylim',[0 1])
r=corr(interelectrodedist(valid_gridpoints),ispc_lap(valid_gridpoints),'type','s');
legend([ 'R^2 = ' num2str(r^2) ])
axis square

%% Figure 22.8

figure

subplot(121)
topoplot(ispc_pre(48,:),EEG.chanlocs,'maplimits',[0 .8],'plotrad',.53);
title([ 'ISPC_r_a_w at ' num2str(time2use) ' ms, ' num2str(freq2use) ' Hz' ])

subplot(122)
topoplot(ispc_lap(48,:),EEG.chanlocs,'maplimits',[0 .8],'plotrad',.53);
title([ 'ISPC_l_a_p at ' num2str(time2use) ' ms, ' num2str(freq2use) ' Hz' ])

%% end.
