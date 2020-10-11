function imfs=emdx(data,varargin)
% EMDX - perform empirical mode decomposition on channels x time x trial data
%
% Usage:
%  imfs = emdx(data[,maxorder,standdev]);
%
% Inputs:
%   data     = [chans time trials] input data (can also input [time trials]
%               or [time] data)
%   maxorder = (optional) maximum order for EMD (default is 30)
%   standdev = (optional) standard deviation for sifting stopping critia (default is .5)
%
% Outputs:
%    imfs    = [chans modes time trials] matrix of intrinsic mode functions
%
%
% Notes:
%    1) It may be useful to low-pass filter the data somewhat above the highest
%       frequency of interest.
%    2) Instantaneous frequency (in Hz) can be estimated from the phase derivative:
%       f = srate*diff(unwrap(angle(hilbert(imfs(1,1,:,1)))))/(2*pi);
%    3) EMD is slow (e.g., compared to time-frequency decomposition). Use data from
%       restricted time windows, and downsample the data if possible.

% mikexcohen@gmail.com

switch nargin
    case 0
        help emdx
        error('No inputs. See help file.');
    case 1
        maxorder = 30;
        maxstd   = .5;
    case 2
        maxorder = varargin{1};
        maxstd   = .5;
    case 3
        maxorder = varargin{1};
        maxstd   = varargin{2};
    otherwise
        help emdx
        error('Too many inputs. See help file.');
end

maxiter = 1000;

%% chaek data size and adjust if necessary

if isvector(data) % one trial
    data = reshape(data,[1 numel(data) 1]);
elseif ismatrix(data) % assume time X trials
    data = reshape(data,[1 size(data)]);
end

[nchans,npnts,ntrials] = size(data);
time = 1:npnts;

% initialize
imfs = zeros(nchans,maxorder,npnts,ntrials);
imforders = zeros(1,ntrials);

%% use griddedInterpolant if exists (much faster than interp1)

% griddedInterpolat should always be preferred when your Matlab version supports it.

if exist('griddedInterpolant','file')
    dofast = true;
else dofast = false;
end

%% loop over channels

for chani=1:nchans
    
    %% loop over trials
    
    for triali=1:ntrials
        
        % data from this trial (must be a row vector)
        imfsignal = squeeze(data(chani,:,triali));
        
        %% loop over IMF order
        
        imforder = 1;
        stop     = false;
        
        while ~stop
            
            %% iterate over sifting process
            
            % initializations
            standdev = 10;
            numiter  = 0;
            signal   = imfsignal;
            
            % "Sifting" means iteratively identifying peaks/troughs in the
            % signal, interpolating across peaks/troughs, and then recomputing
            % peaks/troughs in the interpolated signal. Sifting ends when
            % variance between interpolated signal and previous sifting
            % iteration is minimized.
            while standdev>maxstd && numiter<maxiter
                
                % identify local min/maxima
                localmin  = [1 find(diff(sign(diff(signal)))>0)+1 npnts];
                localmax  = [1 find(diff(sign(diff(signal)))<0)+1 npnts];
                
                % create envelopes as cubic spline interpolation of min/max points
                if dofast % faster method, but works only on recent Matlab versions
                    FL = griddedInterpolant(localmin(:),signal(localmin)','spline');
                    FU = griddedInterpolant(localmax(:),signal(localmax)','spline');
                    env_lower = FL(time);
                    env_upper = FU(time);
                else % backup method, just in case
                    env_lower = interp1(localmin,signal(localmin),time);
                    env_upper = interp1(localmax,signal(localmax),time);
                end
                
                % compute residual and standard deviation
                prevsig   = signal;
                signal    = signal - (env_lower+env_upper)./2;
                standdev  = sum( ((prevsig-signal).^2) ./ (prevsig.^2+eps) ); % eps prevents NaN's
                
                % not too many iterations
                numiter = numiter+1;
                
            end % end sifting
            
            % imf is residual of signal and min/max average (already redefined as signal)
            imfs(chani,imforder,:,triali) = signal;
            imforder = imforder+1;
            
            %% residual is new signal
            
            imfsignal = imfsignal-signal;
            
            %% stop when few points are left
            
            if numel(localmax)<5 || imforder>maxorder
                stop=true;
            end
            
        end % end imf for this trial
        
        imforders(triali) = imforder;
        
    end % end trials
end % end channels

%% clean up

% if no max imfs requested, cut down to size
if nargin==1
    imfs(:,min(imforders):end,:,:) = [];
end

imfs = squeeze(imfs); % in case no trials and/or channels

%% end.
% 
% figure
% hold on
% colorz='rrrrbgmkyrbgmky';
% 
% allfs=zeros(8,numel(imf_times)-1);
% origf=allfs;
% 
% for fi=1:6
%     
%     f  = squeeze(diff(unwrap(angle(hilbert(imfs{end}(fi,:,:))))*ALLEEG(1).srate)/(2*pi));
%     zf = zscore(diff(f));
%     for triali=1:size(f,2)
%         f(find(abs(zf(:,triali))>3)+1,triali)=NaN;
%     end
% 
%     f = nanmean(f,2);
%     
%     allfs(fi,:) = (f*100)./mean(f(dsearchn(imf_times',-400):dsearchn(imf_times',-100)))' - 100;
%     origf(fi,:) = f;
%     
%     plot(imf_times(1:end-1),allfs,colorz(fi))
%     set(gca,'xlim',[-200 1000])
% 
% end
