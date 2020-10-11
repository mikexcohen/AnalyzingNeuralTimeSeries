function psi = data2psiX(data,srate,freqbins,permtest)
% calculates phase slope index (PSI) as formulated in the paper:
%    Nolte G, Ziehe A, Nikulin VV, Schlogl A, Kramer N, Brismar T, Muller KR.
%    Robustly estimating the flow direction of information in complex physical systems.
%    Physical Review Letters. To appear.
%    (for further information:    http://doc.ml.tu-berlin.de/causality/ )
%
% Usage:
%   psi = data2psiX(data,srate,freqbins,permtest);
%
% Input:
%     data:  MxNxT matrix for M channels, N timepoints, and T trials
%    srate:  data sampling rate in Hz
% freqbins:  KxQ matrix containing frequency boundaries in rows K. 
% permtest:  Permutation test (boolean). If true, psi output values are in standardizes Z values
%               relative to a null hypothesis distribution
%
% Output:
%      psi:  phase-slope-index values. For M channels PSI is an MxM matrix if 
%               one frequency bin, or MxMxK if freqbins has K rows (with K>1).
%               psi(i,j) is the directed connectivity from channel i to
%               channel j, (e.g., channel i is the sender if psi(i,j) is
%               positive; channel i is the receiver if negative)
%
%
%  (This function was modified from the original by Mike X Cohen)


% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

% modified from the original dat2psi.m by mikexcohen@gmail.com

%% check inputs and initialize

if nargin<2
    help data2psiX
    error('Read help file.');
elseif nargin==2
    permtest = 0;
    freqbins = [];
elseif nargin==3
    permtest = 0;
end

n_permutes = 1000;


% data dimensions
[nchan,npnts,ntrials]=size(data);

% define frequencies
hz  = linspace(0,srate/2,floor(npnts/2)+1);
nhz = length(hz);

% initialize
cs  = zeros(nchan,nchan,nhz);
psi = zeros(nchan,nchan,size(freqbins,1));

%% compute cross-spectral density

datafft = fft(bsxfun(@times,data,hanning(npnts)'),[],2);

% trial-average cross-spectral density
for triali=1:ntrials
    for freqi=1:nhz
        
        tempfftdat = squeeze(datafft(:,freqi,triali))';
        cs(:,:,freqi) = cs(:,:,freqi) + tempfftdat'*conj(tempfftdat);
    end
end
cs = cs./triali;

%% compute PSI

for freqbini=1:size(freqbins,1)
    
    % find FFT indices of requested frequency bands
    freqidx = dsearchn(hz',freqbins(freqbini,1)) : dsearchn(hz',freqbins(freqbini,2));
    nfidx   = length(freqidx);
    
    if nfidx<4
        warning('There are fewer than four frequency bins. Consider using longer time windows or wider frequency bands.')
    end
    
    % temporary phase-frequency matrix from this frequency band
    pp = zeros(nchan,nchan,nfidx);
    
    for fi=1:nfidx
        pp(:,:,fi) = cs(:,:,freqidx(fi))./sqrt(diag(cs(:,:,freqidx(fi)))*diag(cs(:,:,freqidx(fi)))');
    end
    
    % average phase slope for each frequency band
    psi(:,:,freqbini) = sum(imag(conj(pp(:,:,1:end-1)).*pp(:,:,2:end)),3);
    
    
    %% optional permutation testing
    
    if permtest
        
        nulldist = zeros(nchan,nchan,n_permutes);
        for permi=1:n_permutes
            nulldist(:,:,permi) = sum(imag(conj(pp(:,:,randperm(nfidx))).*pp(:,:,randperm(nfidx))),3);
        end
        
        psi(:,:,freqbini) = ( psi(:,:,freqbini)-mean(nulldist,3) ) ./ std(nulldist,[],3);
    end
    
    % zero-out diagonals
    tmp = squeeze(psi(:,:,freqbini));
    tmp(logical(eye(nchan))) = 0;
    psi(:,:,freqbini) = tmp;
    
end

%% end

