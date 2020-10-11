function [entropy,fd_bins] = entropyx(x,fd_bins)
% ENTROPYX   Compute entropy
%  
%   Inputs:
%       x       :  data matrix
%
%   Optional inputs:
%       bins    :  number of bins to use for distribution discretization
%
%   Outputs:
%       entropy :  entropy of input variable x
%       nbins   :  number of bins used for discretization
%                  (based on Freedman-Diaconis rule)
%
%  Mike X Cohen (mikexcohen@gmail.com)

%% determine the optimal number of bins for each variable

% vectorize in the case of matrices
x=x(:);

if nargin<2 || isempty(fd_bins)
    n            = length(x);
    maxmin_range = max(x)-min(x);
    fd_bins      = ceil(maxmin_range/(2.0*iqr(x)*n^(-1/3))); % Freedman-Diaconis
end

%% compute entropies

% recompute entropy with optimal bins for comparison
hdat1 = hist(x,fd_bins);
hdat1 = hdat1./sum(hdat1);

% convert histograms to probability values
entropy = -sum(hdat1.*log2(hdat1+eps));

%%
