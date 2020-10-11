function [mi entropy fd_bins] = mutualinformationx(x,y,fd_bins,permtest)
% MUTUALINFORMATIONX   Compute mutual information between two vectors
%  
%   Inputs:
%       x,y     :  data matrices of equal size
%
%   Optional inputs:
%       bins    :  number of bins to use for distribution discretization
%      permtest :  perform permutation test and return mi in standard-Z values
%
%   Outputs:
%       mi      :  mutual information in bits
%       entropy :  entropy of x, y, and joint
%       nbins   :  number of bins used for discretization
%                  (based on Freedman-Diaconis rule)
%
%  Mike X Cohen (mikexcohen@gmail.com)

if nargin<2, error('Specify two inputs.'); end
if length(x)~=length(y), error('X and Y must have equal length'); end

%% determine the optimal number of bins for each variable

% vectorize in the case of matrices
x=x(:); y=y(:);

if nargin<3 || isempty(fd_bins)
    n            = length(x);
    maxmin_range = max(x)-min(x);
    fd_bins1     = ceil(maxmin_range/(2.0*iqr(x)*n^(-1/3))); % Freedman-Diaconis
    
    n            = length(y);
    maxmin_range = max(y)-min(y);
    fd_bins2     = ceil(maxmin_range/(2.0*iqr(y)*n^(-1/3)));
    
    % and use the average...
    fd_bins = ceil((fd_bins1+fd_bins2)/2);
end

%% bin data

edges = linspace(min(x),max(x),fd_bins+1);
[nPerBin1,bins1] = histc(x,edges);

edges = linspace(min(y),max(y),fd_bins+1);
[nPerBin2,bins2] = histc(y,edges);

%% compute entropies

% recompute entropy with optimal bins for comparison
hdat1 = hist(x,fd_bins);
hdat1 = hdat1./sum(hdat1);
hdat2 = hist(y,fd_bins);
hdat2 = hdat2./sum(hdat2);

% convert histograms to probability values
for i=1:2
    eval([ 'entropy(' num2str(i) ') = -sum(hdat' num2str(i) '.*log2(hdat' num2str(i) '+eps));' ]);
end

%% compute joint probabilities

jointprobs = zeros(fd_bins);
for i1=1:fd_bins
    for i2=1:fd_bins
        jointprobs(i1,i2) = sum(bins1==i1 & bins2==i2);
    end
end
jointprobs=jointprobs./sum(jointprobs(:));

entropy(3) = -sum(jointprobs(:).*log2(jointprobs(:)+eps));

%% mutual information

mi = sum(entropy(1:2)) - entropy(3);

%% optional permutation testing

if nargin==4
    
    npermutes = 500;
    n = length(bins2);
    
    perm_mi = zeros(1,npermutes);
    
    for permi=1:npermutes
        
        jointprobs = zeros(fd_bins);
        
        % shuffle bins
        binbreak = randsample(round(n*.8),1,1)+round(n*.1);
        switch mod(permi,4)
            case 0, bins2 = [ bins2(binbreak:end);    bins2(1:binbreak-1) ];
            case 1, bins2 = [ bins2(end:-1:binbreak); bins2(1:binbreak-1) ];
            case 2, bins2 = [ bins2(binbreak:end);    bins2(binbreak-1:-1:1) ];
            case 3, bins2 = [ bins2(end:-1:binbreak); bins2(binbreak-1:-1:1) ];
        end
        
        for i1=1:fd_bins
            for i2=1:fd_bins
                jointprobs(i1,i2) = sum(bins1==i1 & bins2==i2);
            end
        end
        jointprobs=jointprobs./sum(jointprobs(:));
        
        perm_jentropy = -sum(jointprobs(:).*log2(jointprobs(:)+eps));
        
        % mutual information
        perm_mi(permi) = sum(entropy(1:2)) - perm_jentropy;
    end
    
    mi = (mi-mean(perm_mi))/std(perm_mi);
end

%%
% simplified replacement for randsample
function y = randsample(x,n,junk)
y=randperm(x);
y=y(1:n);
