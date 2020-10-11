%% Analyzing Neural Time Series Data
% Matlab code for Chapter 33
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% Figure 33.1

figure

subplot(121)
plot(normpdf(-4:.001:4))
axis tight

subplot(122)
a = randn(1000,1);
hist(a,50)
% this histogram was also used in figure 33.3

disp([ 'p_n = ' num2str(sum(a>2)/1000) ])
disp([ 'p_z = ' num2str(1-normcdf(2))  ])

%% Figure 33.5/6

% These figures are generated in the code for chapter 34

%% Figure 33.8

% introduction to bwlabeln and bwconncomp

% create 2D smoothing kernel
[xi,yi] = meshgrid(-10:10,-10:10);
zi = xi.^2+yi.^2;
zi = 1-(zi/max(zi(:)));

% create a random smoothed map
map = conv2(randn(100),zi,'same');

% threshold map at an arbitrary value
mapt = map;
mapt(abs(map)<range(map(:))/4) = 0;

% get labeled map via bwlabeln (in image processing toolbox)
[mapl,nblobs] = bwlabeln(mapt);
% mapl is a map of size <map> that contains numbers at each cluster. 
% To extract information from clusters:
clustcount = zeros(1,nblobs);
clustsum   = zeros(1,nblobs);
for i=1:nblobs
    clustcount(i) = sum(mapl(:)==i);
    clustsum(i)   = sum(map(mapl(:)==i));
end

% bwconncomp works slightly differently, but will give similar information
blobinfo = bwconncomp(mapt);
% blobinfo is a structure that contains coordinates for each cluster rather than a map
% To extract information from clusters:
clustcount = zeros(1,nblobs);
clustsum   = zeros(1,nblobs);
for i=1:nblobs
    clustcount(i) = numel(blobinfo.PixelIdxList{i});
    clustsum(i)   = sum(map(blobinfo.PixelIdxList{i}));
end
% cluster count can be done faster using cellfun:
clustercount = cellfun(@numel,blobinfo.PixelIdxList);


figure
subplot(131)
imagesc(map), axis square
title('original')

subplot(132)
imagesc(mapt), axis square
title('thresholded')

subplot(133)
imagesc(mapl), axis square
title([ 'labeled (' num2str(nblobs) ' clusters in total)' ])

%% Figure 33.10

% This script calls the function fdr.m, which is included in the online
% code and was downloaded (and slightly modified) in summer 2012 from http://www-personal.umich.edu/~nichols/FDR/FDR.m

nsigs = round(linspace(1,500,80));
nnons = round(linspace(1,500,100));

fdrpvals = zeros(20,length(nsigs),length(nnons));

for iteri=1:20
    for i=1:length(nsigs)
        for j=1:length(nnons)
            
            pvals = [ rand(1,nsigs(i))*.05 rand(1,nnons(j))*.5+.05 ];
            temp = fdr(pvals,.05);
            if isempty(temp)
                fdrpvals(iteri,i,j) = NaN;
            else
                fdrpvals(iteri,i,j) = temp;
            end
        end
    end
end

fdrpvals = squeeze(nanmean(fdrpvals));

figure
imagesc(fdrpvals)
set(gca,'clim',[0 .05])
xlabel('number of non-significant p-values')
ylabel('number of significant p-values')

figure
subplot(211)
plot(nanmean(fdrpvals,1))
xlabel('number of non-significant p-values')
ylabel('critical p-value')
set(gca,'ylim',[0 .05])

subplot(212)
plot(nanmean(fdrpvals,2))
xlabel('number of significant p-values')
ylabel('critical p-value')
set(gca,'ylim',[0 .05])

%% end.
