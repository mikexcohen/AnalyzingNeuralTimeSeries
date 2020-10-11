%% Analyzing Neural Time Series Data
% Matlab code for Chapter 10
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% dot products

% two vectors of random numbers
a = randn(10,1);
b = randn(10,1);

% initialize temporary matrix.
pointwise_result=zeros(size(a));

for i=1:length(a)
    pointwise_result(i) = a(i) * b(i);
end
dotproduct = sum(pointwise_result);


% The above code is useful if you are unfamiliar with 
% how a dot product works. Following is a bit more elegant: 
dotproduct = sum(a.*b);

% The most elegant way to compute the dot product, 
% however, is to use Matlab's short-cut: 
dotproduct = a'*b;

% This requires the first vector to be a row vector and the 
% second vector to be a column vector. Otherwise, it will 
% either crash (a*b) or give you the outer product (a*b'). 
% When in doubt, use sum(a.*b). 

%% figure 10.2

% impulse function (all zeros; 1 in the middle)
impfun = zeros(1,100);
impfun(50)=1;
% the figure in the book actually uses the following line, which creates a
% wider boxcar function rather than strictly an impulse function.
impfun(45:55)=1;

kernel = [1 .8 .6 .4 .2]; 

% matlab's convolution function
matlab_conv_result = conv(impfun,kernel,'same');

figure
% plot the signal (impulse or boxcar)
subplot(311)
plot(impfun)
set(gca,'ylim',[-.1 1.1])

% plot the kernel
subplot(312)
plot(kernel,'.-')
set(gca,'xlim',[0 100],'ylim',[-.1 1.1])

% plot the result of convolution
subplot(313)
plot(matlab_conv_result)
set(gca,'xlim',[0 100],'ylim',[-.1 3.6])

%% figure 10.4

% data that we'll use for convolution (must be zero-padded).
dat4conv = [zeros(1,length(kernel)-1) impfun zeros(1,length(kernel)-1) ];

% used for cutting the result of convolution
half_of_kernel_size = ceil((length(kernel)-1)/2);

% initialize convolution output
convolution_result = zeros(1,length(impfun)+length(kernel)-1);

% run convolution (note that kernel is flipped backwards)
for ti=1:length(convolution_result)-half_of_kernel_size
    convolution_result(ti) = sum(dat4conv(ti:ti+length(kernel)-1).*kernel(end:-1:1));
end

% cut off edges
convolution_result = convolution_result(half_of_kernel_size+1:end-half_of_kernel_size);
% Note: In the book figure there was a bug on the previous line ("+1" was typed 
%       as "-1"), which incorrectly showed the convolution as being a few points
%       too far to the right. Thanks to Thijs Perenboom for catching that!

figure
plot(impfun)
hold on
plot(convolution_result,'g')
plot(convolution_result./sum(kernel),'r')
plot(matlab_conv_result./sum(kernel),'ko')
set(gca,'xlim',[0 100],'ylim',[-.1 3.1])
legend({'original timeseries';'unscaled convolution';'manual wavelet convolution';'matlab conv function'})

%% end
