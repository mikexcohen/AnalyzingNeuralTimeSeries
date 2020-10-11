%INTERPOLATE_PERRINX   Interpolate EEG electrodes.
% interp_data = interpolate_perrinX(data,x,y,z[,leg_order,smoothing]);
% 
% INPUTS     : 
%       data : EEG data (can be N-D, but first dimension must be electrodes)
%              Input all data, not only good electrodes! 
%      x,y,z : x,y,z coordinates of electrode positions (e.g., [EEG.chanlocs.X])
%  bad_elecs : index of electrodes to interpolate
%
%  (optional inputs)
%  leg_order : order of Legendre polynomial (default is 10)
%  smoothing : G smoothing parameter (lambda), set to 1e-5 by default
%
%
% OUTPUTS    :
%       data : entire input data with bad electrodes replaced by
%              interpolated data
% (optional outputs)
%       G,Gi : G matrices for good electrodes and interpolated electrodes
% 
%  This is an implementation of algorithms described by 
%    Perrin, Pernier, Bertrand, and Echallier (1989). PubMed #2464490 

% mikexcohen@gmail.com

function [data,G,Gi] = interpolate_perrinX(data,x,y,z,bad_electrodes,varargin) % vararg order: leg_order,smoothing

if nargin<5
    help interpolate_perrinX
    error('Read help file!')
end

%% separate the goodies from the baddies

good_electrodes = true(size(x));
good_electrodes(bad_electrodes) = false;

gx = x(good_electrodes);
gy = y(good_electrodes);
gz = z(good_electrodes);
bx = x(bad_electrodes);
by = y(bad_electrodes);
bz = z(bad_electrodes);

numelectrodes  = numel(gx);
numelectrodesi = numel(bx);

m=2;
if numel(varargin)>0 && ~isempty(varargin{1})
    leg_order=varargin{1};
else
    leg_order=30;
end

% scale XYZ coordinates to unit sphere
[junk,junk,spherical_radii] = cart2sph(x,y,z);
maxrad = max(spherical_radii);
gx = gx./maxrad;
gy = gy./maxrad;
gz = gz./maxrad;
bx = bx./maxrad;
by = by./maxrad;
bz = bz./maxrad;

%% compute G matrix for good electrodes

% initialize
G=zeros(numelectrodes);
cosdist=zeros(numelectrodes);

for i=1:numelectrodes
    for j=i+1:numelectrodes
        cosdist(i,j) = 1 - (( (gx(i)-gx(j))^2 + (gy(i)-gy(j))^2 + (gz(i)-gz(j))^2 ) / 2 );
    end
end
cosdist = cosdist+cosdist' + eye(numelectrodes);


% compute Legendre polynomial
legpoly = zeros(leg_order,numelectrodes,numelectrodes);
for ni=1:leg_order
    temp = legendre(ni,cosdist);
    legpoly(ni,:,:) = temp(1,:,:);
end

% precompute electrode-independent variables
twoN1  = 2*(1:leg_order)+1;
gdenom = ((1:leg_order).*((1:leg_order)+1)).^m;

for i=1:numelectrodes
    for j=i:numelectrodes
        
        g=0;
        for ni=1:leg_order
            % compute G and H terms
            g = g + (twoN1(ni) * legpoly(ni,i,j)) / gdenom(ni);
        end
        G(i,j) =  g/(4*pi);
    end
end
% symmetric matrix
G=G+G';
G = G-eye(numelectrodes)*G(1)/2;

% add smoothing constant to G diagonal
if numel(varargin)==2
    smoothing=varargin{2};
else
    smoothing=1e-5;
end

% add smoothing constant to diagonal 
% (change G so output is unadulterated)
Gs = G + eye(numelectrodes)*smoothing;

%% compute G matrix for to-be-interpolated electrodes


Gi=zeros(numelectrodesi,numelectrodes);
cosdisti=zeros(numelectrodesi);

for i=1:numelectrodesi
    for j=1:numelectrodes
        cosdisti(i,j) = 1 - (( (bx(i)-gx(j))^2 + (by(i)-gy(j))^2 + (bz(i)-gz(j))^2 ) / 2 );
    end
end

% compute Legendre polynomial
legpolyi = zeros(leg_order,numelectrodesi,numelectrodes);
for ni=1:leg_order
    temp = legendre(ni,cosdisti);
    legpolyi(ni,:,:) = temp(1,:,:);
end

for i=1:numelectrodesi
    for j=1:numelectrodes
        
        g=0;
        for ni=1:leg_order
            % compute G and H terms
            g = g + (twoN1(ni) * legpolyi(ni,i,j)) / gdenom(ni);
        end
        Gi(i,j) =  g/(4*pi);
    end
end

%% interpolate

% reshape data to electrodes X time/trials
orig_data_size = squeeze(size(data));
if any(orig_data_size==1)
    data=data(:);
else
    data = reshape(data,orig_data_size(1),prod(orig_data_size(2:end)));
end

% interpolate and reshape to original data size
data(bad_electrodes,:) = Gi*(pinv(G)*data(good_electrodes,:));
data = reshape(data,orig_data_size);

%% end
