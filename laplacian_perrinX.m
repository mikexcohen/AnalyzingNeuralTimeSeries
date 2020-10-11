%LAPLACIAN_PERRINX   Compute surface Laplacian of EEG data.
% [surf_lap,G,H] = laplacian_perrinX(data,x,y,z[,leg_order,smoothing]);
% 
% INPUTS    : 
%      data : EEG data (can be N-D, but first dimension must be electrodes)
%     x,y,z : x,y,z coordinates of electrode positions (e.g., [EEG.chanlocs.X])
%
% (optional inputs)
% leg_order : order of Legendre polynomial (default is 20 [40 for >100 electrodes])
% smoothing : G smoothing parameter (lambda), set to 1e-5 by default
%
%
% OUTPUTS   :
%  surf_lap : the surface Laplacian (second spatial derivative)
% (optional outputs)
%       G,H : G and H matrices
% 
%  This is an implementation of algorithms described by 
%    Perrin, Pernier, Bertrand, and Echallier (1989). PubMed #2464490

% mikexcohen@gmail.com

function [surf_lap,G,H] = laplacian_perrinX(data,x,y,z,varargin) % vararg order: leg_order,smoothing

numelectrodes = numel(x);

if nargin<4
    help laplacian_perrinX
    error('Read help file!')
end

%% compute G and H matrices

% initialize
G=zeros(numelectrodes);
H=zeros(numelectrodes);
cosdist=zeros(numelectrodes);

% default parameters for +/- 100 electrodes
if numelectrodes>100
    m=3; leg_order=40;
else
    m=4; leg_order=20;
end

if numel(varargin)>0 && ~isempty(varargin{1})
    leg_order=varargin{1};
end

% scale XYZ coordinates to unit sphere
[junk,junk,spherical_radii] = cart2sph(x,y,z);
maxrad = max(spherical_radii);
x = x./maxrad;
y = y./maxrad;
z = z./maxrad;

for i=1:numelectrodes
    for j=i+1:numelectrodes
        cosdist(i,j) = 1 - (( (x(i)-x(j))^2 + (y(i)-y(j))^2 + (z(i)-z(j))^2 ) / 2 );
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
hdenom = ((1:leg_order).*((1:leg_order)+1)).^(m-1);

for i=1:numelectrodes
    for j=i:numelectrodes
        
        g=0; h=0;
        
        for ni=1:leg_order
            % compute G and H terms
            g = g + (twoN1(ni)*legpoly(ni,i,j)) / gdenom(ni);
            h = h - (twoN1(ni)*legpoly(ni,i,j)) / hdenom(ni);
        end
        G(i,j) =  g/(4*pi);
        H(i,j) = -h/(4*pi);
    end
end

% mirror matrix
G=G+G'; H=H+H';

% correct for diagonal-double
G = G-eye(numelectrodes)*G(1)/2;
H = H-eye(numelectrodes)*H(1)/2;

%% compute laplacian

% reshape data to electrodes X time/trials
orig_data_size = squeeze(size(data));
if any(orig_data_size==1)
    data=data(:);
else
    data = reshape(data,orig_data_size(1),prod(orig_data_size(2:end)));
end

% smoothing constant
if numel(varargin)==2
    smoothing=varargin{2};
else
    smoothing=1e-5;
end

% add smoothing constant to diagonal 
% (change G so output is unadulterated)
Gs = G + eye(numelectrodes)*smoothing;

% compute C matrix
GsinvS = sum(inv(Gs));
dataGs = data'/Gs;
C      = dataGs - (sum(dataGs,2)/sum(GsinvS))*GsinvS;

% compute surface Laplacian (and reshape to original data size)
surf_lap = reshape((C*H')',orig_data_size);

%% end
