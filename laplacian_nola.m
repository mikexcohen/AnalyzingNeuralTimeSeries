function surf_lap = laplacian_nola(x,y,z,data,smoothing)
%LAPLACIAN_NOLA   Compute surface Laplacian of EEG data via New Orleans method
% surf_lap = laplacian_nola(x,y,z,data[,smoothing]);
% 
% INPUTS    : 
%     x,y,z : x,y,z coordinates of electrode positions
%      data : EEG data (can be N-D, but first dimension must be electrodes)
% (optional inputs)
% smoothing : smoothing parameter
%
%
% OUTPUTS   :
%  surf_lap : the surface Laplacian (second spatial derivative)
% (optional outputs)
% 
% NOTES
%
%   (1) This script is rewritten from appedix J of the Nunez and Srinivasan book 
%       (2nd edition) by Mike X Cohen (algorithms are unchanged)


if nargin<5
    smoothing = 100;
end

n = numel(x);

% budge zero values
x(x==0) = .001;
y(y==0) = .001;
z(z==0) = .001;

%% compute K

k = zeros(n);

for i=1:n
    for j=i+1:n
        s = x(i) - x(j);
        t = y(i) - y(j);
        r = z(i) - z(j);
        str = s.^2 + t.^2 + r.^2;
        k(i,j) = ((str+smoothing).^2) * log(str+smoothing);
    end
end

k = k+k';
kinv = pinv(k);

%% compute E and A

e = zeros(n,10);

e(:,1)  = 1;
e(:,2)  = x;
e(:,3)  = y;
e(:,4)  = x.^2;
e(:,5)  = x.*y;
e(:,6)  = y.^2;
e(:,7)  = z;
e(:,8)  = z.*x;
e(:,9)  = z.*y;
e(:,10) = z.^2;

ke = kinv*e;

et = e';

a = et*ke;

ainv = pinv(a);

%% compute laplacian over data

orig_data_size = squeeze(size(data));
if any(orig_data_size==1)
    data=data(:);
end

data = reshape(data,orig_data_size(1),prod(orig_data_size(2:end)));
surf_lap = zeros(size(data));

for ti=1:size(data,2)
    
    kv = kinv*data(:,ti);
    ev = et*kv;
    
    q = ainv*ev;
    
    eq = e*q;
    keq = kinv*eq;
    
    % compute p
    p = kv - keq;
    
    %% compute laplacian
    
    [az,el,r] = cart2sph(x,y,z);
    el        = pi*ones(size(el))/2 - el;
    
    % trig functions
    st = sin(el);
    ct = cos(el);
    sp = sin(az);
    cp = cos(az);
    
    
    uuxyz = 2*q(4)+2*q(6)+2*q(10) - (2*st.*(q(2)*cp+q(3)*sp)./r +2*q(7)*ct./r+6*(st.^2).*(q(4)*(cp.^2)+q(6)*(sp.^2)+q(5)*sp.*cp)+6*st.*ct.*(q(8)*cp+q(9)*sp)+6*q(10)*ct.^2);
    
    ttcomp = zeros(size(st));
    rrcomp = zeros(size(st));
    
    for j = 1:n
        a = r(j)*(st.*cp-sin(el(j))*cos(az(j))*ones(size(st,1),size(st,2)));
        b = r(j)*(st.*sp-sin(el(j))*sin(az(j))*ones(size(st,1),size(st,2)));
        c = r(j)*(ct-cos(el(j))*ones(size(st,1),size(st,2)));
        
        str  = a.^2+b.^2+c.^2;
        strw = str+smoothing*smoothing*ones(size(st,1),size(st,2));
        
        comterm  = 4*str./strw-((str./strw).^2)+2*log(strw);
        comterm2 = 2*(2*str.*log(strw)+(str.^2)./strw);
        
        tcomp = 3*comterm2+4*str.*comterm;
        dr    = 2*(a.*st.*cp+b.*st.*sp+c.*ct);
        
        rcomp  = dr.*comterm2+2*r(j)*comterm2/2+r(j)*(dr.^2).*comterm;
        ttcomp = ttcomp + p(j)*tcomp;
        rrcomp = rrcomp + p(j)*rcomp/r(j);
    end;
    surf_lap(:,ti) = -(ttcomp+uuxyz-rrcomp);
    
end

surf_lap = reshape(surf_lap,orig_data_size);

%% end.

