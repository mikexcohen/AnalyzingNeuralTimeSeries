function Plength = pathlength(A)

% PATHLENGTH Calculate minimum pathlengths for a given adjacency
%           matrix.
%
%   Input   A: n by n adjacency matrix (symmetric).
%
%   Output  Plength: n by n matrix of pathlengths. Element in
%                    position (i,j) is pathlength from node i to node j.
%                    If no path exists, inf is returned.
%
%   Description:    Powers up the adjacency matrix until either there are
%                   no elements equal to zero or the (n-1)st power has been
%                   reached. Records the first power at which (i,j) element
%                   became nonzero.
%
%   Example: Plength = pathlength(A);
%

% This function was downloaded from the following reference:
% Taylor, A., and D.J. Higham. "Contest: A Controllable Test Matrix Toolbox for % Matlab." ACM Transactions on Mathematical Software (TOMS) 35, no. 4 
% (2009): 1-16.
% 

Anew = A;
n = length(A);
power = 1;
Plength = sign(A + eye(n,n)); % record all paths of length one, including diagonal

while any(any(Anew==0)) && power <= (n-1)
    power = power + 1;
    Anew = Anew*A;
    Plength = Plength + ( (Plength == 0) & (Anew > 0) )*power;
end

Plength(Plength==0) = inf;                  % reset zeros to inf

Plength = Plength - diag(diag(Plength));     % reset diagonal to zero