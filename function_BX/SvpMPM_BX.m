function [Svp,dSvp] = SvpMPM_BX(xp,xv,h)

%1D material point basis functions
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   09/02/2016
% Description:
% Function to determine the one dimensional MPM shape functions based on
% global coordinates.
%
%--------------------------------------------------------------------------
% [Svp,dSvp] = SVPMPM(xp,xv,L)
%--------------------------------------------------------------------------
% Input(s):
% xp    - particle position
% xv    - grid node position
% h     - element length
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

S1 = (-h<(xp-xv)) & ((xp-xv)<=0);
S2 = (0<(xp-xv)) & ((xp-xv)<=h);

SIZE = max(size(xp),size(xv));

Svp = zeros(SIZE);
dSvp = zeros(SIZE);

if sum(S1,'all')>0
    % Svp = 1+(xp-xv)/h;
    % dSvp= 1/h;
    LS = 1+(xp-xv)./h;
    dLS = 1./h; dLS = ones(SIZE).*dLS;
    Svp(S1) = LS(S1);
    dSvp(S1) = dLS(S1);
end

if sum(S2,'all')>0
    % Svp = 1-(xp-xv)/h;
    % dSvp=-1/h;
    LS = 1-(xp-xv)./h;
    dLS = -1./h; dLS = ones(SIZE).*dLS;
    Svp(S2) = LS(S2);
    dSvp(S2) = dLS(S2);
end
end