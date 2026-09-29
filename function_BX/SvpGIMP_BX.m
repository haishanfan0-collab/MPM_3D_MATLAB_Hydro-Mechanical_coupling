function [Svp,dSvp] = SvpGIMP_BX(xp,xv,h,lp)

%1D generalised interpolation material point basis functions
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/05/2015
% Description:
% Function to determine the one dimensional GIMP shape functions based on
% the paper: 
% S.G. Bardenhagen, E.M. Kober, The generalized interpolation material 
% point method, Computer Modeling in Eng. & Sciences 5 (2004) 477-496.
%
% The letters on the different cases refer to the regions show in the
% AMPLE paper.  
%--------------------------------------------------------------------------
% [Svp,dSvp] = SVPGIMP(xp,xv,h,lp)
%--------------------------------------------------------------------------
% Input(s):
% xp    - particle position
% xv    - grid node position
% h     - grid spacing 
% lp    - particle half width 
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function 
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------
S1 = ((-h-lp)<(xp-xv)) & ((xp-xv)<=(-h+lp));
S2 = ((-h+lp)<(xp-xv)) & ((xp-xv)<=(  -lp));
S3 = ((  -lp)<(xp-xv)) & ((xp-xv)<=(   lp));
S4 = ((   lp)<(xp-xv)) & ((xp-xv)<=( h-lp));
S5 = (( h-lp)<(xp-xv)) & ((xp-xv)<=( h+lp));

SIZE = max(size(xp),size(xv));
Svp = zeros(SIZE);
dSvp = zeros(SIZE);

if sum(S1,'all')>0
    % Svp = (h+lp+(xp-xv))^2/(4*h*lp);
    % dSvp= (h+lp+(xp-xv))/(2*h*lp);
    LS = (h+lp+(xp-xv)).^2./(4.*h.*lp);
    dLS= (h+lp+(xp-xv))./(2.*h.*lp);
    Svp(S1) = LS(S1);
    dSvp(S1) = dLS(S1);
end

if sum(S2,'all')>0
    % Svp = 1+(xp-xv)/h;
    % dSvp= 1/h;
    LS = 1+(xp-xv)./h;
    dLS = 1./h; dLS = ones(SIZE).*dLS;
    Svp(S2) = LS(S2);
    dSvp(S2) = dLS(S2);
end

if sum(S3,'all')>0
    % Svp = 1-((xp-xv)^2+lp^2)/(2*h*lp);
    % dSvp=-(xp-xv)/(h*lp);
    LS = 1-((xp-xv).^2+lp.^2)./(2.*h.*lp);
    dLS = -(xp-xv)./(h.*lp);
    Svp(S3) = LS(S3);
    dSvp(S3) = dLS(S3);
end

if sum(S4,'all')>0
    % Svp = 1-(xp-xv)/h;
    % dSvp=-1/h;
    LS = 1-(xp-xv)./h;
    dLS = -1./h; dLS = ones(SIZE).*dLS;
    Svp(S4) = LS(S4);
    dSvp(S4) = dLS(S4);
end

if sum(S5,'all')>0
    % Svp = (h+lp-(xp-xv))^2/(4*h*lp);
    % dSvp=-(h+lp-(xp-xv))/(2*h*lp);
    LS = (h+lp-(xp-xv)).^2./(4.*h.*lp);
    dLS = -(h+lp-(xp-xv))./(2.*h.*lp);
    Svp(S5) = LS(S5);
    dSvp(S5) = dLS(S5);
end

end