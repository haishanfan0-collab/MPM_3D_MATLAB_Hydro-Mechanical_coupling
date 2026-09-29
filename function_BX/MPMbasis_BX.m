function [Svp,dSvp] = MPMbasis_BX(mesh,mpData,node)

%Basis functions for the material point method
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to determine the multi-dimensional MPM shape functions from the
% one dimensional MPM functions.  The function includes both the standard
% and generalised interpolation material point methods. 
%
%--------------------------------------------------------------------------
% [Svp,dSvp] = MPMBASIS(coord,mpC,L)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh data structured array. Function requires:
%           - coord  : nodal coordinates  
%           - h      : grid spacing
%
% mpData - material point structured array.  Function requires:
%           - mpC    : material point coordinates (single point)
%           - lp     : particle domain lengths
%           - mpType : material point type (1 or 2)
%
% node   - background mesh node number
%--------------------------------------------------------------------------
% Ouput(s);
% Svp   - particle characteristic function
% dSvp  - gradient of the characterstic function 
%--------------------------------------------------------------------------
% See also:
% SVPMPM    - MPM basis functions in 1D (mpType = 1
% SVPGIMP   - GIMPM basis functions in 1D (mpType = 2)
%--------------------------------------------------------------------------

coord  = mesh.coord(node,:);                                                % node coordinates
h      = mesh.h;                                                            % grid spacing
mpL = size(mpData,2);
nD     = length(mpData(1).mpC);                                                     % number of dimensions

if mpL==1
    mpC    = mpData.mpC;                                                    % material point coordinates
    lp     = mpData.lp;                                                     % material point domain length
    mpType = mpData.mpType;                                                 % material point type (MPM or GIMPM)
    
    dSvp=zeros(size(coord));                                                    % zero vectors used in calcs

    if mpType == 1
        [S,dS] = SvpMPM_BX(mpC,coord,h);                                        % 1D MPM functions
    elseif mpType == 2 || mpType == 3
        [S,dS] = SvpGIMP_BX(mpC,coord,h,lp);                                    % 1D GIMPM functions
    end

else
    mpC    = reshape([mpData.mpC]',2,mpL)';                                     % material point coordinates
    lp     = reshape([mpData.lp]',2,mpL)';                                      % material point domain length
    mpType = [mpData.mpType]';                                                  % material point type (MPM or GIMPM)

    mpType_GI_YES = (mpType==1);

    S = zeros(length(mpC),nD);
    dS = zeros(length(mpC),nD);

    dSvp=zeros(size(mpC));                                                      % zero vectors used in calcs

    if any(mpType_GI_YES)
        [S0,dS0] = SvpMPM_BX(mpC,coord,h);                                        % 1D MPM functions
        S(mpType_GI_YES,:)=S0(mpType_GI_YES,:);
        dS(mpType_GI_YES,:)=dS0(mpType_GI_YES,:);
    end
    [S0,dS0] = SvpGIMP_BX(mpC,coord,h,lp);                                    % 1D GIMPM functions
    S(~mpType_GI_YES,:)=S0(~mpType_GI_YES,:);
    dS(~mpType_GI_YES,:)=dS0(~mpType_GI_YES,:);
end

if nD == 1
    indx = [];                                                              % index for basis derivatives (1D)
elseif nD == 2
    indx = [2; 1];                                                          % index for basis derivatives (2D)
elseif nD == 3
    indx = [2 3; 1 3; 1 2];                                                 % index for basis derivatives (3D)
end

Svp=prod(S,2);                                                                % basis function
for i=1:nD                                                                  
    dSvp(:,i)=dS(:,i).*prod(S(:,indx(i,:)),2);                                       % gradient of the basis function
end

Svp = Svp';
dSvp = dSvp';

end