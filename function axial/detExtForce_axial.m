function [fext] = detExtForce_axial(nodes,nD,g,mpData,g_wd,Contact,BS)

%Global external force determination  
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the external forces at nodes based on body forces
% and point forces at material points.
%
%--------------------------------------------------------------------------
% [fbdy,mpData] = DETEXTFORCE(coord,etpl,g,eIN,mpData)
%--------------------------------------------------------------------------
% Input(s):
% nodes  - number of nodes (total in mesh)
% nD     - number of dimensions
% g      - gravity
% mpData - material point structured array. Function requires:
%           mpM   : material point mass
%           nIN   : nodes linked to the material point
%           Svp   : basis functions for the material point
%           fp    : point forces at material points
%--------------------------------------------------------------------------
% Ouput(s);
% fext   - external force vector (nodes*nD,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------
if nargin<5
   g_wd = 'z';
   Contact = zeros(size(mpData,2),nD);
   BS=1.0;
end

if nargin<6
   Contact = zeros(size(mpData,2),nD);
   BS=1.0;
end

if nargin<7
   BS=1.0;
end

if g_wd == 'x'
    g_W = 1;
elseif g_wd == 'y'
    g_W = 2;
elseif g_wd == 'z'
    g_W = 3;
else
    error('g_wd only can be ''x'', ''y'', ''z''!')
end

nmp  = size(mpData,2);                                                      % number of material points & dimensions 
fext = zeros(nodes*nD,1);                                                   % zero the external force vector
grav = zeros(nD,1); grav(g_W) = -g;                                         % gavity vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                                    % nodes associated with MP
   nn  = length(nIN);                                                       % number of nodes influencing the MP
   Svp = mpData(mp).Svp;                                                    % basis functions
   fp  = 2*pi*mpData(mp).mpC(1)*(mpData(mp).mpM*grav + BS * mpData(mp).fp + Contact(mp,:)')*Svp;         % material point body & point nodal forces
   ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                     % nodel degrees of freedom 
   fext(ed) = fext(ed) + fp;                                                % combine into external force vector
end
end