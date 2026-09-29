function [mpData,mesh] = pavement_structure

%Problem setup information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Problem setup for a large deformation elasto-plastic collapse.
%
%--------------------------------------------------------------------------
% [lstps,g,mpData,mesh] = SETUPGRID2D
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
% lstps  - number of loadsteps (1)
% g      - gravity (1)
% mpData - structured array with the following fields:
%           - mpType : material point type (1 = MPM, 2 = GIMPM)
%           - cmType : constitutive model type (1 = elastic, 2 = vM)
%           - mpC    : material point coordinates
%           - vp     : material point volume
%           - vp0    : initial material point volume
%           - mpM    : material point mass
%           - nIN    : nodes linked to the material point
%           - eIN    : element associated with the material point
%           - Svp    : basis functions for the material point
%           - dSvp   : basis function derivatives (at start of lstp)
%           - Fn     : previous deformation gradient
%           - F      : deformation gradient
%           - sig    : Cauchy stress
%           - epsEn  : previous logarithmic elastic strain
%           - epsE   : logarithmic elastic strain
%           - mCst   : material constants (or internal state parameters)
%           - fp     : force at the material point
%           - u      : material point displacement
%           - lp     : material point domain lengths
%           - lp0    : initial material point domain lengths
%
% mesh   - structured array with the following fields
%           - coord : mesh nodal coordinates (nodes,nD)
%           - etpl  : element topology (nels,nen)
%           - bc    : boundary conditions (*,2)
%           - h     : background mesh size (nD,1)
%--------------------------------------------------------------------------
% See also:
% FORMCOORD2D - background mesh generation
% DETMPPOS    - local material point positions
% SHAPEFUNC   - background grid basis functions
%--------------------------------------------------------------------------
global mp_number
%% Analysis parameters
mcSM_93 = [2.04	0.18	1.94...
  	0.18	0.002	-2.1	1.48	1.0	1.0	5.8	3.56	...
    0.081	2.0	1730.0	2660.0	0.93	0.35	1.42	...
    0.44	0.88	-1.27	3.37	-0.06	-4.0	0.67 ...	
    1.5	128000.0	80000.0	2e-05	0.2047	0.2247	0.01];  
mcSM_94 = [2.04	0.18	1.94...
  	0.18	0.002	-2.1	1.48	1.0	1.0	5.8	3.56	...
    0.081	2.0	1730.0	2660.0	0.94	0.35	1.42	...
    0.44	0.88	-1.27	3.37	-0.06	-4.0	0.67 ...	
    1.5	128000.0	80000.0	2e-05	0.1847	0.2047	0.01]; 
mcSM_96 = [2.04	0.18	1.94...
  	0.18	0.002	-2.1	1.48	1.0	1.0	5.8	3.56	...
    0.081	2.0	1730.0	2660.0	0.96	0.35	1.42	...
    0.44	0.88	-1.27	3.37	-0.06	-4.0	0.67 ...	
    1.5	128000.0	80000.0	2e-05	0.1647	0.1847	0.01]; 
g=9.81;                                                                     % gravity
rho=1790;                                                                   % material density
a = 10;                                                                     % element multiplier
nelsx=12*a;                                                                  % number of elements in the x direction
nelsy=4*a;                                                                  % number of elements in the y direction
ly=8;  lx=24;                                                                % domain dimensions
mp_number=3;                                                                       % number of material points in each direction per element
mpType = 2;                                                                 % material point type: 1 = MPM, 2 = GIMP
cmType = 5;                                                                 % constitutive model: 1 = elastic, 2 = vM plasticity

%% Mesh generation
[etpl,coord] = formCoord2D(2*nelsx,1.25*nelsy,2*lx,1.25*ly);                % background mesh generation
[nels,nen]   = size(etpl);                                                  % number of elements and nodes per element
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly]./[nelsx nelsy];                                      % element lengths in each direction

%% Boundary conditions on backgroun mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,1)==0                                                       % roller (x=0)
    bc(node*2-1,:)=[node*2-1 0];    
  end
  if coord(node,2)==0                                                       % roller (y=0)
    bc(node*2  ,:)=[node*2   0];
  end
end
bc = bc(bc(:,1)>0,:);                                                       % remove empty part of bc                                     % remove empty part of bc

%% Boundary Flow conditions on backgroun mesh
bc_Flow = zeros(nodes,2);                                                   % generate empty bc matrix
% for node=1:nodes                                                            % loop over nodes
%   if (coord(node,1) == 0.0 && coord(node,2) == 0.0)     % roller (y=0)
%     bc_Flow(node  ,:)=[node   0];
%   end
% end
bc_Flow = bc_Flow(bc_Flow(:,1)>0,:);                                        % remove empty part of bc

%% Element limits for MP-element searching
eMin = zeros(nels,nD);                                                      % element lower coordinate limit
eMax = zeros(nels,nD);                                                      % element upper coordainte limit 
for i = 1:nD
    ci = coord(:,i);                                                        % nodal coordinates in current i direction
    c  = ci(etpl);                                                          % reshaped element coordinates in current i direction
    eMin(:,i) = min(c,[],2);                                                % element lower coordinate limit 
    eMax(:,i) = max(c,[],2);                                                % element upper coordainte limit 
end

%% Mesh data structure generation
mesh.etpl    = etpl;                                                        % element topology
mesh.coord   = coord;                                                       % nodal coordinates
mesh.bc      = bc;                                                          % boundary conditions
mesh.bc      = bc;                                                          % boundary conditions
mesh.h       = h;                                                           % mesh size
mesh.eMin    = eMin;                                                        % element lower coordinate limit 
mesh.eMax    = eMax;                                                        % element upper coordainte limit 
mesh.g = g;
mesh.bc_Flow = bc_Flow;

ngp    = mp_number^nD;                                                             % number of material points per element
GpLoc  = detMpPos(mp_number,nD);                                                   % local MP locations
N      = shapefunc(nen,GpLoc,nD);                                           % basis functions for the material points
[etplmp,coordmp] = formCoord2D(nelsx,nelsy,lx,ly);                          % mesh for MP generation
nelsmp = size(etplmp,1);                                                    % no. elements populated with material points
nmp    = ngp*nelsmp;                                                        % total number of mterial points

mpC=zeros(nmp,nD);                                                          % zero MP coordinates
for nel=1:nelsmp
  indx=(nel-1)*ngp+1:nel*ngp;                                               % MP locations within mpC
  eC=coordmp(etplmp(nel,:),:);                                              % element coordinates
  mpPos=N*eC;                                                               % global MP coordinates
  mpC(indx,:)=mpPos;                                                        % store MP positions
end
lp = zeros(nmp,2);                                                          % zero domain lengths
lp(:,1) = h(1)/(2*mp_number);                                                      % domain half length x-direction
lp(:,2) = h(2)/(2*mp_number);                                                      % domain half length y-direction
vp      = 2^nD*lp(:,1).*lp(:,2);                                            % volume associated with each material point

%% Material point structure generation
for mp = nmp:-1:1                                                           % loop backwards over MPs so array doesn't change size
  mpData(mp).mpType = mpType;                                               % material point type: 1 = MPM, 2 = GIMP
  mpData(mp).cmType = cmType;                                               % constitutive model: 1 = elastic, 2 = vM plasticity
  mpData(mp).mpC    = mpC(mp,:);                                            % material point coordinates
  mpData(mp).vp     = vp(mp);                                               % material point volume
  mpData(mp).vp0    = vp(mp);                                               % material point initial volume
  mpData(mp).mpM    = vp(mp)*rho;                                           % material point mass
  mpData(mp).nIN    = zeros(nen,1);                                         % nodes associated with the material point
  mpData(mp).eIN    = 0;                                                    % element associated with the material point
  mpData(mp).Svp    = zeros(1,nen);                                         % material point basis functions
  mpData(mp).dSvp   = zeros(nD,nen);                                        % derivative of the basis functions
  mpData(mp).Gvp  = zeros(nD,nen);                                                 % basis function derivatives
  mpData(mp).Tvp  = zeros(1,nen);
  mpData(mp).Fn     = eye(3);                                               % previous deformation gradient
  mpData(mp).F      = eye(3);                                               % deformation gradient
  mpData(mp).sig    = zeros(6,1);                                           % Cauchy stress
  mpData(mp).sign    = zeros(6,1);                                          % Cauchy stress
  mpData(mp).epsEn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsE   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).epsPn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsP   = zeros(6,1);                                           % elastic strain (logarithmic)
  if mpC(mp,2)>=ly - 0.80 -0.78
      mpData(mp).mCst   = mcSM_96;
  elseif mpC(mp,2)>=ly - 1.50 -0.78
      mpData(mp).mCst   = mcSM_94;
  else
      mpData(mp).mCst   = mcSM_93;
  end
  mpData(mp).fp     = zeros(nD,1);                                          % point forces at material points
  mpData(mp).u      = zeros(nD,1);                                          % material point displacements
  mpData(mp).du     = zeros(nD,1);
  mpData(mp).ddu    = zeros(nD,1);
  mpData(mp).STATEV    = zeros(19,1);
  mpData(mp).STATEVn    = zeros(19,1);
  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                           % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                           % initial material point domain lengths (GIMP)
  else
    mpData(mp).lp     = zeros(1,nD);                                        % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                        % initial material point domain lengths (MPM)
  end
  mpData(mp).Scution_Flow = 1.0;
  mpData(mp).Gravity_Water_Content = 0.0;
  mpData(mp).dH_sum = 0.0;
  mpData(mp).dH_sum_strain=0.0;
  mpData(mp).modulus=0.0;
end

% boundary_mp_generate
mpData = boundary_mp_generate(mpData,nen,mcSM_96,1,rho);
LS=[];
nmp  = length(mpData);
for mp = nmp:-1:1
    if (mpData(mp).mpC(1)<=1e-3 || mpData(mp).mpC(1)>=lx || mpData(mp).mpC(2)<=1e-3) || ...
       (mpData(mp).mpC(1)*2/3+mpData(mp).mpC(2)>48/3)
        LS(end+1) = mp;                                         % point forces at material points
    end

    if mpData(mp).mpC(1)>=1.00 && mpData(mp).mpC(1)<=12-0.75 && mpData(mp).mpC(2)>= ly - 0.78
        mpData(mp).cmType = 1;
        mpData(mp).mCst = [1800e6,0.25];
        % mpData(mp).modulus = mpData(mp).mCst(1)/1e6;
    end

end
mpData(LS)=[];

l11 = 3.625-0.9-0.15;
l12 = 3.625-0.9+0.15;
l21 = 3.625+0.9-0.15;
l22 = 3.625+0.9+0.15;

l31 = 7.375-0.9-0.15;
l32 = 7.375-0.9+0.15;
l41 = 7.375+0.9-0.15;
l42 = 7.375+0.9+0.15;

nmp  = length(mpData);
mpC = reshape([mpData.mpC],nD,nmp)';                                        % all material point coordinates (nmp,nD)

A11 = mpC(:,1)>=l11;
A12 = mpC(:,1)<=l12;

A21 = mpC(:,1)>=l21;
A22 = mpC(:,1)<=l22;

A31 = mpC(:,1)>=l31;
A32 = mpC(:,1)<=l32;

A41 = mpC(:,1)>=l41;
A42 = mpC(:,1)<=l42;

A2 = mpC(:,2)>=ly;

A1 = A11.*A12+A21.*A22+A31.*A32+A41.*A42;
A1 = A1>0;

A = sum(A2.*A1);
F = 200000/A;

for mp = nmp:-1:1
    if ((mpC(mp,1)>=l11 && mpC(mp,1)<=l12) || ...
        (mpC(mp,1)>=l21 && mpC(mp,1)<=l22) || ...
        (mpC(mp,1)>=l31 && mpC(mp,1)<=l32) || ...
        (mpC(mp,1)>=l41 && mpC(mp,1)<=l42)) &&...
        mpC(mp,2)>=ly
        mpData(mp).fp     = [0;-F];
    else
        mpData(mp).fp     = zeros(nD,1);                                    % point forces at material points
    end
end

end