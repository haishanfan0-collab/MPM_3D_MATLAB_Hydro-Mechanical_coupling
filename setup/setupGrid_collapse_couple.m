function [mpData,mesh] = setupGrid_collapse_couple

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

%% Analysis parameters
global mp_number
E=20e6;   v=0.3;   fc=[5,2];                                              % Young's modulus, Poisson's ratio, yield strength   
g=10;                                                                       % gravity
rho=1790;                                                                   % material density
a = 10;                                                                     % element multiplier
nelsx=4*a;                                                                  % number of elements in the x direction
nelsy=4*a;                                                                  % number of elements in the y direction
ly=8;  lx=8;                                                                % domain dimensions
mp_number=4;                                                                       % number of material points in each direction per element
mpType = 2;                                                                 % material point type: 1 = MPM, 2 = GIMP
cmType = 4;                                                                 % constitutive model: 1 = elastic, 2 = vM plasticity

%% Mesh generation
[etpl,coord] = formCoord2D(4*nelsx,1.25*nelsy,4*lx,1.25*ly);                % background mesh generation
[nels,nen]   = size(etpl);                                                  % number of elements and nodes per element
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly]./[nelsx nelsy];                                      % element lengths in each direction

%% Boundary conditions on backgroun mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,1)==0 || coord(node,1)==4*lx                                % roller (x=0)
    bc(node*2-1,:)=[node*2-1 0];    
  end
  if coord(node,2)==0                                                       % roller (y=0)
    bc(node*2  ,:)=[node*2   0];
  end
end
bc = bc(bc(:,1)>0,:);                                                       % remove empty part of bc                                     % remove empty part of bc

%% Boundary Flow conditions on backgroun mesh
bc_Flow = zeros(nodes,2);                                                   % generate empty bc matrix
% for node=1:nodes                                                          % loop over nodes
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

%% Material point generation
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
JL=[];
FHS_AHA3 = max(mpC(:,2));
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
  mpData(mp).Fn     = eye(3);                                               % previous deformation gradient
  mpData(mp).F      = eye(3);                                               % deformation gradient
  mpData(mp).sig    = zeros(6,1);                                           % Cauchy stress
  mpData(mp).sign    = zeros(6,1);                                          % Cauchy stress
  mpData(mp).epsEn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsE   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).epsPn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsP   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).mCst   = [E*2 v fc*2];
   if (mpC(mp,1)+mpC(mp,2)/1.5)>(lx+ly)*0.5*0.85
       mpData(mp).mCst   = [E v fc];
   end
  mpData(mp).fp     = zeros(nD,1);                                          % point forces at material points
  mpData(mp).u      = zeros(nD,1);                                          % material point displacements
  mpData(mp).du     = zeros(nD,1);
  mpData(mp).ddu    = zeros(nD,1);

  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                           % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                           % initial material point domain lengths (GIMP)
  else
    mpData(mp).lp     = zeros(1,nD);                                        % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                        % initial material point domain lengths (MPM)
  end

  %% 非饱和渗流参数
  mpData(mp).Flow_P = 0.0;
  mpData(mp).Flow_Type = 1;                                                 % 基质吸力模型，1-Van Genuchten模型
  mpData(mp).dH     = 0.0;
  mpData(mp).Flow_CP = [0.016, 0.500, 2.00];                                % 基质吸力模型参数
  
  FHS_c = (mpC(mp,2))/FHS_AHA3;
  FHS_c1 = 4*FHS_c^2-8*FHS_c+4;
  mpData(mp).Gravity_Water_Content = 0.17+FHS_c1/200;                        % 初始质量含水率 0.14~0.18

  mpData(mp).Gs = 2770;                                                     % 比重
  mpData(mp).Roumax = 1850;                                                 % 最大干密度
  rhoA = rho*(1+mpData(mp).Gravity_Water_Content);
  mpData(mp).mpM = vp(mp)*rhoA;                                             % 重新获取质量
  mpData(mp).e0 = mpData(mp).Gs/rho-1;                                      % 初始空隙率
  mpData(mp).e = mpData(mp).e0;                                             % 孔隙率
  mpData(mp).Sr = mpData(mp).Gravity_Water_Content*mpData(mp).Gs/1000/mpData(mp).e; % 饱和度计算
  mpData(mp).Ksat = [3e-7;1e-7];                 % 饱和渗透系数
  if (mpC(mp,1)+mpC(mp,2)/1.5)>=(lx+ly)*0.5*0.85-0.20 && (mpC(mp,1)+mpC(mp,2)/1.5)<=(lx+ly)*0.5*0.85+0.20
       mpData(mp).Ksat = [5e-5;2e-5];
  end
  mpData(mp).Ksat0 = mpData(mp).Ksat;                                       % 饱和渗透系数
  mpData(mp).eps_Flow = zeros(2,1);                                         % 饱和渗透系数
  mpData(mp).sig_Flow = zeros(2,1);                                         % 饱和渗透系数
  mpData(mp).epsn_Flow = zeros(2,1);                                        % 饱和渗透系数
  mpData(mp).sign_Flow = zeros(2,1);                                        % 饱和渗透系数
  Fa = mpData(mp).Flow_CP(1);
  Fm = mpData(mp).Flow_CP(2);
  Fn = mpData(mp).Flow_CP(3);
  mpData(mp).Scution_Flow = ((mpData(mp).Sr^(-1/Fm) - 1)^(1/Fn))/Fa;        % 饱和渗透系数
  mpData(mp).H0      = mpData(mp).Scution_Flow/-9.81;
  mpData(mp).dH_sum = 0.0;
  mpData(mp).dH_sum_strain = 0.0;
  mpData(mp).Roud   = rho;
  mpData(mp).Saturated_YES   = 0;
  %% 绘图
  % Sr = 0.5:0.01:1.0;
  % Scution = ((Sr.^(-1/Fm) - 1).^(1/Fn))./Fa;
  % plot(Sr,Scution)

  % 删除
  if (mpC(mp,1)+mpC(mp,2)/2)>lx
      JL(end+1)=mp;
  end

end
mpData(JL)=[];

% boundary_mp_generate
mpData = boundary_mp_generate(mpData,nen,[E v fc],FHS_AHA3,rho);

% Flow_boundary
% boundary_mp_generate
mpData = boundary_Flow(mpData,{'Slope','Up'},2*10^-5,'Right');
end