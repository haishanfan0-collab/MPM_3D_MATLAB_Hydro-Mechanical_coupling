function [lstps,g,mpData,mesh] = setupGrid3D

% Problem setup information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Problem setup for the compaction of a 2D column under self weight with
% roller boundary conditions applied to the sides and the base of the
% domain.
%
%--------------------------------------------------------------------------
% [lstps,g,mpData,mesh] = SETUPGRID
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
% lstps  - number of loadsteps (1) / 迭代上限
% g      - gravity (1)             / 重力加速度
% mpData - structured array with the following fields:
%           - mpType : material point type (1 = MPM, 2 = GIMPM)             /物质点类别 （1-标准MPM，2-GIMPM<物质点有限体积方法>）
%           - cmType : constitutive model type (1 = elastic, 2 = vM plas.)  /本构模型类型
%           - mpC    : material point coordinates                           /物质点的坐标
%           - vp     : material point volume                                /物质点的体积（后续自动求解）
%           - vp0    : initial material point volume                        /vp(t=0)
%           - mpM    : material point mass                                  /物质点质量
%           - nIN    : nodes linked to the material point                   /物质点所属单元对应的节点序列（自动识别）
%           - eIN    : element associated with the material point           /物质点所数单元（自动识别）
%           - nSMe   : number of stiffness matrix entries for the MP        /刚度矩阵维度（自动识别）
%           - Svp    : basis functions for the material point               /材料点的基函数（自动求解）
%           - dSvp   : basis function derivatives (at start of lstp)        /基函数导数求解（自动求解）
%           - Fn     : previous deformation gradient                        /上一步的变形梯度（自动求解）
%           - F      : deformation gradient                                 /当前变形梯度（自动求解）
%           - sig    : Cauchy stress                                        /柯西应力（自动求解）
%           - epsEn  : previous logarithmic elastic strain                  /上一步的弹性应变（自动求解）
%           - epsE   : logarithmic elastic strain                           /当前弹性应变（自动求解）
%           - mCst   : material constants (or internal state parameters)    /材料常数矩阵，每个物质点都可以不一样 [E,u]
%           - fp     : force at the material point                          /物质点处外力（自动求解）
%           - u      : material point displacement                          /物质点位移（自动求解）
%           - lp     : material point domain lengths                        /物质点影响范围 [x,y]（后续自动求解）
%           - lp0    : initial material point domain lengths                /初始影响范围 [x,y]
%
% mesh   - structured array with the following fields
%           - coord : mesh nodal coordinates (nodes,nD)                     /网格节点坐标
%           - etpl  : element topology (nels,nen)                           /网格单元构成（顺时针）
%           - bc    : boundary conditions (*,2)                             /边界条件，底部侧面固定位移
%           - h     : background mesh size (nD,1)                           /背景网格尺寸
%--------------------------------------------------------------------------
% See also:
% FORMCOORD2D - background mesh generation
% DETMPPOS    - local material point positions
% SHAPEFUNC   - background grid basis functions
%--------------------------------------------------------------------------

%% Analysis parameters
E      = 10e6;   v = 0.35;                                                  % Young's modulus, Poisson's ratio   
mCst   = [E v];                                                             % material constants
g      = 10;                                                                % gravity
rho    = 2000;                                                              % material density
lstps  = 40;                                                                % number of loadsteps
nelsx  = 4;                                                                 % number of elements in the x direction
nelsy  = 4;                                                                 % number of elements in the y direction
nelsz  = 32;                                                                % number of elements in the y direction
lz     = 50;  lx = lz/(nelsz/nelsx);     ly = lz/(nelsz/nelsy);             % domain dimensions
mp     = 2;                                                                 % number of material points in each direction per element
mpType = 2;                                                                 % material point type: 1 = MPM, 2 = GIMP
cmType = 1;                                                                 % constitutive model: 1 = elastic, 2 = vM plasticity

%% Mesh generation
[etpl,coord] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz);                     % background mesh generation
[nels,nen]   = size(etpl);                                                  % number of elements and nodes per element
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly lz]./[nelsx nelsy nelsz];                             % element lengths in each direction

%% Boundary conditions on background mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,1)==0 || coord(node,1)==lx                                  % roller sides
    bc(node*3-2,:)=[node*3-2 0];    
  end
  if coord(node,2)==0 || coord(node,2)==ly                                  % roller sides
    bc(node*3-1,:)=[node*3-1 0];    
  end
  if coord(node,3)==0                                                       % roller base
    bc(node*3  ,:)=[node*3   0];
  end
end
bc = bc(bc(:,1)>0,:);                                                       % remove empty part of bc

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
mesh.etpl  = etpl;                                                          % element topology
mesh.coord = coord;                                                         % nodal coordinates
mesh.bc    = bc;                                                            % boundary conditions
mesh.h     = h;                                                             % mesh size
mesh.eMin  = eMin;                                                          % element lower coordinate limit 
mesh.eMax  = eMax;                                                          % element upper coordainte limit 

%% Material point generation
ngp    = mp^nD;                                                             % number of material points per element
GpLoc  = detMpPos(mp,nD);                                                   % local MP locations (for each element)
N      = shapefunc(nen,GpLoc,nD);                                           % basis functions for the material points
[etplmp,coordmp] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz);                 % mesh for MP generation
nelsmp = size(etplmp,1);                                                    % no. elements populated with material points
nmp    = ngp*nelsmp;                                                        % total number of mterial points

mpC=zeros(nmp,nD);                                                          % zero MP coordinates
for nel=1:nelsmp
  indx=(nel-1)*ngp+1:nel*ngp;                                               % MP locations within mpC
  eC=coordmp(etplmp(nel,:),:);                                              % element coordinates
  mpPos=N*eC;                                                               % global MP coordinates
  mpC(indx,:)=mpPos;                                                        % store MP positions
end
lp = zeros(nmp,nD);                                                         % zero domain lengths
lp(:,1) = h(1)/(2*mp);                                                      % domain half length x-direction
lp(:,2) = h(2)/(2*mp);                                                      % domain half length y-direction
lp(:,3) = h(3)/(2*mp);                                                      % domain half length z-direction
vp      = 2^nD*lp(:,1).*lp(:,2).*lp(:,3);                                   % volume associated with each material point

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
  mpData(mp).nSMe   = 0;                                                    % number of stiffness entries associated with the material point
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
  mpData(mp).mCst   = mCst;                                                 % material constants (or internal variables) for constitutive model
  mpData(mp).fp     = zeros(nD,1);                                          % point forces at material points
    mpData(mp).u      = zeros(nD,1);                                          % material point displacements
  mpData(mp).du     = zeros(nD,1);
  mpData(mp).ddu     = zeros(nD,1);
  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                           % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                           % initial material point domain lengths (GIMP)
  else
    mpData(mp).lp     = zeros(1,nD);                                        % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                        % initial material point domain lengths (MPM)
  end
end
end