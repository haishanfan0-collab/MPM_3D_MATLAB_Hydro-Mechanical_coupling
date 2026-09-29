function [mpData,mesh] = DAN_Zhidian_TEST(Load,YSD,HSL)

% 【核心-路基-地基耦合结构物质点法前处理模型】
% 基于GIMP物质点法的二维路基-地基耦合体系几何建模、材料赋值、边界条件
% 与初始场设置，涵盖路面结构层、路基填土、天然地基的非饱和-饱和分区、
% 循环荷载施加、渗流边界配置等全过程，为湿-力耦合物质点法分析提供
% 完整初始计算状态
%--------------------------------------------------------------------------
% 作者: FAN Haishan（基于William Coombs原始框架重构）
% 日期: 29/04/2026
% 描述:
% 进行路基-地基耦合体系的物质点法前处理建模，采用分层分区策略实现：
%   1. 路面结构层：沥青面层（线弹性，高模量）+ 基层（线弹性，中模量）
%   2. 路基填土：非饱和区（大周期弹塑性边界面模型，cmType=5）+
%               饱和区（线弹性，cmType=1，考虑浮力效应）
%   3. 天然地基：非饱和区（Mohr-Coulomb模型，cmType=3）+
%               饱和区（线弹性，cmType=1，高饱和度软塑状态）
% 特殊功能：循环车辆荷载自动施加（四组双轮荷载，200kN/组）、
%          降雨边界配置（坡面+顶面入渗）、吸力-饱和度初始场计算、
%          渗流边界条件（底面零水头、坡面/顶面降雨强度）
%
% 几何尺度：路基宽度24m，填高8m，地基深度4m（总计算域36m×14m）
% 路面结构：面层0.18m + 基层0.20m + 路基1.50m（总厚度1.88m）
%--------------------------------------------------------------------------
% [mpData,mesh] = pavement_structure_couple
%--------------------------------------------------------------------------
% 输入: 无（纯参数化建模，所有参数内部定义）
%--------------------------------------------------------------------------
% 输出:
% mpData      - 物质点信息结构体数组 [1×nmp]，每个元素含字段：
%   力学核心字段:
%     mpType        - 物质点类型: 2=GIMP（广义插值物质点）
%     cmType        - 本构模型类型: 1=线弹性, 3=Mohr-Coulomb, 5=边界面模型
%     mpC           - 物质点坐标 [x, y] [m]
%     vp/vp0        - 当前/初始体积 [m³]
%     mpM           - 物质点质量 [kg]（含含水质量）
%     nIN/eIN       - 关联节点/单元编号
%     Svp/dSvp/Gvp/Tvp - 形函数及其导数、GIMP域函数、时间导数
%     Fn/F          - 上一步/当前变形梯度 [3×3]
%     sig/sign      - 当前/上一步Cauchy应力 [6×1] [Pa]（Voigt记法）
%     epsE/epsEn    - 当前/上一步弹性应变 [6×1]
%     epsP/epsPn    - 当前/上一步塑性应变 [6×1]
%     eps/epsn      - 当前/上一步总应变 [6×1]
%     mCst          - 材料常数向量（依cmType维度不同）
%     fp            - 外力向量 [fx, fy] [N]
%     u/du/ddu      - 位移/增量/二次增量 [m]
%     lp/lp0        - 当前/初始域长度 [m]（GIMP特性）
%
%   非饱和渗流字段:
%     Flow_P        - 孔隙水压力 [Pa]
%     Flow_Type     - 渗流模型: 0=不透水, 1=VG模型, 99=饱和区
%     Flow_CP       - SWCC参数 [Fa, Fm, Fn]
%     Flow_Position - 边界位置标记: 'Slope'/'Up'/'Right'等
%     Flow_SIZE     - 边界流量 [m/s]（降雨强度）
%     dH            - 水头增量 [m]
%     H0            - 初始水头 [m]
%     dH_sum        - 累积水头 [m]
%     dH_sum_strain - 应变诱发水头 [m]
%     Gravity_Water_Content      - 当前质量含水率 [-]
%     Gravity_Water_Content_old  - 上一步质量含水率 [-]
%     Gs            - 土粒比重 [-]
%     Roumax        - 最大干密度 [kg/m³]
%     Roud          - 当前干密度 [kg/m³]
%     e0/e          - 初始/当前孔隙比 [-]
%     Sr            - 饱和度 [-]
%     Scution_Flow  - 基质吸力 [kPa]
%     Ksat/Ksat0    - 当前/初始饱和渗透系数 [m/s]（水平;垂直）
%     eps_Flow/sig_Flow      - 渗流应变/应力 [2×1]
%     epsn_Flow/sign_Flow    - 上一步渗流应变/应力 [2×1]
%     Pore_Pressure - 孔隙水压力 [Pa]
%     Saturated_YES - 饱和标志: 0=非饱和, 1=饱和
%
%   状态变量字段:
%     STATEV/STATEVn - 当前/上一步状态变量 [19×1]（边界面模型专用）
%     modulus        - 当前割线模量 [MPa]（用于输出诊断）
%
% mesh        - 背景网格信息结构体，含字段：
%     etpl    - 单元拓扑 [nels×nen]
%     coord   - 节点坐标 [nodes×nD]
%     bc      - 位移边界条件 [nbc×2]（自由度编号, 固定值）
%     bc_Flow - 渗流边界条件 [nbc_f×2]（节点编号, 水头值）
%     h       - 单元尺寸 [nD×1] [m]
%     eMin/eMax - 单元坐标上下限 [nels×nD]
%     g       - 重力加速度 [m/s²]
%--------------------------------------------------------------------------
% 材料参数配置:
% 边界面模型参数 mcSM (29维向量，详见Boundary_module_BX注释):
%   mcSM_93: 压实度K0=0.93（路基核心区，高模量）
%   mcSM_94: 压实度K0=0.94（路基过渡区，中模量）
%   mcSM_96: 压实度K0=0.96（路基边坡区，低模量）
%
% 线弹性材料参数:
%   沥青面层: [E=3000MPa, ν=0.25]（cmType=1）
%   基层:     [E=80MPa, ν=0.35]（cmType=1）
%   饱和地基: [E=20MPa, ν=0.45]（cmType=1，地下水位以下）
%
% Mohr-Coulomb材料参数:
%   非饱和地基: [E=40MPa, ν=0.40, C=12kPa, φ=10°, ψ=-1.5°, Gs=2660]
%--------------------------------------------------------------------------
% 荷载配置:
% 循环车辆荷载: 四组双轮矩形荷载，轮压0.3m×0.3m，中心间距1.8m
%   轮位1: x∈[2.575, 2.875]m, y∈[11.82, 12.12]m（左轮）
%   轮位2: x∈[3.775, 4.075]m, y∈[11.82, 12.12]m（左轮）
%   轮位3: x∈[6.325, 6.625]m, y∈[11.82, 12.12]m（右轮）
%   轮位4: x∈[7.525, 7.825]m, y∈[11.82, 12.12]m（右轮）
% 单轮荷载力: F = 200kN / (4轮×接触面积) = 200000/A [N/m²]（等效均布）
%--------------------------------------------------------------------------
% 边界条件:
% 位移边界:
%   - 左边界(x=0):    x方向固定（滚轴约束）
%   - 右边界(x=36m):  x方向固定（滚轴约束）
%   - 底边界(y=0):    y方向固定（铰支约束）
% 渗流边界:
%   - 底面(y≤1.0m):   零水头边界（排水层）
%   - 坡面:           降雨入渗（Flow_Position='Slope'）
%   - 顶面(路面范围):  降雨入渗（Flow_Position='Up'，路面区域设为极小流量）
%--------------------------------------------------------------------------
% 建模流程:
% 1. 定义三层边界面模型参数（不同压实度）
% 2. 生成背景网格（36m×14m，18a×7a单元，a=7倍加密因子）
% 3. 设置位移/渗流边界条件
% 4. 生成路基物质点（GIMP，3×3/单元，分层赋值cmType=5）
% 5. 删除边坡外多余物质点，生成边界虚拟物质点
% 6. 路面结构层重赋值（沥青/基层线弹性）
% 7. 生成地基物质点（GIMP，2×2/单元，分区赋值cmType=1/3）
% 8. 施加循环车辆荷载（四组双轮，200kN等效）
%--------------------------------------------------------------------------
% 调用关系:
% 被调用: MPM_Solve_Main_BX（主求解器，作为初始配置函数）
% 内部调用: formCoord2D（网格生成）、detMpPos（物质点局部坐标）、
%           shapefunc（形函数计算）、boundary_mp_generate（边界虚拟MP）、
%           boundary_Flow（渗流边界配置）
%--------------------------------------------------------------------------

%% Analysis parameters
% mcSM = [1.80	 0.12	1.70...
%     0.12	0.01	0.8 	1.58	1.0	1.8	8.8	1.56	...
%     0.194	2.0	1730.0	2600.0	0.93	0.35	1.42	...
%     0.44	0.88	-1.27	9.37	-0.06	-4.0	0.67 ...	
%     1.5	103000  103000	2e-05	0.2016	0.2016	0.00];  
mcSM = [1.94	 0.16	1.84...
    0.16	0.01	0.8 	1.68	1.0	2.0	16.8	2.16	...
    0.181	2.0	1730.0	2600.0	0.93	0.35	1.42	...
    0.44	0.88	-1.27	9.37	-0.06	-4.0	0.67 ...	
    1.5	103000  103000	2e-05	0.2016	0.2016	0.00];  
mcSM(16)=YSD;
mcSM(end-2)=HSL/100;
mcSM(end-1)=HSL/100;

ly=1;  lx=1; lz=1;                                                          % domain dimensions
g=9.80;                                                                     % gravity
a = 2;                                                                      % element multiplier
nelsx=a;                                                                    % number of elements in the x direction
nelsy=a;                                                                    % number of elements in the y direction
nelsz=a;                                                                    % number of elements in the y direction
mp_number=2;                                                                % number of material points in each direction per element

%% Mesh generation
[etpl,coord] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz);         % background mesh generation
[nels,nen]   = size(etpl);                                                  % number of elements and nodes per element
[nodes,nD]   = size(coord);                                                 % number of nodes and dimensions
h            = [lx ly lz]./[nelsx nelsy nelsy];                             % element lengths in each direction

%% Boundary conditions on backgroun mesh
bc = zeros(nodes*nD,2);                                                     % generate empty bc matrix
for node=1:nodes                                                            % loop over nodes
  if coord(node,2)<=0.001                                                  % roller (y=0)
    bc(node*3-1  ,:)=[node*3-1   0];
  end

  if coord(node,1)==0.5*lx  && coord(node,2)<=0.001                           % roller (y=0)
      bc(node*3-2  ,:)=[node*3-2   0];
  end
  if coord(node,3)==0.5*lz  && coord(node,2)<=0.001                           % roller (y=0)
      bc(node*3  ,:)=[node*3   0];
  end
end
bc = bc(bc(:,1)>0,:); 

%% Boundary Flow conditions on backgroun mesh
bc_Flow = zeros(nodes,2);                                                   % generate empty bc matrix
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

%% MPM generation
ngp    = mp_number^nD;                                                      % number of material points per element
GpLoc  = detMpPos(mp_number,nD);                                            % local MP locations
N      = shapefunc(nen,GpLoc,nD);                                           % basis functions for the material points
[etplmp,coordmp] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz);                 % mesh for MP generation
% coordmp(:,1) = coordmp(:,1)+0.5*lx;
% coordmp(:,2) = coordmp(:,2)+0.5*ly;
% coordmp(:,3) = coordmp(:,3)+0.5*lz;

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
lp(:,1) = h(1)/(2*mp_number);                                               % domain half length x-direction
lp(:,2) = h(2)/(2*mp_number);                                               % domain half length y-direction
lp(:,3) = h(2)/(2*mp_number);                                               % domain half length y-direction
vp      = 2^nD*lp(:,1).*lp(:,2).*lp(:,3);                                   % volume associated with each material point

%% Material point structure generation
for mp = nmp:-1:1                                                           % loop backwards over MPs so array doesn't change size
  % 基础信息
  mpData(mp).mpType = 2;                                                    % material point type: 1 = MPM, 2 = GIMP
  mpData(mp).cmType = 5;                                                    % constitutive model: 1 = elastic, 2 = vM plasticity
  mpData(mp).mpC    = mpC(mp,:);                                            % material point coordinates
  mpData(mp).vp     = vp(mp);                                               % material point volume
  mpData(mp).vp0    = vp(mp);                                               % material point initial volume
  mpData(mp).nIN    = zeros(nen,1);                                         % nodes associated with the material point
  mpData(mp).eIN    = 0;                                                    % element associated with the material point
  mpData(mp).Svp    = zeros(1,nen);                                         % material point basis functions
  mpData(mp).dSvp   = zeros(nD,nen);                                        % derivative of the basis functions
  mpData(mp).Gvp  = zeros(nD,nen);                                          % basis function derivatives
  mpData(mp).Tvp  = zeros(1,nen);
  mpData(mp).Fn     = eye(3);                                               % previous deformation gradient
  mpData(mp).F      = eye(3);                                               % deformation gradient
  mpData(mp).sig    = zeros(6,1);                                           % Cauchy stress
  mpData(mp).sign    = zeros(6,1);                                          % Cauchy stress
  mpData(mp).epsEn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsE   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).epsPn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
  mpData(mp).epsP   = zeros(6,1);                                           % elastic strain (logarithmic)
  mpData(mp).mCst   = mcSM;
  mpData(mp).Ksat = [2e-8;1e-8];                                            % 饱和渗透系数
  mpData(mp).fp     = zeros(nD,1);                                          % point forces at material points
  mpData(mp).u      = zeros(nD,1);                                          % material point displacements
  mpData(mp).du     = zeros(nD,1);
  mpData(mp).ddu    = zeros(nD,1);
  mpData(mp).STATEV    = zeros(20,1);
  mpData(mp).STATEVn    = zeros(20,1);
  if mpData(mp).mpType == 2
    mpData(mp).lp     = lp(mp,:);                                           % material point domain lengths (GIMP)
    mpData(mp).lp0    = lp(mp,:);                                           % initial material point domain lengths (GIMP)
  else
    mpData(mp).lp     = zeros(1,nD);                                        % material point domain lengths (MPM)
    mpData(mp).lp0    = zeros(1,nD);                                        % initial material point domain lengths (MPM)
  end
  mpData(mp).modulus=0.0;

  %% 非饱和渗流参数
  mpData(mp).Flow_P = 0.0;
  mpData(mp).Flow_Type = 1;                                                 % 基质吸力模型，1-Van Genuchten模型
  mpData(mp).dH     = 0.0;
  % 基质吸力模型参数
  mpData(mp).Flow_CP = [mpData(mp).mCst(27),mpData(mp).mCst(25),mpData(mp).mCst(26)];
  mpData(mp).Gravity_Water_Content = mpData(mp).mCst(end-2);                % 初始质量含水率 0.14~0.18
  try
      mpData(mp).Gravity_Water_Content_old = mpData(mp).Gravity_Water_Content;
  catch
  end
  mpData(mp).Gs = mpData(mp).mCst(15);                                      % 比重
  mpData(mp).Roumax = mpData(mp).mCst(14);                                  % 最大干密度
  rho = mpData(mp).mCst(16)*mpData(mp).Roumax;
  rhoA = rho*(1+mpData(mp).Gravity_Water_Content);
  mpData(mp).mpM = vp(mp)*rhoA;                                             % 重新获取质量
  mpData(mp).e0 = mpData(mp).Gs/rho-1;                                      % 初始空隙率
  mpData(mp).e = mpData(mp).e0;                                             % 孔隙率
  mpData(mp).Sr = mpData(mp).Gravity_Water_Content*mpData(mp).Gs/1000/mpData(mp).e; % 饱和度计算
  mpData(mp).eps_Flow = zeros(2,1);                                         % 饱和渗透系数
  mpData(mp).sig_Flow = zeros(2,1);                                         % 饱和渗透系数
  mpData(mp).epsn_Flow = zeros(2,1);                                        % 饱和渗透系数
  mpData(mp).sign_Flow = zeros(2,1);                                        % 饱和渗透系数
  Fa = mpData(mp).Flow_CP(1)/1000;
  Fm = mpData(mp).Flow_CP(2);
  Fn = mpData(mp).Flow_CP(3);
  if mpData(mp).Sr>1.0
      error('   Error')
  end
  mpData(mp).Scution_Flow = ((mpData(mp).Sr^(-1/Fm) - 1)^(1/Fn))*Fa/mpData(mp).e^(1/(Fm*Fn));   % 饱和渗透系数
  mpData(mp).H0      = mpData(mp).Scution_Flow/-9.81;
  mpData(mp).dH_sum = 0.0;
  mpData(mp).dH_sum_strain = 0.0;
  mpData(mp).Roud   = rho;
  mpData(mp).Saturated_YES   = 0;
  mpData(mp).STATEVn(7) = mpData(mp).Scution_Flow;
  mpData(mp).Pore_Pressure = 0.0;
end

mpC_mp = reshape([mpData.mpC],3,length(mpData))';

Mian_SxUp = (mpC_mp(:,1)>=(lx-1e-4-lp(:,1)));
Mian_SxDown = (mpC_mp(:,1)<=(0+1e-4+lp(:,1)));

Mian_SyUp = (mpC_mp(:,2)>=(ly-1e-4-lp(:,2)));

Mian_SzUp = (mpC_mp(:,3)>=(lz-1e-4-lp(:,3)));
Mian_SzDown = (mpC_mp(:,3)<=(0+1e-4+lp(:,3)));        

% 荷载施加
Load_WY = Load(1)*1000*(lx*ly*2+ly*lz*2+lx*lz)/...
    (sum(Mian_SxUp)+sum(Mian_SxDown)+sum(Mian_SyUp)+sum(Mian_SzUp)+sum(Mian_SzDown));
Load_JC = Load(2)*1000*(lx*lz)/sum(Mian_SyUp);
Load_YH = Load(3)*1000*(lx*lz)/sum(Mian_SyUp);

Contact = zeros(size(mpData,2),nD);

for i = 1:size(mpData,2)
    
    % 右表面
    if Mian_SxUp(i)
        Contact(i,1) = -Load_WY;
    end

    % 左表面
    if Mian_SxDown(i)
        Contact(i,1) = Load_WY;
    end

    % 后表面
    if Mian_SzUp(i)
        Contact(i,3) = -Load_WY;
    end

    % 前表面
    if Mian_SzDown(i)
        Contact(i,3) = Load_WY;
    end

    % 上表面
    if Mian_SyUp(i)
        mpData(i).fp = [0;-Load_YH;0];
        Contact(i,2) = -Load_WY-Load_JC;
    end
    
end

mesh.Contact_Force = Contact;

nmp  = length(mpData);
for mp = nmp:-1:1
    mpData(mp).Pore_Pressure = 0.0;
    mpData(mp).eps = zeros(6,1);
    mpData(mp).epsn = zeros(6,1);
    mpData(mp).Ksat0 = mpData(mp).Ksat;
end
end
