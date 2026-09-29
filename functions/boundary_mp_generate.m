function  mpData = boundary_mp_generate(mpData,nen,mCst,pavement_layer,rho)

% 【核心-GIMP边界虚拟物质点生成器】
% 基于物质点域边界轮廓识别的高精度边界虚拟物质点自动生成，用于消除GIMP物质点法
% 在自由边界处的形函数截断误差，提升边界应力/位移计算精度与接触问题处理可靠性
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 03/03/2026
% 描述:
% 进行GIMP物质点法边界虚拟物质点（Boundary Material Points, BMP）的自动识别与生成。
% 核心思想源自文献：A high-fidelity material point method for frictional contact problems
% (arXiv:2403.13534v1 [math.NA] 20 Mar 2024)。
% 
% 算法流程：
%   1. 边界轮廓识别：遍历所有物质点域角点（1D:2个, 2D:8个, 3D:26个），提取边界坐标
%   2. 冗余点删除：通过unique去重，按坐标排序后识别极值点（上下/左右包络）
%   3. 边界点精简：删除近似重合点，保留有效边界控制点
%   4. 虚拟物质点生成：在边界控制点位置生成mpType=3的虚拟物质点，继承内部点属性
% 
% 虚拟物质点特性：
%   - 不参与应力-应变计算（无本构更新），仅用于边界形函数完整性修正
%   - 质量/体积极小（约为内部点的0.1%），避免对整体质量矩阵产生显著影响
%   - 继承同层材料参数，确保边界物理属性一致性
%   - 支持渗流场耦合（若内部点含Flow_Type等渗流字段）
%
%--------------------------------------------------------------------------
% mpData = boundary_mp_generate(mpData, nen, mCst, FHS_AHA3, rho)
%--------------------------------------------------------------------------
% 输入:
% mpData    - 物质点信息结构体数组 [1×mpM]（输入为内部物质点，输出追加虚拟点）
% nen       - 单元节点数 [-]（用于初始化形函数数组维度）
% mCst      - 材料参数矩阵 [nLayer×nPar] 或 [1×nPar]
%             若为多层：按y坐标分层赋值（第1行底层, 第2行中层, 第3行上层）
%             若为单层：所有边界点统一赋值
% FHS_AHA3  - 预留参数（当前未使用，保持接口兼容性）
% rho       - 参考密度 [kg/m³]（用于虚拟点质量计算）
%--------------------------------------------------------------------------
% 输出:
% mpData    - 更新后的物质点结构体数组 [1×(mpM+nBMP)]，追加字段：
%             新增虚拟物质点索引：mpM+1 : mpM+nBMP
%             虚拟点关键字段：
%               mpType = 3          — 边界虚拟物质点标识
%               mpC                  — 边界控制点坐标 [x,y] 或 [x,y,z]
%               vp/vp0 = mpvv        — 极小体积（约为内部点总体积的0.1%/nBMP）
%               mpM = mpvv*rho       — 极小质量（避免质量矩阵奇异）
%               lp/lp0 = mpll       — 虚拟域长度（由各向同性体积反推）
%               cmType               — 继承内部点本构类型（同层）
%               mCst                 — 继承内部点材料参数（同层分层赋值）
%               其余力学/渗流字段   — 继承或零初始化（详见代码）
%--------------------------------------------------------------------------
% 维度支持:
% 1D: 每物质点生成2个边界点（左右端点）
% 2D: 每物质点生成8个边界点（4边中点+4角点），经包络提取后精简
% 3D: 每物质点生成26个边界点（6面心+12棱中点+8角点），经包络提取后精简
%
% 2D边界识别算法细节:
%   - 按x排序提取每个x坐标对应的y_max/y_min（上下包络）
%   - 按y排序提取每个y坐标对应的x_max/x_min（左右包络）
%   - 合并上下左右包络点，去重后删除近似重合点（容差1e-5）
%   - 最终保留的边界控制点即为虚拟物质点位置
%
% 材料分层赋值逻辑（2D/3D）:
%   若mCst为多层矩阵（size(mCst,1)>1），按y坐标分层：
%     y ≥ 1.5*ly - 0.80 - 0.78  →  mCst(3,:)  （上层，路面结构区）
%     y ≥ 1.5*ly - 1.50 - 0.80 - 0.78  →  mCst(2,:)  （中层，路基核心区）
%     其他  →  mCst(1,:)  （下层，地基区）
%   其中ly=8m为路基高度，分层厚度对应路面结构层配置
%--------------------------------------------------------------------------
% 关键参数:
% mpV_all   - 所有内部物质点总体积 [m³]
% mpvv      - 单个虚拟物质点体积 = mpV_all × 0.001 / nBMP [m³]
% mpll      - 虚拟域半长度 = (mpvv / 2^nD)^(1/nD) [m]
% nBMP      - 生成的虚拟物质点数量（由边界几何复杂度决定）
%--------------------------------------------------------------------------
% 渗流场兼容性:
% 若输入mpData(1)含Flow_Type字段，虚拟点自动继承渗流参数：
%   Flow_Type, Flow_CP, Gravity_Water_Content, Gs, Roumax, Ksat, Ksat0
%   初始吸力Scution_Flow由SWCC公式基于初始饱和度计算
%   H0 = Scution_Flow / -9.81（初始水头，吸力转正水头）
%--------------------------------------------------------------------------
% 调用关系:
% 被调用: pavement_structure_couple（路基建模主函数，用于边坡边界处理）
% 内部调用: 无（纯MATLAB基础函数：unique, sortrows, vecnorm等）
%--------------------------------------------------------------------------
mpM = size(mpData,2);
% 计算总的物质点体积
mpV_all = sum([mpData.vp]);
ly=8;
% 物质点边缘坐标识别
nD = length(mpData(1).mpC);

% 一维节点识别
if nD == 1
    boundary_MAT = zeros(mpM*2,nD);
    for mp = 1:mpM
        boundary_MAT((mp-1)*2+1,1) = mpData(mp).mpC-mpData(mp).lp0;
        boundary_MAT((mp-1)*2+2,1) = mpData(mp).mpC+mpData(mp).lp0;
    end
end

% 二维节点识别
if nD == 2
    boundary_MAT = zeros(mpM*8,nD);
    for mp = 1:mpM
        % 上下左右
        boundary_MAT((mp-1)*8+1,:) = mpData(mp).mpC-[0,mpData(mp).lp0(2)];
        boundary_MAT((mp-1)*8+2,:) = mpData(mp).mpC+[0,mpData(mp).lp0(2)];
        boundary_MAT((mp-1)*8+3,:) = mpData(mp).mpC-[mpData(mp).lp0(1),0];
        boundary_MAT((mp-1)*8+4,:) = mpData(mp).mpC+[mpData(mp).lp0(1),0];
        % 角点
        boundary_MAT((mp-1)*8+5,:) = mpData(mp).mpC+[mpData(mp).lp0(1),mpData(mp).lp0(2)];
        boundary_MAT((mp-1)*8+6,:) = mpData(mp).mpC+[mpData(mp).lp0(1),-mpData(mp).lp0(2)];
        boundary_MAT((mp-1)*8+7,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),mpData(mp).lp0(2)];
        boundary_MAT((mp-1)*8+8,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),-mpData(mp).lp0(2)];
    end
end

% 三维节点识别
if nD == 3
    boundary_MAT = zeros(mpM*26,nD);
    for mp = 1:mpM
        % 1
        boundary_MAT((mp-1)*26+1,:) =  mpData(mp).mpC-[mpData(mp).lp0(1),0,0];
        boundary_MAT((mp-1)*26+2,:) =  mpData(mp).mpC+[mpData(mp).lp0(1),0,0];
        boundary_MAT((mp-1)*26+3,:) =  mpData(mp).mpC-[0,mpData(mp).lp0(2),0];
        boundary_MAT((mp-1)*26+4,:) =  mpData(mp).mpC+[0,mpData(mp).lp0(2),0];
        boundary_MAT((mp-1)*26+5,:) =  mpData(mp).mpC-[0,0,mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+6,:) =  mpData(mp).mpC+[0,0,mpData(mp).lp0(3)];

        % 2
        boundary_MAT((mp-1)*26+7,:) =  mpData(mp).mpC+[ mpData(mp).lp0(1), mpData(mp).lp0(2),0];
        boundary_MAT((mp-1)*26+8,:) =  mpData(mp).mpC+[-mpData(mp).lp0(1), mpData(mp).lp0(2),0];
        boundary_MAT((mp-1)*26+9,:) =  mpData(mp).mpC+[ mpData(mp).lp0(1),-mpData(mp).lp0(2),0];
        boundary_MAT((mp-1)*26+10,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),-mpData(mp).lp0(2),0];

        boundary_MAT((mp-1)*26+11,:) = mpData(mp).mpC+[ mpData(mp).lp0(1),0, mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+12,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),0, mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+13,:) = mpData(mp).mpC+[ mpData(mp).lp0(1),0,-mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+14,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),0,-mpData(mp).lp0(3)];

        boundary_MAT((mp-1)*26+15,:) = mpData(mp).mpC+[0, mpData(mp).lp0(2), mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+16,:) = mpData(mp).mpC+[0,-mpData(mp).lp0(2), mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+17,:) = mpData(mp).mpC+[0, mpData(mp).lp0(2),-mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+18,:) = mpData(mp).mpC+[0,-mpData(mp).lp0(2),-mpData(mp).lp0(3)];

        % 3
        boundary_MAT((mp-1)*26+19,:) = mpData(mp).mpC+[ mpData(mp).lp0(1), mpData(mp).lp0(2), mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+20,:) = mpData(mp).mpC+[-mpData(mp).lp0(1), mpData(mp).lp0(2), mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+21,:) = mpData(mp).mpC+[ mpData(mp).lp0(1),-mpData(mp).lp0(2), mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+22,:) = mpData(mp).mpC+[ mpData(mp).lp0(1), mpData(mp).lp0(2),-mpData(mp).lp0(3)];

        boundary_MAT((mp-1)*26+23,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),-mpData(mp).lp0(2), mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+24,:) = mpData(mp).mpC+[-mpData(mp).lp0(1), mpData(mp).lp0(2),-mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+25,:) = mpData(mp).mpC+[ mpData(mp).lp0(1),-mpData(mp).lp0(2),-mpData(mp).lp0(3)];
        boundary_MAT((mp-1)*26+26,:) = mpData(mp).mpC+[-mpData(mp).lp0(1),-mpData(mp).lp0(2),-mpData(mp).lp0(3)];
    end
end

% 删除重复坐标
boundary_MAT = unique(boundary_MAT, 'rows');

% 二维节点边界识别
if nD == 2
    % 按x找最大最小y
    LS1 = find_ChongFu_2D(boundary_MAT,nD,1);
    % 按y找最大最小x
    LS2 = find_ChongFu_2D(boundary_MAT,nD,2);
    % 合并
    LS3 = [LS1;LS2];
    LS3 = unique(LS3, 'rows');
    
    % 重复点删除
    LS_Check = [(1:length(LS3))',vecnorm(LS3, 2, 2),LS3(:,1)];
    LS_Check = sortrows(LS_Check,2);
    LS_Check(:,2) = [0;LS_Check(2:end,2)-LS_Check(1:end-1,2)];
    LS_Check(:,3) = [0;LS_Check(2:end,3)-LS_Check(1:end-1,3)];
    LS_Check = [LS_Check,vecnorm(LS_Check(:,2:3), 2, 2)];
    LS_Check = sortrows(LS_Check,4);
    LS_Check(1,:) = [];
    LS3(LS_Check(LS_Check(:,4)<1e-5,1),:)=[];
end

% 三维节点边界识别
if nD == 3
    boundary_MAT_sum = (1+boundary_MAT(:,1)).^(0.1+boundary_MAT(:,2)+min(boundary_MAT(:,2))).*boundary_MAT(:,3);
    [~, IA, ~] = uniquetol(boundary_MAT_sum,1e-12);
    boundary_MAT = boundary_MAT(IA,:);
    boundary_MAT = unique(boundary_MAT, 'rows');
    % x,z → y
    LS1 = find_ChongFu_3D(boundary_MAT,nD,2);
    % y,z → x
    LS2 = find_ChongFu_3D(boundary_MAT,nD,1);
    % x,y → z
    LS3 = find_ChongFu_3D(boundary_MAT,nD,3);
    LS3 = [LS1;LS2;LS3];
    % 删除重复项
    LS3_sum = (1+LS3(:,1)).^(0.1+LS3(:,2)+min(LS3(:,2))).*LS3(:,3);
    [~, IA, ~] = uniquetol(LS3_sum,1e-12);
    LS3 = LS3(sortrows(IA,1),:);
end

%% 增加物质点
mpvv = mpV_all*0.001/size(LS3,1);
mpll = (mpvv/(2^nD))^(1/nD);
for mp = mpM+size(LS3,1):-1:mpM+1                                           % loop backwards over MPs so array doesn't change size
    mpData(mp).mpM = mpvv*rho;
    mpData(mp).mpType = 3;                                                    % material point type: 1 = MPM, 2 = GIMP, 3 = boundryMP
    mpData(mp).cmType = mpData(1).cmType;                                     % constitutive model: 1 = elastic, 2 = vM plasticity
    mpData(mp).mpC    = LS3(mp-mpM,:);                                        % material point coordinates
    mpData(mp).vp     = mpvv;                                                 % material point volume
    mpData(mp).vp0    = mpvv;                                                 % material point initial volume
    mpData(mp).nIN    = zeros(nen,1);                                         % nodes associated with the material point
    mpData(mp).eIN    = 0;                                                    % element associated with the material point
    mpData(mp).Svp    = zeros(1,nen);                                         % material point basis functions
    mpData(mp).dSvp   = zeros(nD,nen);                                        % derivative of the basis functions
    try
        mpData(mp).Gvp  = mpData(1).Gvp;                                                 % basis function derivatives
        mpData(mp).Tvp  = mpData(1).Tvp;
    catch
    end
    mpData(mp).Fn     = eye(3);                                               % previous deformation gradient
    mpData(mp).F      = eye(3);                                               % deformation gradient
    mpData(mp).sig    = zeros(6,1);                                           % Cauchy stress
    mpData(mp).sign   = zeros(6,1);                                           % Cauchy stress
    mpData(mp).epsEn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
    mpData(mp).epsE   = zeros(6,1);                                           % elastic strain (logarithmic)
    mpData(mp).epsPn  = zeros(6,1);                                           % previous elastic strain (logarithmic)
    mpData(mp).epsP   = zeros(6,1);                                           % elastic strain (logarithmic)
    if size(mCst,1)==1
        mpData(mp).mCst   = mCst;
        mpData(mp).Ksat = [2e-8;1e-8];                     % 饱和渗透系数
    else
        if mpData(mp).mpC(2)>=1.50*ly - pavement_layer(1) - pavement_layer(2)
            mpData(mp).mCst   = mCst(4,:);
            mpData(mp).Ksat = [1e-7;2e-7];                 % 饱和渗透系数
        elseif mpData(mp).mpC(2)>=1.50*ly - pavement_layer(1) - pavement_layer(2) - pavement_layer(3)
            mpData(mp).mCst   = mCst(3,:);
            mpData(mp).Ksat = [5e-7;1e-6];                 % 饱和渗透系数
        else
            mpData(mp).mCst   = mCst(2,:);
            mpData(mp).Ksat = [1e-6;5e-6];                 % 饱和渗透系数
        end

        if (mpData(mp).mpC(1)*2/3+(mpData(mp).mpC(2)-ly*0.5)>48/3-13^0.5/6)
            mpData(mp).mCst   = mCst(1,:);
            mpData(mp).Ksat = [5e-6;1e-5];                 % 饱和渗透系数
        end

        mpData(mp).Roumax = mpData(mp).mCst(14);
        rho = mpData(mp).mCst(16)*mpData(mp).Roumax;
    end
    mpData(mp).fp     = zeros(nD,1);                                          % point forces at material points
    mpData(mp).u      = zeros(nD,1);                                          % material point displacements
    mpData(mp).du     = zeros(nD,1);
    mpData(mp).ddu    = zeros(nD,1);

    mpData(mp).lp     = mpll*ones(1,nD);                                      % material point domain lengths (MPM)
    mpData(mp).lp0    = mpll*ones(1,nD);                                      % initial material point domain lengths (MPM)
    mpData(mp).Flow_P = 0.0;
    mpData(mp).STATEV    = mpData(1).STATEV;
    mpData(mp).STATEVn    = mpData(1).STATEVn;
    mpData(mp).Scution_Flow = 1.0;
    mpData(mp).Gravity_Water_Content = 0.0;
    mpData(mp).dH_sum = 0.0;
    mpData(mp).dH_sum_strain=0.0;
    mpData(mp).modulus=0.0;

    % 渗流参数
    Name_mp = fieldnames(mpData);
    if sum(strcmp(Name_mp,'Flow_Type'))>0
        mpData(mp).Flow_P = 0.0;
        mpData(mp).Flow_Type = 1;                                                 % 基质吸力模型，1-Van Genuchten模型
        mpData(mp).dH     = 0.0;
        mpData(mp).Flow_CP = [mpData(mp).mCst(27),mpData(mp).mCst(25),mpData(mp).mCst(26)];                                   % 基质吸力模型参数
        mpData(mp).Gravity_Water_Content = mpData(mp).mCst(end-2);                % 初始质量含水率 0.14~0.18
        mpData(mp).Gravity_Water_Content_old = mpData(mp).Gravity_Water_Content;
        
        mpData(mp).Gs = mpData(mp).mCst(15);                                      % 比重
        mpData(mp).Roumax = mpData(mp).mCst(14);                                  % 最大干密度
        rhoA = rho*(1+mpData(mp).Gravity_Water_Content);
        mpData(mp).mpM = mpvv*rhoA;                                               % 重新获取质量
        mpData(mp).e0 = mpData(mp).Gs/rho-1;                                      % 初始空隙率
        mpData(mp).e = mpData(mp).e0;                                             % 孔隙率
        mpData(mp).Sr = mpData(mp).Gravity_Water_Content*mpData(mp).Gs/1000/mpData(mp).e; % 饱和度计算
        mpData(mp).eps_Flow = zeros(2,1);                                         % 饱和渗透系数
        mpData(mp).sig_Flow = zeros(2,1);                                         % 饱和渗透系数
        mpData(mp).epsn_Flow = zeros(2,1);                                        % 饱和渗透系数
        mpData(mp).sign_Flow = zeros(2,1);                                        % 饱和渗透系数
        Fa = mpData(mp).Flow_CP(1);
        Fm = mpData(mp).Flow_CP(2);
        Fn = mpData(mp).Flow_CP(3);
        mpData(mp).Scution_Flow = ((mpData(mp).Sr^(-1/Fm) - 1)^(1/Fn))*Fa/mpData(mp).e^(1/(Fm*Fn))/1000;   % 饱和渗透系数
        mpData(mp).H0      = mpData(mp).Scution_Flow/-9.81;
        mpData(mp).dH_sum = 0.0;
        mpData(mp).dH_sum_strain = 0.0;
        mpData(mp).Roud   = rho;
        mpData(mp).Saturated_YES   = 0;
        mpData(mp).Pore_Pressure = 0.0;
    end
end
end

% 2维按照固定一个维度对另一个维度进行排序
function LS1 = find_ChongFu_2D(boundary_MAT,nD,FX)
    boundary_MAT = sortrows(boundary_MAT,FX);
    LS = zeros(size(boundary_MAT,1),2);
    LS(:,1) = 1:size(boundary_MAT,1);
    LS(2:end,2) = boundary_MAT(2:end,FX)-boundary_MAT(1:end-1,FX);
    LS = LS(abs(LS(:,2))>1e-4,:);
    LS = [1,1;LS]; LS = [LS;size(boundary_MAT,1)+1,1];

    LS1 = zeros(2*(size(LS,1)-1),nD);
    for i = 1:size(LS,1)-1
        if FX ==1
            LS1(2*(i-1)+1,:) = [boundary_MAT(LS(i),1), max(boundary_MAT(LS(i):LS(i+1)-1,2))];
            LS1(2*(i-1)+2,:) = [boundary_MAT(LS(i),1), min(boundary_MAT(LS(i):LS(i+1)-1,2))];
        elseif FX ==2
            LS1(2*(i-1)+1,:) = [max(boundary_MAT(LS(i):LS(i+1)-1,1)),boundary_MAT(LS(i),2)];
            LS1(2*(i-1)+2,:) = [min(boundary_MAT(LS(i):LS(i+1)-1,1)),boundary_MAT(LS(i),2)];
        end
    end
end

% 3维按照固定一个维度对另一个维度进行排序
function LS1 = find_ChongFu_3D(boundary_MAT,nD,FX)
if FX==1
    FX1 = 2;
    FX2 = 3;
elseif FX==2
    FX1 = 1;
    FX2 = 3;
elseif FX==3
    FX1 = 1;
    FX2 = 2;
end
% 排序，cell分解
LS_Hang = 0;
boundary_MAT = sortrows(boundary_MAT,FX1);
[~, IA1, ~] = uniquetol(boundary_MAT(:,FX1),1e-3);
IA1(end+1) = size(boundary_MAT,1)+1;
CELL_X{length(IA1)-1,1}=[];
CELL_X1{length(IA1)-1,1}=[];
for i=1:length(IA1)-1
    CELL_X1{i} = boundary_MAT(IA1(i):IA1(i+1)-1,:);
    boundary_LS = sortrows(CELL_X1{i},FX2);

    [~, IA2, ~] = uniquetol(boundary_LS(:,FX2),1e-3);
    IA2(end+1) = size(boundary_LS,1)+1;
    CELL_X2{length(IA2)-1,1}=[];
    for j=1:length(IA2)-1
        CELL_X2{j} = boundary_LS(IA2(j):IA2(j+1)-1,:);
    end
    LS_Hang = LS_Hang + length(IA2) - 1;
    CELL_X{i}=CELL_X2;
end

% 组装 A
LS1 = zeros(LS_Hang*2,nD); AHA = 1;
for i = 1:size(CELL_X,1)
    for j = 1:size(CELL_X{i},1)
        [~,N]=max(CELL_X{i}{j}(:,FX));
        LS1(AHA,:) = CELL_X{i}{j}(N,:);
        [~,N]=min(CELL_X{i}{j}(:,FX));
        LS1(LS_Hang+AHA,:) = CELL_X{i}{j}(N,:);
        AHA = AHA+1;
    end
end


end
