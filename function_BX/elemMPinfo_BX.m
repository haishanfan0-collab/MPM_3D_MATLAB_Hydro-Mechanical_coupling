function [mesh,mpData,VA_old,bounday_imfo] = elemMPinfo_BX(mesh,mpData,method)

%Determine the basis functions for material points 
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to determine the basis functions and spatial derivatives of each
% material point.  The function works for regular background meshes with
% both the standard and generalised interpolation material point methods.
% The function also determines, and stores, the elements associated with
% the material point and a unique list of nodes that the material point
% influences.  The number of stiffness matrix entries for each material
% point is determined and stored. 
%
%--------------------------------------------------------------------------
% [fbdy,mpData] = ELEMMPINFO(mesh,mpData)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - mesh structured array. Function requires: 
%           - coord : coordinates of the grid nodes (nodes,nD)
%           - etpl  : element topology (nels,nen) 
%           - h     : background mesh size (nD,1)
% mpData - material point structured array.  Function requires:
%           - mpC   : material point coordinates
%--------------------------------------------------------------------------
% Ouput(s);
% mesh   - mesh structured array. Function modifies:
%           - eInA  : elements in the analysis 
% mpData - material point structured array. Function modifies:
%           - nIN   : nodes linked to the material point
%           - eIN   : element associated with the material point
%           - Svp   : basis functions for the material point
%           - dSvp  : basis function derivatives (at start of lstp)
%           - nSMe  : number stiffness matrix entries for the MP
%--------------------------------------------------------------------------
% See also:
% ELEMFORMP         - find elements for material point
% NODESFORMP        - nodes associated with a material point 
% MPMBASIS          - MPM basis functions
%--------------------------------------------------------------------------
[nodes,nD] = size(mesh.coord); 
nmp      = size(mpData,2);                                                  % number of material points & dimensions
nels = size(mesh.etpl,1);                                                   % numerb of elements in mesh
mpC  = reshape([mpData.mpC],nD,nmp).';                                      % all material point coordinates (nmp,nD)
lp   = reshape([mpData.lp] ,nD,nmp).';                                      % all domain lengths
eInA = zeros(nels,1);                                                       % zero elements taking part in the analysis
nDoF = nodes*nD; 

dudt_old  = zeros(nDoF,1);                                                  % zero internal dudt_old vector
du2dt2_old  = zeros(nDoF,1);                                                % zero internal du2dt2_old vector
mesh.HVA_old.dHdt = zeros(nodes,1);
mesh.HVA_old.H = zeros(nodes,1);
Name_mp = fieldnames(mpData);

FHS_JL{nodes,3}  = [];

LS=[];
FHS_AHA=0;

%% 网格区域信息
bounday_imfo.mesh_in=zeros(nels,1);
bounday_imfo.mesh_boundry=zeros(nels,1);
bounday_imfo.mesh_out=(1:nels)';
bounday_imfo.node_boundry=zeros(nDoF,nD+1);
bounday_imfo.node_boundry_mesh={};

FHS_JS=1;
etpl_S2 = size(mesh.etpl,2);
for mp = 1:nmp
    try
        eIN  = elemForMP_BX(mesh,mpC(mp,:),lp(mp,:),nD);                              % elements connected to the material point
        nIN = unique(sort(reshape(mesh.etpl(eIN,:),length(eIN)*etpl_S2,1),1))';
        nn   = length(nIN);                                                     % number of nodes influencing the MP
        
        [Svp,dSvp] = MPMbasis_BX(mesh,mpData(mp),nIN);                          % basis function and spatial derivatives                                      % basis function derivatives

        FHS_JL(nIN) = cellfun(@(x) [x, mp-FHS_AHA], FHS_JL(nIN), 'UniformOutput', false);
        
        % 进行网格信息采集
        if mpData(mp).mpType ~=3
            bounday_imfo.mesh_in(eIN)=eIN;
            bounday_imfo.mesh_out(eIN)=0;
        else
            bounday_imfo.mesh_boundry(eIN)=eIN;
            bounday_imfo.mesh_out(eIN)=0;
            bounday_imfo.mesh_in(eIN)=0;
            bounday_imfo.node_boundry(FHS_JS,:) = [mp-FHS_AHA,mpData(mp).mpC];
            bounday_imfo.node_boundry_mesh{end+1,1} = eIN;
            FHS_JS = FHS_JS+1;
        end

        mpData(mp).nIN  = nIN;                                                  % nodes associated with material point
        mpData(mp).eIN  = eIN;                                                  % elements associated with material point
        mpData(mp).Svp  = Svp/sum(Svp);                                                 % basis functions
        mpData(mp).dSvp = dSvp;                                                 % basis function derivatives
        mpData(mp).nSMe = (nn*nD)^2;                                            % number stiffness matrix components
        mpData(mp).nn = nn;
        if mp==1
            mpData(mp).nn_size = [1;nn];
        else
            mpData(mp).nn_size = [mpData(mp-1).nn_size(end)+1;mpData(mp-1).nn_size(end)+nn];
        end
        eInA(eIN) = 1;                                                          % identify elements in the analysis

    catch
        LS(end+1)=mp;
        FHS_AHA=FHS_AHA+1;
    end

end
mesh.eInA = eInA;                                                           % store eInA to mesh structured array
% 删除异常物质点 /保险
mpData(LS)=[];
mesh.Contact_Force(LS,:) = [];

%% 动力模式下进行速度和加速度插值
% 第一步，获取平面散点
nmp      = size(mpData,2);
mpC = reshape([mpData.mpC]',nD,nmp)';

mp_number = 4;
CZ_Pos = -1e3 * ones(nodes,nD);
node_all = zeros(nodes,1);
if strcmp(method, 'dynamic')
    for node = 1:nodes
        if ~isempty(FHS_JL{node,1})
            nn   = length(FHS_JL{node,1});

            mp = FHS_JL{node,1};                                                      % current node
            [Svp,~] = MPMbasis_BX(mesh,mpData(mp),node);                            % basis function and spatial derivatives

            FHS_JL{node,2}=Svp;
            FHS_JL{node,3}=sum(Svp);

            dudt_old(nD*(node-1)+1:nD*node) = sum([Svp',Svp'].*reshape([mpData(FHS_JL{node,1}).du],nD,nn).')/max(1.0,FHS_JL{node,3});
            du2dt2_old(nD*(node-1)+1:nD*node) = sum([Svp',Svp'].*reshape([mpData(FHS_JL{node,1}).ddu],nD,nn).')/max(1.0,FHS_JL{node,3});
            
            if sum(strcmp(Name_mp,'Flow_Type'))>0
                if FHS_JL{node,3}>=mp_number^nD*0.99

                    mesh.HVA_old.dHdt(node) = Svp*[mpData(FHS_JL{node,1}).dH]'/max(1.0,FHS_JL{node,3});
                    % 总水头 = 初始吸力水头+位置水头+水头累积
                    mesh.HVA_old.H(node) = Svp*([mpData(FHS_JL{node,1}).H0] + ...
                        [mpData(FHS_JL{node,1}).dH_sum] + ...
                        [mpData(FHS_JL{node,1}).dH_sum_strain])'/max(1.0,FHS_JL{node,3}) + mesh.coord(node,2);
                else
                    CZ_Pos(node,:) = mesh.coord(node,:);
                    node_all(node) = node;
                end
            end
        end
    end
    CZ_Pos(CZ_Pos(:,1)<0,:) = [];
    node_all(node_all==0) = [];
    mesh.HVA_old.dHdt(node_all) = interpMPtoNode(CZ_Pos, mpC, [mpData.dH]', mesh.h(1));
    % 总水头 = 初始吸力水头+位置水头+水头累积
    mesh.HVA_old.H(node_all) = interpMPtoNode(CZ_Pos, mpC, [mpData.H0]'+[mpData.dH_sum]'+[mpData.dH_sum_strain]', mesh.h(1)) + CZ_Pos(:,2);
end


%% 边界
dudt_old(mesh.bc(:,1))=0;
du2dt2_old(mesh.bc(:,1))=0;

%% 进行网格区域识别
bounday_imfo.mesh_in(bounday_imfo.mesh_in(:,1)==0,:)=[];
bounday_imfo.mesh_boundry(bounday_imfo.mesh_boundry(:,1)==0,:)=[];
bounday_imfo.mesh_out(bounday_imfo.mesh_out(:,1)==0,:)=[];
bounday_imfo.node_boundry(bounday_imfo.node_boundry(:,1)==0,:)=[];

VA_old.dudt = dudt_old;
VA_old.du2dt2 = du2dt2_old;

%% 渗流识别-更新
nmp = size(mpData,2); LS = (1:nmp)';
if sum(strcmp(Name_mp,'Flow_P'))>0
    mpData_Select = mpData([mpData.Flow_P]~=0);
    LS = LS([mpData.Flow_P]~=0);
    if ~isempty(mpData_Select)
        Name = unique([mpData_Select.Flow_Position]);
        for i = 1:length(Name)
            mpData_SelectLS = mpData_Select(strcmp([mpData_Select.Flow_Position],Name(i)));
            LS_LS = [LS(strcmp([mpData_Select.Flow_Position],Name(i))),...
                (reshape([mpData_SelectLS.mpC]',2,size(mpData_SelectLS,2)))',...
                [mpData_SelectLS.Flow_SIZE]'];
            LS_LS = sortrows(LS_LS,2);
            for j = 1:size(LS_LS,1)
                if j == 1
                    mpData(LS_LS(j,1)).Flow_P = LS_LS(j,end)*abs(LS_LS(j+1,2)-LS_LS(j,2));
                elseif j == size(LS_LS,1)
                    mpData(LS_LS(j,1)).Flow_P = LS_LS(j,end)*abs(LS_LS(j,2)-LS_LS(j-1,2));
                else
                    mpData(LS_LS(j,1)).Flow_P = LS_LS(j,end)*abs(LS_LS(j+1,2)-LS_LS(j-1,2))*0.50;
                end
            end
        end
    end
end
end




function B = interpMPtoNode(CZ_Pos, mpC, A, cell_size, varargin)
% INTERPBOUNDARYNODEFAST 极速版边界插值（向量化优化）
%   优化策略：
%   1. 分块向量化距离计算（避免逐节点循环）
%   2. 移除昂贵的cond()计算，改用rcond()或简单行列式检查
%   3. 简化角部检测逻辑（点少才检测，点多直接IDW）
%   4. 预分配所有数组，避免动态扩容

    %% 参数设置
    p = inputParser;
    addParameter(p, 'FieldType', 'none');
    addParameter(p, 'MinPoints', 4);
    addParameter(p, 'AngleThresh', pi);
    addParameter(p, 'ChunkSize', 1000);  % 分块大小，根据内存调整
    parse(p, varargin{:});
    
    fieldType = p.Results.FieldType;
    minPts = p.Results.MinPoints;
    angleThresh = p.Results.AngleThresh;
    chunkSize = p.Results.ChunkSize;
    
    [nNodes, ~] = size(CZ_Pos);
    % nMP = size(mpC, 1);
    B = zeros(nNodes, 1);
    
    % 搜索半径
    searchRad = 2.5 * cell_size;
    searchRad2 = searchRad^2;
    h2 = (0.8 * cell_size)^2;  % 高斯权重参数
    
    %% 分块向量化处理（避免内存爆炸）
    nChunks = ceil(nNodes / chunkSize);
    
    for iChunk = 1:nChunks
        idxStart = (iChunk-1)*chunkSize + 1;
        idxEnd = min(iChunk*chunkSize, nNodes);
        nCurrent = idxEnd - idxStart + 1;
        
        nodesChunk = CZ_Pos(idxStart:idxEnd, :);
        
        % === 核心优化：向量化距离矩阵计算 ===
        % 结果: dist2Chunk(i,j) = 第i个节点到第j个物质点的距离平方
        % 使用广播机制，避免显式循环
        dx = nodesChunk(:,1) - mpC(:,1)';  % nCurrent × nMP
        dy = nodesChunk(:,2) - mpC(:,2)';
        dist2Chunk = dx.^2 + dy.^2;
        
        % 对每个节点快速处理（此时内层循环开销极小）
        for i = 1:nCurrent
            globalIdx = idxStart + i - 1;
            dist2Row = dist2Chunk(i, :);
            
            % 快速筛选（逻辑索引，比find快）
            nearbyMask = dist2Row < searchRad2;
            nNearby = sum(nearbyMask);
            
            if nNearby < 3
                % 极少点：最近邻
                [~, nearest] = min(dist2Row);
                B(globalIdx) = A(nearest);
                continue;
            end
            
            % 提取局部数据（只提取需要的列）
            idxNearby = find(nearbyMask);
            x_loc = mpC(idxNearby, 1) - nodesChunk(i, 1);
            y_loc = mpC(idxNearby, 2) - nodesChunk(i, 2);
            A_loc = A(idxNearby);
            d_loc = sqrt(dist2Row(idxNearby))';
            
            % === 快速几何分类 ===
            % 优化：只有点很少（<8）时才做角部检测（计算atan2较贵）
            if nNearby <= 8
                angles = atan2(y_loc, x_loc);  % 注意：这里应该是y_loc, x_loc
                covAngle = fastCoverageAngle(angles);
                isCorner = (covAngle < angleThresh);
            else
                isCorner = false;  % 点多默认为内部或一般边界
            end
            
            % === 快速插值分支 ===
            if isCorner
                % 角部：最近4点IDW + 硬截断
                [~, sortIdx] = mink(d_loc, min(4, nNearby));  % mink比sort快
                w = 1./(d_loc(sortIdx).^2 + 0.01*cell_size^2);
                val = (w' * A_loc(sortIdx)) / sum(w);
                % 值域截断（防止过冲）
                localMin = min(A_loc(sortIdx));
                localMax = max(A_loc(sortIdx));
                val = max(localMin, min(localMax, val));
                
            elseif nNearby >= minPts
                % 一般边界：快速MLS（简化权重计算）
                val = fastMLS(x_loc, y_loc, A_loc, d_loc, cell_size);
            else
                % 内部或点不足：高斯IDW
                w = exp(-d_loc.^2 / h2);
                val = (w' * A_loc) / sum(w);
            end
            
            B(globalIdx) = val;
        end
    end
    
    %% 物理约束（可选，向量化处理）
    if ~strcmp(fieldType, 'none')
        B = applyConstraintVectorized(B, fieldType);
    end
end

%% 快速覆盖角计算（不排序所有角度，只找最大间隙）
function covAngle = fastCoverageAngle(angles)
    if length(angles) < 3
        covAngle = 2*pi; 
        return;
    end
    % 快速检查：如果点分布在四个象限， coverage一定大
    % 精确计算：只找最大间隙
    s = sort(angles);
    gaps = diff([s; s(1)+2*pi]);  % 环形差分
    maxGap = max(gaps);
    covAngle = 2*pi - maxGap;
end

%% 极速MLS（避免cond，使用rcond或简单检查）
function val = fastMLS(x, y, A, dist, cell_size)
    n = length(x);
    
    % 设计矩阵
    P = [ones(n,1), x, y];
    
    % 权重（高斯）
    h = cell_size * 0.6;
    w = exp(-(dist/h).^2);
    
    % 加权正规方程 P'WP * coeff = P'WA
    % 直接计算，不构建稀疏W
    WP = P .* w;  % 广播乘法，比diag(w)*P快
    lhs = P' * WP;
    rhs = P' * (w .* A);
    
    % 快速病态检查（比cond快100倍）
    % 方法1：检查对角线
    diagLhs = abs(diag(lhs));
    if min(diagLhs) < 1e-12 * max(diagLhs)
        % 病态：退化IDW
        val = (w' * A) / sum(w);
        return;
    end
    
    % 方法2：尝试求解，捕获警告
    try
        coeff = lhs \ rhs;
        val = coeff(1);
    catch
        % 求解失败，退化IDW
        val = (w' * A) / sum(w);
    end
end

%% 向量化约束处理
function B = applyConstraintVectorized(B, fieldType)
    switch lower(fieldType)
        case 'suction'
            B = max(0, B);  % 非负
        case 'saturation'
            B = max(0, min(1, B));
        case 'displacement'
            % 位移限制在全局极值（软约束）
            % 注意：这里假设B整体合理，只做硬边界
            B = max(-1e6, min(1e6, B));
    end
end