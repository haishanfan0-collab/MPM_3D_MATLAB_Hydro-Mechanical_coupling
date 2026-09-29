function [mesh,mpData,VA_old,bounday_imfo] = elemMPinfo(mesh,mpData,Cal_par)

% 确定物质点形函数及渗流信息
%--------------------------------------------------------------------------
% 作者: FAN Haishan（改动范围超过60%） 
% 原始代码作者：William Coombs
% 【改动范围】 新增速度、加速度、水头插值；新增渗流流量插值，支持动力、渗流计算
% 日期: 28/04/2026

% 描述:
% 该子函数分为4个步骤。（1）确定物质点所属背景网格单元，读取对应的节点序列，计算形函数及
% 形函数导数。（2）识别背景网格所属类别，对背景网格进行分类，为接触碰撞算法做准备。（3）
% 进行物质点→背景网格插值，插值信息包括：速度、加速度、水头。（4）完成渗流信息计算，根据
% 流量和渗流边界物质点，进行物质点流量计算
%
%--------------------------------------------------------------------------
% [mesh,mpData,VA_old,bounday_imfo] = elemMPinfo(mesh,mpData,method)
%--------------------------------------------------------------------------
% 输入:
% mesh   - 网格结构体参数. 包含: 
%           - coord : 节点坐标 (nodes,nD)
%           - etpl  : 各单元节点信息 (nels,nen) 
%           - h     : 背景网格尺寸 (nD,1)
% mpData - 物质点结构体参数.  包含:
%           - mpC   : 物质点坐标
% method - 求解方法，分为"dynamic"和"static"
%--------------------------------------------------------------------------
% 输出:
% mesh   - 网格结构体参数. 包含:
%           - eInA  : 激活的背景网格节点信息 
%           - HVA_old: 节点水头信息（由于后续开发，子函数输出变量原因，直接嵌入mesh
% mpData - 物质点结构体参数.  包含:
%           - nIN   : 物质点关联背景网格节点
%           - eIN   : 物质点关联背景网格
%           - Svp   : 物质点形函数，用于相互插值计算
%           - dSvp  : 物质点形函数导数，用于应变、水头梯度计算
%           - nSMe  : 物质点关联单元的节点刚度矩阵尺寸
% VA_old - 速度加速度插值信息，结构体参数，包含 dudt du2dt2两个参数
% bounday_imfo - 背景网格所属信息，按照背景网格包含物质点分为：边界网格、内部网格和外
%                部网格，用于后续接触碰撞算法的接触对识别
%--------------------------------------------------------------------------
% 该函数包含子函数:
% elemForMP         - 物质点所属单元搜索
% nodesForMP        - 背景网格单元节点搜索
% MPMbasis          - MPM 形函数
% interpMPtoNode    - 插值代码，完成物质点 → 背景网格投影计算
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

for mp = 1:nmp
    try
        eIN  = elemForMP(mesh,mpC(mp,:),lp(mp,:));                              % elements connected to the material point
        nIN  = nodesForMP(mesh.etpl,eIN).';                                     % unique list of nodes associated with elements
        nn   = length(nIN);                                                     % number of nodes influencing the MP
        Svp  = zeros(1,nn);                                                     % zero basis functions
        dSvp = zeros(nD,nn);                                                    % zero basis function derivatives
        for i = 1:nn
            node = nIN(i);                                                      % current node
            [S,dS] = MPMbasis(mesh,mpData(mp),node);                            % basis function and spatial derivatives
            Svp(i) = Svp(i) + S;                                                % basis functions for all nodes
            dSvp(:,i) = dSvp(:,i) + dS;                                         % basis function derivatives
            FHS_JL{node,1}(end+1)=mp-FHS_AHA;
        end
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

% du1 = [mpC,du(:,1)];
% du2 = [mpC,du(:,2)];
% ddu1 = [mpC,ddu(:,1)];
% ddu2 = [mpC,ddu(:,2)];

% 核心：选cubic核（适配局部突变+全局光滑）
% rbf_du1 = scatteredInterpolant(du1(:,1), du1(:,2), du1(:,3),'natural','boundary'); 
% rbf_du2 = scatteredInterpolant(du2(:,1), du2(:,2), du2(:,3),'natural','boundary'); 
% rbf_ddu1 = scatteredInterpolant(ddu1(:,1), ddu1(:,2), ddu1(:,3),'natural','boundary'); 
% rbf_ddu2 = scatteredInterpolant(ddu2(:,1), ddu2(:,2), ddu2(:,3),'natural','boundary'); 

% if sum(strcmp(Name_mp,'Flow_Type'))>0
    % H0 = [mpC,[mpData.H0]'];
    % dH = [mpC,[mpData.dH]'];
    % dH_sum = [mpC,[mpData.dH_sum]'];
    % dH_sum_strain = [mpC,[mpData.dH_sum_strain]'];
    % rbf_H0 = scatteredInterpolant(H0(:,1), H0(:,2), H0(:,3),'natural','boundary');
    % rbf_dH = scatteredInterpolant(dH(:,1), dH(:,2), dH(:,3),'natural','boundary');
    % rbf_dH_sum = scatteredInterpolant(dH_sum(:,1), dH_sum(:,2), dH_sum(:,3),'natural','boundary');
    % rbf_dH_sum_strain = scatteredInterpolant(dH_sum_strain(:,1), dH_sum_strain(:,2), dH_sum_strain(:,3),'natural','boundary');
% end

% if strcmp(method, 'dynamic')
%     for node = 1:nodes
%         if ~isempty(FHS_JL{node,1})
%             CZ_Pos = mesh.coord(node,:);
%             dudt_old(nD*(node-1)+1) = rbf_du1(CZ_Pos(1),CZ_Pos(2));
%             dudt_old(nD*(node-1)+2) = rbf_du2(CZ_Pos(1),CZ_Pos(2));
%             du2dt2_old(nD*(node-1)+1) = rbf_ddu1(CZ_Pos(1),CZ_Pos(2));
%             du2dt2_old(nD*(node-1)+2) = rbf_ddu2(CZ_Pos(1),CZ_Pos(2));
%             if sum(strcmp(Name_mp,'Flow_Type'))>0
%                 mesh.HVA_old.dHdt(node) = rbf_dH(CZ_Pos(1),CZ_Pos(2));
%                 % 总水头 = 初始吸力水头+位置水头+水头累积
%                 mesh.HVA_old.H(node) = rbf_H0(CZ_Pos(1),CZ_Pos(2)) + ...
%                                        rbf_dH_sum(CZ_Pos(1),CZ_Pos(2)) + ...
%                                        rbf_dH_sum_strain(CZ_Pos(1),CZ_Pos(2));
%             end
%         end
%     end
% end
mp_number = 4;
CZ_Pos = -1e100 * ones(nodes,nD);
node_all = zeros(nodes,1);

for node = 1:nodes
    if ~isempty(FHS_JL{node,1})
        nn   = length(FHS_JL{node,1});
        Svp  = zeros(1,nn);
        for i = 1:nn
            mp = FHS_JL{node,1}(i);                                                      % current node
            [S,~] = MPMbasis(mesh,mpData(mp),node);                            % basis function and spatial derivatives
            Svp(i) = Svp(i) + S;                                                % basis functions for all nodes
        end
        FHS_JL{node,2}=Svp;
        FHS_JL{node,3}=sum(Svp);

        dudt_old(nD*(node-1)+1:nD*node) = sum(repmat(Svp', 1, nD).*reshape([mpData(FHS_JL{node,1}).du],nD,nn).')/max(1.0,FHS_JL{node,3});
        du2dt2_old(nD*(node-1)+1:nD*node) = sum(repmat(Svp', 1, nD).*reshape([mpData(FHS_JL{node,1}).ddu],nD,nn).')/max(1.0,FHS_JL{node,3});

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
                % dudt_old(nD*(node-1)+1) = interpMPtoNode(CZ_Pos, mpC, du(:,1), mesh.h(1));
                % dudt_old(nD*(node-1)+2) = interpMPtoNode(CZ_Pos, mpC, du(:,2), mesh.h(1));
                % du2dt2_old(nD*(node-1)+1) = interpMPtoNode(CZ_Pos, mpC, ddu(:,1), mesh.h(1));
                % du2dt2_old(nD*(node-1)+2) = interpMPtoNode(CZ_Pos, mpC, ddu(:,2), mesh.h(1));
                % mesh.HVA_old.dHdt(node) = interpMPtoNode(CZ_Pos, mpC, [mpData.dH]', mesh.h(1));
                % 总水头 = 初始吸力水头+位置水头+水头累积
                % mesh.HVA_old.H(node) = interpMPtoNode(CZ_Pos, mpC, [mpData.H0]', mesh.h(1)) + ...
                %     interpMPtoNode(CZ_Pos, mpC, [mpData.dH_sum]', mesh.h(1)) + ...
                %     interpMPtoNode(CZ_Pos, mpC, [mpData.dH_sum_strain]', mesh.h(1)) + mesh.coord(node,2);
            end
        end
    end
end

CZ_Pos(CZ_Pos(:,1)<-1e99,:) = [];
node_all(node_all==0) = [];
mesh.HVA_old.dHdt(node_all) = interpMPtoNode(CZ_Pos, mpC, [mpData.dH]', mesh.h(1));
% 总水头 = 初始吸力水头+位置水头+水头累积
mesh.HVA_old.H(node_all) = interpMPtoNode(CZ_Pos, mpC, [mpData.H0]'+[mpData.dH_sum]'+[mpData.dH_sum_strain]', mesh.h(1)) + CZ_Pos(:,2);



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

Q = interpFlow(max(0.0,Cal_par.Calculate_time - Cal_par.Gravity_Time));
if abs(Q)<1e-50
    Q = 1e-50;
end

if sum(strcmp(Name_mp,'Flow_P'))>0
    if nD == 2
        mpData_Select = mpData([mpData.Flow_P]~=0);
        LS = LS([mpData.Flow_P]~=0);
        if ~isempty(mpData_Select)
            Name = unique([mpData_Select.Flow_Position]);
            for i = 1:length(Name)
                mpData_SelectLS = mpData_Select(strcmp([mpData_Select.Flow_Position],Name(i)));
                LS_LS = [LS(strcmp([mpData_Select.Flow_Position],Name(i))),...
                    (reshape([mpData_SelectLS.mpC]',2,size(mpData_SelectLS,2)))'];
                LS_LS = sortrows(LS_LS,2);
                for j = 1:size(LS_LS,1)
                    if j == 1
                        mpData(LS_LS(j,1)).Flow_P = Q*abs(LS_LS(j+1,2)-LS_LS(j,2));
                    elseif j == size(LS_LS,1)
                        mpData(LS_LS(j,1)).Flow_P = Q*abs(LS_LS(j,2)-LS_LS(j-1,2));
                    else
                        mpData(LS_LS(j,1)).Flow_P = Q*abs(LS_LS(j+1,2)-LS_LS(j-1,2))*0.50;
                    end
                end
            end
        end
    elseif nD==3
        % 三维流量更新代码（后加，不影响平面应变代码）
        mpData_mpC = [mpData([mpData.mpType] == 3).mpC];
        mpData_mpC = reshape(mpData_mpC,nD,length(mpData_mpC)/nD)';
        for mp = 1:nmp
            % 邻域半径（半宽）
            dx = 0.5 * lp(mp,1);
            dz = 0.5 * lp(mp,3);

            % 找出 X、Y 均在邻域范围内的点（含自身）
            inRange = abs(mpC(:,1) - mpC(mp,1)) <= dx & ...
                abs(mpC(:,3) - mpC(mp,3)) <= dz;

            % 判断当前点的 Z 是否为邻域内最大值（含自身用 >=）
            if mpC(mp,2) >= max(mpC(inRange,2))
                mpData(mp).mpType = 3;
                mpData(mp).Flow_Position = 'Slope';
            else
                mpData(mp).mpType = 2;
                mpData(mp).Flow_Position = [];
            end

            if mpData(mp).mpType == 3
                W_area = lp(mp,1)*lp(mp,3);
                mpData(mp).Flow_P = Q*W_area;
                mpData(mp).Flow_SIZE = Q;
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
    addParameter(p, 'ChunkSize', 500);  % 分块大小，根据内存调整
    parse(p, varargin{:});
    
    fieldType = p.Results.FieldType;
    minPts = p.Results.MinPoints;
    angleThresh = p.Results.AngleThresh;
    chunkSize = p.Results.ChunkSize;
    
    [nNodes, ~] = size(CZ_Pos);
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
                angles = atan2(y_loc, y_loc);  % 注意：这里应该是y_loc, x_loc
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

function Q = interpFlow(time)
% 【流量时程线性插值】
% 输入: Flow_SIZE — N×2 矩阵，第一列为时间，第二列为流量
%       time      — 当前时刻（浮点数）
% 输出: Q         — 插值得到的当前流量（标量或向量，取决于time维度）
%
% 处理逻辑:
%   1) time 超出范围时，外推取端点值（不报错）
%   2) time 恰好匹配节点时，直接返回对应流量
%   3) 多 time 输入时，支持向量化输出
    A = load('Flow_SIZE');
    % 提取时间和流量列
    T = A.Flow_SIZE(:,1);
    Qv = A.Flow_SIZE(:,2);

    % 确保时间单调递增（如果不是，先排序）
    if ~issorted(T)
        [T, idx] = sort(T);
        Qv = Qv(idx);
    end

    % 处理标量或向量输入
    time = time(:);  % 统一为列向量

    % 向量化插值
    Q = interp1(T, Qv, time, 'linear', 'extrap');

end