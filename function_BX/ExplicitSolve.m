function [mpData, dU, V, a, drct] = ExplicitSolve(mpData, fint, fext, M, C, dt, V_old, bc_dofs, isFirstStep, fd)
%EXPLICITSOLVE 显式中心差分单步求解器 (含约束反力与物质点覆盖节点筛选)
%
%   输入:
%     fint        - 节点内力向量 [ndof x 1]
%     fext        - 节点外力向量 [ndof x 1] (不含反力)
%     M           - 集中质量向量 [ndof x 1] 或对角矩阵 [ndof x ndof]
%     C           - 阻尼矩阵 [ndof x ndof] (支持稀疏/稠密非对角)
%     dt          - 时间步长
%     V_old       - 上一步速度 [ndof x 1] (首步 v_0, 常规步 v_{n-1/2})
%     bc_dofs     - (可选) 约束自由度索引 [N1 x 2]，第一列为自由度编号
%     isFirstStep - (可选) true=首步, false=常规步 (默认false)
%     fd          - (可选) 有物质点的自由自由度索引 [1 x N2] 或 [N2 x 1]
%
%   输出:
%     dU   - 位移增量 [ndof x 1], 主程序执行 U = U + dU
%     V    - 更新后速度 [ndof x 1]
%     a    - 当前步加速度 [ndof x 1]
%     drct - 约束反力 [ndof x 1], 非约束自由度为 0

%% 质量矩阵对角化处理
if ~isvector(M)
    M = full(diag(M));
end
M = M(:);
nDoF = length(M);

%% 阻尼力
F_damp = C * V_old;

%% 残余力 (不平衡力)
Fres = fext - fint - F_damp;
drct = zeros(nDoF, 1);

%% 约束处理: 反力计算与残余力修正
if nargin >= 7 && ~isempty(bc_dofs)
    bc_idx = bc_dofs(:,1);              % 提取第一列自由度索引 [N1 x 1]
    drct(bc_idx) = -Fres(bc_idx);       % 反力 = -残余力
    fext_eff = fext + drct;
    Fres_eff = fext_eff - fint - F_damp;
else
    Fres_eff = Fres;
end

%% 加速度
a = Fres_eff ./ M;

%% 速度更新 (中心差分)
if nargin >= 8 && ~isempty(isFirstStep) && isFirstStep
    V = V_old + 0.5 * dt * a;           % 首步: v_{1/2} = v_0 + 0.5*dt*a_0
else
    V = V_old + dt * a;                 % 常规: v_{n+1/2} = v_{n-1/2} + dt*a_n
end

%% 位移增量
dU = dt * (V_old + V)/2;

%% 强制约束边界置零
if nargin >= 7 && ~isempty(bc_dofs)
    a(bc_idx)  = 0;
    V(bc_idx)  = 0;
    dU(bc_idx) = 0;
end

%% 非 fd 节点全部置零 (无物质点覆盖的背景网格节点冻结)
if nargin >= 9 && ~isempty(fd)
    fd_idx = fd(:);                     % 统一转为列向量 [N2 x 1]
    non_fd = true(nDoF, 1);
    non_fd(fd_idx) = false;             % fd 自由度保留，其余标记为 true
    a(non_fd)  = 0;
    V(non_fd)  = 0;
    dU(non_fd) = 0;
end

%% 更新F
nmp   = length(mpData);                                                     % number of material points
ddF   = zeros(nmp,3,3);                                                     % derivative of duvw wrt. spatial position
nD = length(mpData(1).mpC);

if nD==1                                                                    % 1D case
    fPos=1;                                                                 % deformation gradient positions                                                             % Cauchy stress components for internal force
elseif nD==2                                                                % 2D case (plane strain & stress)
    fPos=[1 5 4 2];
else                                                                        % 3D case
    fPos=[1 5 9 4 2 8 6 3 7];
end

nIN_all = [mpData.nIN];
nn_all = [mpData.nn]';
ed  = repmat((nIN_all-1)*nD,nD,1)+repmat((1:nD).',1,sum(nn_all));           % degrees of freedom of nodes (matrix form)
ed  = reshape(ed,1,sum(nn_all)*nD);                                         % degrees of freedom of nodes (vector form)
dNx_all = [mpData.dSvp];
nn_sizeall = [mpData.nn_size]';

if nD==1                                                                    % 1D case
    G_all = dNx_all;                                                        % strain-displacement matrix
elseif nD==2                                                                % 2D case (plane strain & stress)
    G_all = zeros(4,nD*sum(nn_all));                                        % zero the strain-disp matrix (2D)
    G_all([1 3],1:nD:end)=dNx_all;                                          % strain-displacement matrix
    G_all([4 2],2:nD:end)=dNx_all;
else                                                                        % 3D case
    G_all=zeros(9,nD*sum(nn_all));                                          % zero the strain-disp matrix (3D)
    G_all([1 4 9],1:nD:end)=dNx_all;                                        % strain-displacement matrix
    G_all([5 2 6],2:nD:end)=dNx_all;
    G_all([8 7 3],3:nD:end)=dNx_all;
end

% 根据节点位移计算节点应变-一次性
% 1. ddF计算
% 通过索引一次性提取所有需要的 uvw 值
% 结果：values(j) = uvw( ed_all(j) )
values = dU(ed);  % 189796 x 1
% 逐元素加权（广播机制）
% G_scaled(:,j) = G_all(:,j) * values(j)
G_scaled = G_all .* values';  % 4 x 189796
% 列方向累积求和（前缀和）
Prefix = cumsum(G_scaled, 2);  % 4 x 189796
% 通过索引相减提取各段结果（向量化，无循环）
s = (nn_sizeall(:,1)-1)*nD + 1;
e = nn_sizeall(:,2)*nD;
% 初始化结果矩阵 4 x 23530
Result = Prefix(:, e);
% 处理起始位置大于1的段：减去前一个前缀和
% 利用逻辑掩码处理 s=1 的边界情况
mask = (s > 1);
if any(mask)
    s_prev = s - 1;
    Result(:, mask) = Result(:, mask) - Prefix(:, s_prev(mask));
end

ddF(:,fPos) = Result';

% 2. dF 计算（保持不变）
I = zeros(nmp, 3, 3);
for k = 1:3, I(:,k,k) = 1; end

dF = I + ddF;

% 3. F 计算（若需要，补充缺失的乘法）
Fn_p = reshape([mpData.Fn], 3, 3, nmp);
dF_p = permute(dF, [2, 3, 1]);
F = pagemtimes(dF_p, Fn_p);

for mp=1:nmp
    mpData(mp).F    = F(:,:,mp);                                            % store deformation gradient
end

end