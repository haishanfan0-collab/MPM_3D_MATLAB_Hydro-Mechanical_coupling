function [mesh,mpData,HVA,HVA_old,oobf_Flow,fd_Flow,H,frct_Flow,Cal_par] = ...
    AMPLE_Solve_Flow_Main_BX(mesh,mpData,HVA_old,...  % 物质点/网格参数
    oobf_Flow,fd_Flow,...                          % 外力数据
    H,frct_Flow,...                                % 核心出装
    Cal_par)                                       % 计算参数

% 【核心-渗流场迭代求解代码】，内置渗流 Newton-Raphson 迭代控制、水头梯度更新、线性
% 方程组求解、渗流内力与刚度矩阵组装、收敛判断全过程
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行物质点法渗流场迭代求解，涵盖物质点水头梯度计算、整体渗流刚度/质量矩阵线性求解、
% 节点水头增量与反力更新、渗流不平衡力组装及 Newton-Raphson 收敛判断全过程，为湿-力
% 耦合分析中渗流场计算核心主函数
%
%--------------------------------------------------------------------------
% [mesh, mpData, HVA, HVA_old, oobf_Flow, fd_Flow, H, frct_Flow, Cal_par] =
%    AMPLE_Solve_Flow_Main_BX(mesh, mpData, HVA_old,
%                             oobf_Flow, fd_Flow,
%                             H, frct_Flow,
%                             Cal_par)
%--------------------------------------------------------------------------
% 输入:
% mesh      - 网格信息结构体参数（含 coord, bc_Flow, Flow.Kt/Ct 等渗流矩阵）
% mpData    - 物质点信息结构体参数
% HVA_old   - 上一时间步节点水头、水头速率信息结构体
% oobf_Flow - 渗流节点不平衡力向量
% fd_Flow   - 渗流需要求解的节点自由度序列
% H         - 当前节点总水头向量
% frct_Flow - 渗流节点约束反力向量
% Cal_par   - 求解控制参数结构体（含 NRit_Flow, dt, fErr_Flow 等）
%--------------------------------------------------------------------------
% 输出:
% mesh      - 更新后的网格信息结构体参数
% mpData    - 更新后的物质点信息结构体参数（含水头梯度 epsn_Flow 等）
% HVA       - 当前步节点水头速率信息结构体
% HVA_old   - 上一时间步节点水头速率信息结构体
% oobf_Flow - 更新后的渗流不平衡力向量
% fd_Flow   - 渗流需要求解的节点自由度序列
% H         - 更新后的节点总水头向量
% frct_Flow - 更新后的渗流节点约束反力向量
% Cal_par   - 更新后的求解控制参数结构体（迭代计数、误差等）
%--------------------------------------------------------------------------
% 此子函数包含子函数
% model_Solve_Flow - 渗流线性方程组求解器，执行边界条件处理、等效刚度矩阵组装、
%                    水头增量与反力计算
% detMPs_Flow_BX   - 物质点渗流本构及内部流量/刚度矩阵计算（向量版）
% detExtFlaw       - 渗流外部源项（降雨/蒸发等）节点力计算
%--------------------------------------------------------------------------

[nodes,~] = size(mesh.coord);
%% 关键：获取节点deltaH / 考虑由于变形引起的水头变化，每一次的 dH/dx 和 dH/dy 都是现算的
if Cal_par.NRit_Flow==0
    nmp   = length(mpData);                                                     % number of material points
    for mp=1:nmp                                                                % material point loop
        %------------------- 通过H求解deltaH，类似于通过位移增量更新 ---------------------------------------------------
        ed = mpData(mp).nIN;        % 物质点关联的节点（类比位移求解的节点集）
        G = mpData(mp).dSvp;        % 形函数导数 ∇N (nD × nn)（类比位移求解的B矩阵）
        mpData(mp).epsn_Flow = G*HVA_old.H(ed);  % 根据总水头获取水头变化率
        % mpData(mp).epsn_Flow(2) = 1 + mpData(mp).epsn_Flow(2);
    end
end

%% 求解渗流
[dH,dQ,HVA.dHdt] = model_Solve_Flow(mesh.bc_Flow, mesh.Flow.Kt, mesh.Flow.Ct, oobf_Flow,...
    Cal_par.NRit_Flow, fd_Flow, Cal_par.dt,...
    H, HVA_old.dHdt);            % linear solver
dH(abs(dH/Cal_par.dt)>0.10) = 0.0;
Cal_par.NRit_Flow = Cal_par.NRit_Flow+1;                                           % increment the NR counter
H  = H+dH;                                                                         % update displacements
frct_Flow = frct_Flow+dQ;                                                          % update reaction forces
[Qint,mpData,mesh] = detMPs_Flow_BX(H,mpData,mesh);                                   % global stiffness & internal force

Q_fext = detExtFlaw(nodes,mpData);
Q_FORCE = Q_fext+frct_Flow;

fint_all = (Qint+mesh.Flow.Ct*HVA.dHdt);
%
% check = [mesh.coord,Q_FORCE,Qint,fint_all,HVA.dHdt,H];
% check_Down=check(check(:,2)~=0,:);
% check_Left=check(check(:,1)~=0,:);
Cal_par.fErr_Flow = norm(Q_FORCE-fint_all)/norm(Q_FORCE+eps);                      % normalised oobf error
% disp(num2str( Cal_par.fErr_Flow ))
if Cal_par.NRit_Flow<5
    oobf_Flow = Q_FORCE-fint_all;
else
    oobf_Flow = oobf_Flow*0.1+(Q_FORCE-fint_all)*0.9;
end

% fprintf(1,'%s %2i %s: %8.3e, %s: %8.3e, %s: %8.3e \n','  iteration',Cal_par.NRit_Flow,' NR error',Cal_par.fErr_Flow,'max(oobf)',norm(abs(Q_FORCE-fint_all)), 'max(H)', norm(abs(H)));   % text output to screen (NR error)
end


%% 求解主程序
function [dH,dQ,dHdt] = model_Solve_Flow(bc, Kt, Ct, oobf, NRit, fd, dt, H, dHdt_old)
% 【渗流线性方程组求解器】，内置等效刚度矩阵组装、边界条件处理、水头增量与反力计算
%--------------------------------------------------------------------------
% 作者: FAN Haishan（基于 William Coombs 原始框架重构）
% 日期: 28/04/2026
% 描述:
% 进行渗流场线性方程组求解，涵盖等效刚度矩阵 K_hat = Ct/(γ·dt) + Kt 组装、非零水头
% 边界条件施加、自由自由度线性求解、节点水头速率更新、约束自由度反力计算全过程。
% 第 0 次迭代（NRit=0）仅返回上一时间步速率，用于初始化刚度矩阵；第 1 次及以上迭代
% 执行完整线性求解与反力更新。
%
%--------------------------------------------------------------------------
% [dH, dQ, dHdt] = model_Solve_Flow(bc, Kt, Ct, oobf, NRit, fd, dt, H, dHdt_old)
%--------------------------------------------------------------------------
% 输入:
% bc        - 渗流边界条件矩阵 [nbc×2]，第1列为约束自由度编号，第2列为边界水头值
% Kt        - 渗流整体传导刚度矩阵 [nDoF×nDoF]（稀疏）
% Ct        - 渗流整体容量/质量矩阵 [nDoF×nDoF]（稀疏）
% oobf      - 渗流节点不平衡力（源项-内力）向量 [nDoF×1]
% NRit      - Newton-Raphson 迭代计数器（从 0 开始）
% fd        - 自由自由度索引向量 [nfree×1]
% dt        - 时间步长 [标量]
% H         - 当前时间步节点总水头向量 [nDoF×1]（用于计算速率）
% dHdt_old  - 上一时间步节点水头变化率 [nDoF×1]（NRit=0 时直接返回）
%--------------------------------------------------------------------------
% 输出:
% dH        - 节点水头增量向量 [nDoF×1]
% dQ        - 节点约束反力增量向量 [nDoF×1]（仅边界自由度非零）
% dHdt      - 更新后的节点水头变化率 [nDoF×1]
%--------------------------------------------------------------------------
% 关键算法:
%   1. 等效刚度: K_hat = Ct/(γ·dt) + Kt  （γ=1.0，向后差分隐式格式）
%   2. 边界处理: dH(bc) = 2·bc(:,2) （首步）或 0 （后续步，保持固定）
%   3. 自由求解: dH(fd) = K_hat(fd,fd) \ [oobf - K_hat(:,bc)·dH(bc)](fd)
%   4. 反力计算: dQ(bc) = K_hat(bc,:)·dH - oobf(bc)
%   5. 速率更新: dHdt = (H + dH) / (γ·dt)
%--------------------------------------------------------------------------
% 参见:
%   AMPLE_Solve_Flow_Main_BX - 渗流场迭代主控函数
%   detMPs_Flow_BX           - 渗流刚度/内力矩阵组装
%--------------------------------------------------------------------------

nDoF = length(oobf);                                                        % number of degrees of freedom
dH = zeros(nDoF,1);                                                         % zero reaction increment
dQ = zeros(nDoF,1);

if (NRit)>0
    ga = 1.00;
    dH(bc(:,1))=(1+sign(1-NRit))*bc(:,2);                                   % apply non-zero boundary conditions

    % --- 1. 等效刚度矩阵 ---
    K_hat = Ct/(ga*dt) + Kt;

    % --- 3. 边界条件处理 ---
    F_hat_mod = oobf - K_hat(:, bc(:,1)) * dH(bc(:,1));

    % 求解自由自由度上的位移
    dH(fd) = K_hat(fd, fd) \ F_hat_mod(fd);

    % --- 5. 更新加速度和速度 (必须使用求得的 u_new) ---
    dHdt = (H+dH)/(ga*dt);
    % --- 6. 计算约束自由度上的总反力 (可选) ---
    dQ(bc(:,1))=Kt(bc(:,1),:)*dH-oobf(bc(:,1));                          % determine reaction forces
else
    dHdt = dHdt_old;                                                        % v_{n+1}
end

end