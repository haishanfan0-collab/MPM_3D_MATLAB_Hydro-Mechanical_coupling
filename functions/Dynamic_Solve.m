function [duvw,dudt,du2dt2,drct] = Dynamic_Solve(bc,gK,gM,gC,oobf,NRit,fd,uvw_old,dudt_old,du2dt2_old,dt)
% 【隐式动力线性方程组求解器】，内置 Newmark-β 等效刚度矩阵组装、边界条件处理、
% 位移增量求解、速度/加速度更新、约束反力计算全过程
%--------------------------------------------------------------------------
% 作者: FAN Haishan（基于 William Coombs 原始框架重构）
% 日期: 28/04/2026
% 描述:
% 进行隐式动力分析线性方程组求解，采用 Newmark-β 平均加速度法（γ=0.5, β=0.25，无条件稳定）。
% 涵盖等效刚度矩阵 K_hat = a0·M + a1·C + K 组装、非零位移边界条件施加、自由自由度线性求解、
% 节点加速度与速度更新、约束自由度反力计算全过程。第 0 次迭代（NRit=0）仅返回上一时间步
% 速度/加速度，用于初始化刚度矩阵；第 1 次及以上迭代执行完整隐式动力求解。
%
%--------------------------------------------------------------------------
% [duvw, dudt, du2dt2, drct] = 
%    Dynamic_Solve(bc, gK, gM, gC, oobf, NRit, fd, uvw_old, dudt_old, du2dt2_old, dt)
%--------------------------------------------------------------------------
% 输入:
% bc           - 位移边界条件矩阵 [nbc×2]，第1列为约束自由度编号，第2列为边界位移值
% gK           - 整体切线刚度矩阵 [nDoF×nDoF]（稀疏）
% gM           - 整体质量矩阵 [nDoF×nDoF]（稀疏）
% gC           - 整体阻尼矩阵 [nDoF×nDoF]（稀疏）
% oobf         - 节点不平衡力向量 [nDoF×1]
% NRit         - Newton-Raphson 迭代计数器（从 0 开始）
% fd           - 自由自由度索引向量 [nfree×1]
% uvw_old      - 上一时间步节点位移向量 [nDoF×1]
% dudt_old     - 上一时间步节点速度向量 [nDoF×1]
% du2dt2_old   - 上一时间步节点加速度向量 [nDoF×1]
% dt           - 时间步长 [标量]
%--------------------------------------------------------------------------
% 输出:
% duvw         - 节点位移增量向量 [nDoF×1]
% dudt         - 更新后的节点速度向量 [nDoF×1]（v_{n+1}）
% du2dt2       - 更新后的节点加速度向量 [nDoF×1]（a_{n+1}）
% drct         - 节点约束反力增量向量 [nDoF×1]（仅边界自由度非零）
%--------------------------------------------------------------------------
% 关键算法（Newmark-β 平均加速度法，γ=0.5, β=0.25）:
%   1. 积分常数: a0=1/(β·dt²), a1=γ/(β·dt), a2=1/(β·dt), a3=1/(2β)-1,
%               a6=(1-γ)·dt, a7=γ·dt
%   2. 等效刚度: K_hat = a0·M + a1·C + K
%   3. 边界处理: duvw(bc) = 2·bc(:,2)（首步）或 0（后续步，保持固定）
%   4. 自由求解: duvw(fd) = K_hat(fd,fd) \ [oobf - K_hat(:,bc)·duvw(bc)](fd)
%   5. 加速度更新: a_{n+1} = a0·(u_n + Δu) - a2·v_n - a3·a_n
%   6. 速度更新:   v_{n+1} = v_n + a6·a_n + a7·a_{n+1}
%   7. 反力计算:   drct(bc) = K_hat(bc,:)·Δu - oobf(bc)
%--------------------------------------------------------------------------
% 参见:
%   AMPLE_Solve_Main_BX - 隐式动力迭代主控函数
%   detMPs_BX           - 动力刚度/质量/阻尼矩阵组装
%--------------------------------------------------------------------------
nDoF = length(oobf);                                                        % number of degrees of freedom 
drct = zeros(nDoF,1);                                                       % zero reaction increment
duvw = zeros(nDoF,1);
if (NRit)>0
    duvw(bc(:,1))=(1+sign(1-NRit))*bc(:,2);                                 % apply non-zero boundary conditions

    %计算整体位移向量
    % 假定迭代参数
    ga = 0.50; bt=0.25;
    % 模型参数计算
    a0 = 1/(bt*dt*dt);a1 = ga/(bt*dt);a2 = 1/(bt*dt);
    a3 = 1/(2*bt)-1; 
    a6 = (1-ga)*dt;  a7 = ga*dt;

    % --- 1. 等效刚度矩阵 ---
    K_hat = a0*gM + a1*gC + gK;

    % --- 3. 边界条件处理 ---
    duvw(fd)=K_hat(fd,fd)\(oobf(fd)-gK(fd,bc(:,1))*duvw(bc(:,1)));

    % --- 5. 更新加速度和速度 (必须使用求得的 u_new) ---
    du2dt2 = a0 * (uvw_old+duvw) - a2*dudt_old - a3*du2dt2_old;                % a_{n+1}
    dudt   = dudt_old + a6*du2dt2_old + a7*du2dt2;                             % v_{n+1}

    % --- 6. 计算约束自由度上的总反力 (可选) ---
    drct(bc(:,1))=gK(bc(:,1),:)*duvw-oobf(bc(:,1));                         % determine reaction forces 
else
    du2dt2 = du2dt2_old;                % a_{n+1}
    dudt   = dudt_old;                    % v_{n+1}
end

end

