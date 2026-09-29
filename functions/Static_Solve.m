function [duvw,drct] = Static_Solve(bc,Kt,oobf,NRit,fd)
% 【隐式静力线性方程组求解器】，内置边界条件处理、直接法位移增量求解、约束反力计算
%--------------------------------------------------------------------------
% 作者: FAN Haishan（基于 William Coombs 原始框架重构）
% 日期: 28/04/2026
% 描述:
% 进行隐式静力分析线性方程组求解，采用直接法（稀疏矩阵左除）求解位移增量。涵盖非零
% 位移边界条件施加、自由自由度线性方程组求解、约束自由度反力计算全过程。第 0 次迭代
% （NRit=0）返回零向量，用于初始化刚度矩阵；第 1 次及以上迭代执行完整静力求解。
%
%--------------------------------------------------------------------------
% [duvw, drct] = Static_Solve(bc, Kt, oobf, NRit, fd)
%--------------------------------------------------------------------------
% 输入:
% bc     - 位移边界条件矩阵 [nbc×2]，第1列为约束自由度编号，第2列为边界位移值
% Kt     - 整体切线刚度矩阵 [nDoF×nDoF]（稀疏）
% oobf   - 节点不平衡力向量 [nDoF×1]
% NRit   - Newton-Raphson 迭代计数器（从 0 开始）
% fd     - 自由自由度索引向量 [nfree×1]
%--------------------------------------------------------------------------
% 输出:
% duvw   - 节点位移增量向量 [nDoF×1]
% drct   - 节点约束反力向量 [nDoF×1]（仅边界自由度非零）
%--------------------------------------------------------------------------
% 关键算法:
%   1. 边界处理: duvw(bc) = 2·bc(:,2)（首步）或 0（后续步，保持固定）
%   2. 自由求解: duvw(fd) = Kt(fd,fd) \ [oobf(fd) - Kt(fd,bc)·duvw(bc)]
%   3. 反力计算: drct(bc) = Kt(bc,:)·duvw - oobf(bc)
%--------------------------------------------------------------------------
% 参见:
%   AMPLE_Solve_Main_BX - 隐式静力迭代主控函数
%   detMPs_BX           - 静力刚度矩阵与内力组装
%--------------------------------------------------------------------------
nDoF = length(oobf);                                                        % number of degrees of freedom 
duvw = zeros(nDoF,1);                                                       % zero displacement increment
drct = zeros(nDoF,1);                                                       % zero reaction increment
if (NRit)>0                                                             
    duvw(bc(:,1))=(1+sign(1-NRit))*bc(:,2);                                 % apply non-zero boundary conditions
    duvw(fd)=Kt(fd,fd)\(oobf(fd)-Kt(fd,bc(:,1))*duvw(bc(:,1)));             % solve for displacements
    drct(bc(:,1))=Kt(bc(:,1),:)*duvw-oobf(bc(:,1));                         % determine reaction forces 
end
end