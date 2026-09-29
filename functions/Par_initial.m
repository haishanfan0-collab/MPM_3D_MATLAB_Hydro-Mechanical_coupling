function [uvw,VA,Stifiness,oobf,fErr,frct,NRit] = Par_initial(nDoF)
% 初始化代码
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 初始化参数，提前分配内存
%
%--------------------------------------------------------------------------
% [uvw,VA,Stifiness,oobf,fErr,frct,NRit] = Par_initial(nDoF)
%--------------------------------------------------------------------------
% 输入:
% nDoF  - 维度，大小为 节点维度*节点个数
%--------------------------------------------------------------------------
% 输出:
% uvw        - 位移初始化
% VA         - 速度加速度初始化
% Stifiness  - 刚度矩阵初始化
% oobf       - 不平衡力初始化
% fErr       - 迭代误差初始化
% frct       - 节点内力初始化
% NRit       - 迭代次数初始化
%--------------------------------------------------------------------------

uvw  = zeros(nDoF,1);                                                       % zeros displacements (for plotting function)
VA.dudt  = zeros(nDoF,1);
VA.du2dt2  = zeros(nDoF,1);
Stifiness.Kt   = sparse(nDoF,nDoF);                                                   % zero global stiffness matrix
Stifiness.Mt   = sparse(nDoF,nDoF);                                                   % zero global mass matrix
Stifiness.Ct   = sparse(nDoF,nDoF);                                                   % zero global damp matrix
oobf = zeros(nDoF,1);                                                       % initial out-of-balance force
fErr = 1;                                                                   % initial error
frct = zeros(nDoF,1);                                                       % zero the reaction forces
NRit = 0;

end