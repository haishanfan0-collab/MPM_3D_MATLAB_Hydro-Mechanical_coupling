function [H,HVA,Flow,oobf_Flow,fErr_Flow,frct_Flow,NRit_Flow] = Par_initial_Flow(nDoF)
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
% H          - 水头初始化
% HVA        - dHdt初始化
% Flow       - 流体整体传导矩阵、容水度矩阵初始化
% oobf_Flow  - 不平衡流量初始化
% fErr_Flow  - 渗流误差参数初始化
% frct_Flow  - 节点内流量初始化
% NRit_Flow  - 渗流误差参数初始化
%--------------------------------------------------------------------------

H  = zeros(nDoF,1);                                                         % 水头矩阵
HVA.dHdt  = zeros(nDoF,1);                                                  % 水头变化率

Flow.Kt   = sparse(nDoF,nDoF);                                              % 整体传导矩阵
Flow.Ct   = sparse(nDoF,nDoF);                                              % 整体容储矩阵
Flow.Mt   = sparse(nDoF,nDoF);                                              % 整体容储矩阵

oobf_Flow = zeros(nDoF,1);                                                  % 不平衡流量
frct_Flow = zeros(nDoF,1);                                                  % 边界流量
fErr_Flow = 1;                                                              % 误差

NRit_Flow = 0;  

end