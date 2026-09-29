function [fd] = detFDoFs_Flow(mesh)
% 确定渗流计算激活背景网格序列，识别不包含物质点的背景网格 / 新增子函数
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期:   28/04/2026
% 描述:
% 根据网格信息，确定需要求解的背景网格序列，用于求解器识别
%
%--------------------------------------------------------------------------
% [fd] = detFDoFs_Flow(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh  - 网格结构体参数. 包含:
%           - etpl   : 各单元节点信息 (nels,nen) 
%           - eInA   : 激活的背景网格节点信息 
%           - bc_Flow: 渗流边界条件，通常指位移约束 (*,2)
%--------------------------------------------------------------------------
% 输出:
% fd    - 激活的背景网格节点序列 (*,1)
%--------------------------------------------------------------------------

nDoF = size(mesh.coord,1);                                                  % no. nodes and dimensions

incN   = unique(mesh.etpl(mesh.eInA>0,:));                                  % unique active node list

fd     = (1:nDoF);                                                          % all degrees of freedom
fd(mesh.bc_Flow(:,1)) = 0;                                                  % zero fixed displacement BCs
fd     = fd(incN);                                                          % only include active DoF 
fd     = fd(fd>0);                                                          % remove fixed displacement BCs
end