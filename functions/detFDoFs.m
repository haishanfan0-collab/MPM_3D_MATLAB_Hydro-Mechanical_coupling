function [fd] = detFDoFs(mesh)
% 确定力学计算激活背景网格序列，识别不包含物质点的背景网格 / 子函数未修改，保留原始代码
%--------------------------------------------------------------------------
% 作者: William Coombs
% 日期:   17/12/2018
% 描述:
% 根据网格信息，确定需要求解的背景网格序列，用于求解器识别
%
%--------------------------------------------------------------------------
% [fd] = detFDoFs(mesh)
%--------------------------------------------------------------------------
% 输入:
% mesh  - 网格结构体参数. 包含:
%           - etpl  : 各单元节点信息 (nels,nen) 
%           - eInA  : 激活的背景网格节点信息 
%           - bc   : 边界条件，通常指位移约束 (*,2)
%--------------------------------------------------------------------------
% 输出:
% fd    - 激活的背景网格节点序列 (*,1)
%--------------------------------------------------------------------------

[nodes,nD] = size(mesh.coord);                                              % no. nodes and dimensions
nDoF   = nodes*nD;                                                          % no. degrees of freedom
incN   = unique(mesh.etpl(mesh.eInA>0,:));                                  % unique active node list
iN     = size(incN,1);                                                      % number of nodes in the list
incDoF = reshape(ones(nD,1)*incN'*nD-(nD-1:-1:0).'*ones(1,iN),1,iN*nD);     % active degrees of freedom
fd     = (1:nDoF);                                                          % all degrees of freedom
fd(mesh.bc(:,1)) = 0;                                                       % zero fixed displacement BCs
fd     = fd(incDoF);                                                        % only include active DoF 
fd     = fd(fd>0);                                                          % remove fixed displacement BCs
end