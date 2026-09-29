function [eIN] = elemForMP(mesh,mpC,lp)

% 搜索物质点所属背景网格单元 / 子函数未修改，保留原始代码
%--------------------------------------------------------------------------
% 作者: William Coombs
% 日期:   06/05/2015
% 描述:
% 用于确定与材料点相关联的元素的函数，假设材料点的域关于粒子位置对称
%
%--------------------------------------------------------------------------
% [eIN] = elemForMP(mesh,mpC,lp)
%--------------------------------------------------------------------------
% 输入:
% mesh   - 结构体参数. 包含
%           - etpl  : 各单元节点信息 (nels,nen) 
%           - eMin  : 各网格最小坐标 (nels,nD)
%           - eMax  : 各网格最大坐标 (nels,nD)
% mpC   - 物质点坐标，按照维度拆分 (1,nD)
% lp    - 物质点尺寸
%--------------------------------------------------------------------------
% 输出:
% eIN   - 物质点关联背景网格序列
%--------------------------------------------------------------------------

Cmin = mesh.eMin;                                                           % element lower coordinate limit 
Cmax = mesh.eMax;                                                           % element upper coordinate limit 
etpl = mesh.etpl;                                                           % element topology
nD   = size(Cmin,2);                                                        % number of dimensions
nels = size(etpl,1);                                                        % number of elements
d    = sqrt(eps)*lp;                                                        % tolerance on the search (avoids zero overlap detection with GIMPM)
Pmin = mpC-lp+d;                                                            % particle domain extents (lower)
Pmax = mpC+lp-d;                                                            % particle domain extents (upper)
a    = true(nels,1);                                                        % initialise logical array
for i=1:nD 
    a = a.*((Cmin(:,i)<=Pmax(i)).*(Cmax(:,i)>=Pmin(i)));                      % element overlap with mp domain
end
eIN = find(a);                                                              % elements overlaps by the domain
end                                                                        