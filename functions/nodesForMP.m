function [nodes] = nodesForMP(etpl,elems)

% 获取背景网格节点 / 子函数未修改，保留原始代码
%--------------------------------------------------------------------------
% 作者: William Coombs
% 日期:   06/05/2015
% 描述:
% 给定背景网格序列，返回背景网格包含节点信息
%
%--------------------------------------------------------------------------
% [nodes] = NODESFORMP(etpl,elems)
%--------------------------------------------------------------------------
% 输入:
% etpl  - 背景网格节点信息
% elems - 背景网格序列
%--------------------------------------------------------------------------
% 输出: 
% nodes - 节点序列
%--------------------------------------------------------------------------

nen=size(etpl,2);                                                           % number of nodes per element
n=size(elems,1)*size(elems,2);                                              % number of elements in group
nn=n*nen;                                                                   % number of nodes (inc. duplicates)
e=sort(reshape(etpl(elems,:),nn,1),1);                                      % list of all nodes (inc. duplicates)
d=[1; e(2:nn)-e(1:nn-1)]>0;                                                 % duplicate removal
nodes=e(d);                                                                 % unique list of nodes
end