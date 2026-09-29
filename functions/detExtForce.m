function [fext] = detExtForce(nodes,nD,g,mpData,g_wd,Contact,BS)

% fext外力计算子函数
%--------------------------------------------------------------------------
% 作者: FAN Haishan（改动范围10%） 
% 原始代码作者：William Coombs
% 【改动范围】 新增接触力，可指定重力维度，支持第一个周期重力线性增加
% 日期: 28/04/2026
% 描述:
% 通过节点力、重力、接触力，计算背景网格各节点的合节点力，用于后续求解器求解和不平衡力迭代
%
%--------------------------------------------------------------------------
% [fext] = detExtForce(nodes,nD,g,mpData,g_wd,Contact,BS)
%--------------------------------------------------------------------------
% 输入:
% nodes  - 节点总数
% nD     - 坐标维度
% g      - 自重应力
% mpData - 物质点信息结构体参数，包含: 
%           mpM   : 物质点质量
%           nIN   : 物质点关联节点序列
%           Svp   : 物质点形函数
%           fp    : 物质点上施加集中力
% g_wd   - 自重应力维度，指定'x','y','z'
% Contact- 接触力
% BS     -自重应力倍数，用于外部调整该分析步自重应力大小
%--------------------------------------------------------------------------
% 输出: 
% fext   - 节点外力序列 (nodes*nD,1)
%--------------------------------------------------------------------------

if nargin<5
   g_wd = 'z';
   Contact = zeros(size(mpData,2),nD);
   BS=zeros(nD,1);
end

if nargin<6
   Contact = zeros(size(mpData,2),nD);
   BS=zeros(nD,1);
end

if nargin<7
   BS=zeros(nD,1);
end

if g_wd == 'x'
    g_W = 1;
elseif g_wd == 'y'
    g_W = 2;
elseif g_wd == 'z'
    g_W = 3;
else
    error('g_wd only can be ''x'', ''y'', ''z''!')
end

nmp  = size(mpData,2);                                                      % number of material points & dimensions 
fext = zeros(nodes*nD,1);                                                   % zero the external force vector
grav = zeros(nD,1); grav(g_W) = -g;                                         % gavity vector
for mp = 1:nmp
   nIN = mpData(mp).nIN;                                                    % nodes associated with MP
   nn  = length(nIN);                                                       % number of nodes influencing the MP
   Svp = mpData(mp).Svp;                                                    % basis functions
   fp  = (mpData(mp).mpM*grav + BS .* mpData(mp).fp + Contact(mp,:)')*Svp;  % material point body & point nodal forces
   ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                     % nodel degrees of freedom 
   fext(ed) = fext(ed) + fp;                                                % combine into external force vector
end
end