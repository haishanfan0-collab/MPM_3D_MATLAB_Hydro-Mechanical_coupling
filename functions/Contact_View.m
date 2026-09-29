function Contact_message = Contact_View(meshA,meshB,mpDataB)

% 边界物质点碰撞检测代码
%--------------------------------------------------------------------------
% Author: FAN Haishan
% Date:   04/03/2026
% 描述:
% 通过A和B的的网格和物质点信息，识别出来碰撞，并形成接触对
%
%--------------------------------------------------------------------------
% Contact_message = Contact_View(meshA,mpDATAA,meshB,mpDATAB)
%--------------------------------------------------------------------------
% Input(s):
% meshA,mpDATAA   物质点团A的信息
% meshB,mpDATAB   物质点团B的信息

%--------------------------------------------------------------------------
% Ouput(s);
% Contact_message - 接触对信息（结构体参数）

%--------------------------------------------------------------------------
if nargin<1
    [~,mpDataA,meshA] = setupGrid_beam;
    [~,mpDataB,meshB] = setupGrid_beam_down;
    mpUA = reshape([mpDataA.mpC]',2,size(mpDataA,2));
    mpUB = reshape([mpDataB.mpC]',2,size(mpDataB,2));
    [meshA,~,~,meshA.bounday_imfo] = elemMPinfo(meshA,mpDataA,'static'); 
    [meshB,mpDataB,~,meshB.bounday_imfo] = elemMPinfo(meshB,mpDataB,'static'); 
    figure(1)
    hold on
    plot(mpUA(1,1441:end),mpUA(2,1441:end),'ro')
    plot(mpUB(1,1441:end),mpUB(2,1441:end),'b.')
end

%% 相邻边界单元检测
Contact_message = [];
mesh_inter = intersect(meshA.bounday_imfo.mesh_boundry, meshB.bounday_imfo.mesh_boundry);
mesh_boundaryA = meshA.bounday_imfo.node_boundry_mesh;
mesh_boundaryB = meshB.bounday_imfo.node_boundry_mesh;

%% MPM-A 处于相邻单元的边界物质点序列搜索
node_ContactA = zeros(length(mesh_boundaryA),3);
Find_ID = 1;
for i1 = 1:length(mesh_boundaryA)
    Find_node = intersect(mesh_boundaryA{i1}, mesh_inter);
    if ~isempty(Find_node)
        node_ContactA(Find_ID,:) = meshA.bounday_imfo.node_boundry(i1,:);
        Find_ID = Find_ID + 1;
    end
end
node_ContactA(node_ContactA(:,1)==0,:) = [];

%% MPM-B 处于相邻单元的边界物质点序列搜索
node_ContactB = zeros(length(mesh_boundaryB),3);
Find_ID = 1;
for i2 = 1:length(mesh_boundaryB)
    Find_node = intersect(mesh_boundaryB{i2}, mesh_inter);
    if ~isempty(Find_node)
        node_ContactB(Find_ID,:) = meshB.bounday_imfo.node_boundry(i2,:);
        Find_ID = Find_ID + 1;
    end
end
node_ContactB(node_ContactB(:,1)==0,:) = [];

%% 组建接触信息结构体数据
mpCB = [(1:size(mpDataB,2))',reshape([mpDataB.mpC],size(mpDataB(1).mpC,2),size(mpDataB,2))'];
mpCB(end-size(mesh_boundaryB,1)+1:end,:) = [];

%% 组装接触信息
for i3 = size(node_ContactA,1):-1:1
    Contact_message(i3).ID = i3;                                 % 接触对ID，A表示从表面，B表示主表面
    Contact_message(i3).Cong_ID = node_ContactA(i3,1);           % 接触对从表面物质点 ID
    Contact_message(i3).Cong_Coord = node_ContactA(i3,2:end);    % 接触对从表面物质点 坐标

    % 搜索距离从节点距离最近的主表面物质点(前3)
    DistanceAB = sortrows([(1:length(node_ContactB))',vecnorm(node_ContactB(:,2:end)-node_ContactA(i3,2:end),2,2)],2);
    Select_B = node_ContactB(DistanceAB(1:3,1),:);
    [~,LS] = min(vecnorm(Select_B(:,2:end)-mean(Select_B(:,2:end)),2,2)); 
    if LS == 1
        Select_B = Select_B([2,1,3],:);
    elseif LS == 3
        Select_B = Select_B([1,3,2],:);
    end

    Contact_message(i3).Main_ID = Select_B(:,1);                 % 接触对从表面物质点 ID
    Contact_message(i3).Main_Coord = Select_B(:,2:end);          % 接触对从表面物质点 坐标

    % 搜索主表面物质点最近的核心物质点序列(y距离最远，x最近） 
    l1 = distToPerpendicularBisector(mpCB(:,2:end),Select_B(1,2:end),Select_B(2,2:end));
    L2 = distToLine(mpCB(:,2:end),Select_B(1,2:end),Select_B(2,2:end));
    LS = sortrows([mpCB(:,1),l1+1./(0.05+L2)],2);
    Contact_message(i3).Main_MPID = mpCB(LS(1),1);
    Contact_message(i3).Main_MPCoord = mpCB(LS(1),2:end);

    % 获得主表面的法线和切线信息
    Contact_message(i3).Main_Tan1 = triangleOutwardNormals([mpCB(LS(1),2:end);Select_B(1:2,2:end)]);
    Contact_message(i3).Main_Tan2 = triangleOutwardNormals([mpCB(LS(1),2:end);Select_B(2:3,2:end)]);

    Contact_message(i3).Main_Nor1 = (Select_B(1,2:end) - Select_B(2,2:end))/norm((Select_B(1,2:end) - Select_B(2,2:end)));
    Contact_message(i3).Main_Nor2 = (Select_B(2,2:end) - Select_B(3,2:end))/norm((Select_B(2,2:end) - Select_B(3,2:end)));


    % 求解最近的两个三角形的距离
    Contact_message(i3).dist1 = pointToTriangleDistance([mpCB(LS(1),2:end);Select_B(1:2,2:end)],Contact_message(i3).Cong_Coord);
    Contact_message(i3).dist2 = pointToTriangleDistance([mpCB(LS(1),2:end);Select_B(2:3,2:end)],Contact_message(i3).Cong_Coord);
end

end

function d = distToLine(C, A, B)
% C: [N,2], A/B: [1,2], 返回C到直线AB的距离 [N,1]
   
   AB = (B - A)';
   AC = C - A;
   
   % 叉积绝对值 / |AB|
   d = abs(AC(:,1).*AB(2) - AC(:,2).*AB(1)) / norm(AB);
end

function d = distToPerpendicularBisector(C, A, B)
% C: [N,2], A/B: [1,2]或[2,1]，返回d: [N,1]
   
   M = (A + B) / 2;           % 中点
   AB = (B - A)';             % 确保列向量 [2,1]
   MC = C - M;                % [N,2]
   
   d = abs(MC * AB) / norm(AB);  % [N,1]
end


%% 辅助函数
function dist = pointToTriangleDistance(A, B)
% pointToTriangleDistance 计算点B到三角形A的最短距离
% 输入: A - 3x2矩阵，三角形的三个顶点坐标 [x1,y1; x2,y2; x3,y3]
%      B - 1x2向量，待测点坐标 [x,y]
% 输出: dist - 最短距离（三角形内部为负值，外部为正值）

    % 提取顶点
    P1 = A(1,:);
    P2 = A(2,:);
    P3 = A(3,:);
    
    % 1. 首先判断点B是否在三角形内部
    if isPointInTriangle(B, P1, P2, P3)
        % 内部：计算到三条边的距离，取最小值，返回负值
        dist = -pointToLineSegment(B, P2, P3);
    else
        % 外部：计算到三条边的距离（点到线段），取最小值，返回正值
        dist = pointToLineSegment(B, P2, P3);
    end
end


%% 辅助函数：判断点是否在三角形内部（使用重心坐标法）
function inside = isPointInTriangle(P, A, B, C)
    % 使用叉积法判断点P是否在三角形ABC内部
    % 计算三个叉积，如果同号则在内部
    
    % 向量AB, BC, CA
    AB = B - A;
    BC = C - B;
    CA = A - C;
    
    % 向量AP, BP, CP
    AP = P - A;
    BP = P - B;
    CP = P - C;
    
    % 计算叉积（二维叉积的z分量）
    cross1 = AB(1)*AP(2) - AB(2)*AP(1);  % AB x AP
    cross2 = BC(1)*BP(2) - BC(2)*BP(1);  % BC x BP
    cross3 = CA(1)*CP(2) - CA(2)*CP(1);  % CA x CP
    
    % 判断符号：如果三个叉积同号（或至少一个为零），则在内部或边上
    inside = (cross1 >= 0 && cross2 >= 0 && cross3 >= 0) || ...
             (cross1 <= 0 && cross2 <= 0 && cross3 <= 0);
end

%% 辅助函数：计算点到线段的距离
function dist = pointToLineSegment(P, A, B)
    % 计算点P到线段AB的最短距离
    
    AB = B - A;
    AP = P - A;
    
    % 计算投影参数t
    len2 = dot(AB, AB);  % |AB|^2
    
    if len2 == 0
        % A和B重合，退化为点到点距离
        dist = norm(P - A);
        return;
    end
    
    % 投影点在线段上的参数 t = [(P-A)·(B-A)] / |B-A|^2
    t = max(0, min(1, dot(AP, AB) / len2));
    
    % 投影点坐标
    projection = A + t * AB;
    
    % 距离
    dist = norm(P - projection);
end

function N = triangleOutwardNormals(A)
    % 计算边向量
    E = [A(2,:)-A(1,:); A(3,:)-A(2,:); A(1,:)-A(3,:)];
    
    % 左侧法向（逆时针旋转90度）
    N_left = [-E(:,2), E(:,1)];
    
    % 通过有符号面积判断顶点顺序
    signed_area = 0.5 * ((A(2,1)-A(1,1))*(A(3,2)-A(1,2)) - ...
                         (A(3,1)-A(1,1))*(A(2,2)-A(1,2)));
    
    % 逆时针时左侧是内侧，需取反；顺时针时左侧就是外侧
    N = N_left * sign(-signed_area);  % signed_area>0(逆时针)时取反
    
    % 单位化
    N = N ./ sqrt(sum(N.^2, 2));
    N = N(2,:);
end