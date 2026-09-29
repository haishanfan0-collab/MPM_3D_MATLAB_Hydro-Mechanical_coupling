function [Contact_message,Force_Cong,Force_Main,Force_ADD,Contact_SL] = Contact_Update(Contact_par,Contact_message,mpDataA,uvwA,mpDataB,uvwB,dt,Contact_Force_Up,Contact_Force_Down,separate_YES)

% 边界物质点应力迭代代码
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
% Contact_par     接触面参数，法向刚度/切向阻尼
% Contact_message 目前接触面信息
% mpDataA uvwA    物质点团A的信息，位移
% mpDataB uvwB    物质点团B的信息，位移
% dt              当前时间步，用于计算切向速度求解切向力
% Contact_Force   目前总的接触力
%--------------------------------------------------------------------------
% Ouput(s);
% Contact_message - 接触对信息（结构体参数）

%--------------------------------------------------------------------------
Contact_SL = 0.0;

Contact_TAN_E = Contact_par(1);
Contact_NOR_U = Contact_par(2);
Contact_NOR_V = Contact_par(3);

if nargin<10
    separate_YES = 'separate_YES';
end

if strcmp(separate_YES, 'separate_YES')
    ALLOW = -1e-6;
else
    ALLOW1 = [mpDataA.mpType];
    ALLOW = [mpDataA.lp];
    ALLOW = ALLOW(2:2:end);
    ALLOW = min(ALLOW(ALLOW1))*0.50;
end


%% 获取从节点ID序列-A
nD  = length(mpDataA(1).mpC);
if ~isempty(Contact_message)
node_ContactA = [[Contact_message.Cong_ID]',(reshape([Contact_message.Cong_Coord]',2,size(Contact_message,2)))'];
% 更新位移
if sum(abs(uvwA))>0
    nD  = length(mpDataA(1).mpC);                                               % number of dimensions
    for i=1:size(node_ContactA,1)
        mp = node_ContactA(i,1);
        nIN = mpDataA(mp).nIN;                                                  % nodes associated with material point
        nn  = length(nIN);                                                      % number nodes
        N   = mpDataA(mp).Svp;                                                  % basis functions
        ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
        node_ContactA(i,2:end)   = node_ContactA(i,2:end) + N*uvwA(ed);          % update material point coordinates
    end
end

%% 获取从节点ID序列-B
node_ContactB = zeros(3*size(Contact_message,2),3);
node_ContactB(:,1) = reshape([Contact_message.Main_ID],3*size(Contact_message,2),1);
LS = [Contact_message.Main_Coord]';
for i = 1:size(Contact_message,2)
    for j =1:3
        node_ContactB((i-1)*3+j,2:end) = [LS((i-1)*2+1:i*2,j)]';
    end
end

% 更新位移
if sum(abs(uvwB))>0
                                                   % number of dimensions
    for i=1:size(node_ContactB,1)
        mp = node_ContactB(i,1);
        nIN = mpDataB(mp).nIN;                                                  % nodes associated with material point
        nn  = length(nIN);                                                      % number nodes
        N   = mpDataB(mp).Svp;                                                  % basis functions
        ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
        node_ContactB(i,2:end)   = node_ContactB(i,2:end) + N*uvwB(ed);          % update material point coordinates
    end
end

%% 获取主要物质点ID序列-B
mpCB = [[Contact_message.Main_MPID]',(reshape([Contact_message.Main_MPCoord]',2,size(Contact_message,2)))'];
% 更新位移
if sum(abs(uvwB))>0
    nD  = length(mpDataA(1).mpC);                                               % number of dimensions
    for i=1:size(mpCB,1)
        mp = mpCB(i,1);
        nIN = mpDataB(mp).nIN;                                                  % nodes associated with material point
        nn  = length(nIN);                                                      % number nodes
        N   = mpDataB(mp).Svp;                                                  % basis functions
        ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
        mpCB(i,2:end)   = mpCB(i,2:end) + N*uvwB(ed);                            % update material point coordinates
    end
end

%% 组装接触信息
Force_Cong = zeros(size(mpDataA,2),nD);
Force_Main = zeros(size(mpDataB,2),nD);
for i3 = 1:size(Contact_message,2)
    Contact_message(i3).Cong_CoordNEW = node_ContactA(i3,2:end);    % 接触对从表面物质点 坐标
    
    Select_B = node_ContactB((i3-1)*3+1:i3*3,:);
    % 获得主表面的法线和切线信息
    Contact_message(i3).Main_CoordNEW = Select_B(:,2:end);          % 接触对从表面物质点 坐标

    % 搜索主表面物质点最近的核心物质点序列
    Contact_message(i3).Main_MPCoordNEW = mpCB(i3,2:end);

    % 获得主表面的法线和切线信息
    Contact_message(i3).Main_Tan1NEW = triangleOutwardNormals([mpCB(i3,2:end);Select_B(1:2,2:end)]);
    Contact_message(i3).Main_Tan2NEW = triangleOutwardNormals([mpCB(i3,2:end);Select_B(2:3,2:end)]);

    Contact_message(i3).Main_Nor1NEW = (Select_B(1,2:end) - Select_B(2,2:end))/norm((Select_B(1,2:end) - Select_B(2,2:end)));
    Contact_message(i3).Main_Nor2NEW = (Select_B(2,2:end) - Select_B(3,2:end))/norm((Select_B(2,2:end) - Select_B(3,2:end)));


    % 求解最近的两个三角形的距离
    Contact_message(i3).dist1NEW = pointToTriangleDistance([mpCB(i3,2:end);Select_B(1:2,2:end)],Contact_message(i3).Cong_CoordNEW);
    Contact_message(i3).dist2NEW = pointToTriangleDistance([mpCB(i3,2:end);Select_B(2:3,2:end)],Contact_message(i3).Cong_CoordNEW);
    
    if ~isfield(Contact_message(i3), 'Contact_YES') || isempty(Contact_message(i3).Contact_YES)
        Contact_message(i3).Contact_YES = 0.0;
        Contact_message(i3).Force_NOR_OLD = 0.0;
        Contact_message(i3).dist1_OLD = 0.0;
        Contact_message(i3).dist2_OLD = 0.0;
        Contact_message(i3).AHA = 1.0;
        Contact_message(i3).Force_NOR_OLD = zeros(1,nD);
        Contact_message(i3).Force_TAN_OLD = zeros(1,nD);
        Contact_message(i3).JS = 0;
    end

    % 求解接触力
    if Contact_message(i3).dist1NEW<-1e-6 || Contact_message(i3).Contact_YES == 1
        Contact_message(i3).Contact_YES = 1;

        if Contact_message(i3).dist1_OLD*Contact_message(i3).dist1NEW<0
            Contact_message(i3).AHA = max(0.01,Contact_message(i3).AHA * 0.50);
        end

        Contact_SL = Contact_SL + 1e-6*Contact_TAN_E*Contact_message(i3).AHA;

        % 求解法向侵入（方向以作用于主面为正）
        B1 = (Contact_message(i3).Main_Tan1NEW+Contact_message(i3).Main_Tan1)/2;
        B2 = (Contact_message(i3).Main_Nor1NEW+Contact_message(i3).Main_Nor1)/2;
        Force_TAN = Contact_TAN_E*Contact_message(i3).AHA*B1*...
            Contact_message(i3).dist1NEW;

        % 目前总的法向力计算(总力+该步法向力,在垂直方向上的投影)
        A = (-Contact_Force_Up(Contact_message(i3).Cong_ID,:)+Force_TAN);
        
        Force_TANLS = abs(dot(A, B1) / norm(B1));
        TAN_ALL = (dot(A,B1) / dot(B1,B1)) * B1;

        % 计算初始时刻对于节点1和2的距离
        dis_old = pointProjection(Contact_message(i3).Cong_Coord,Contact_message(i3).Main_Coord(1:2,:));
        dis_new = pointProjection(Contact_message(i3).Cong_CoordNEW,Contact_message(i3).Main_CoordNEW(1:2,:));
        L_mean = (Contact_message(i3).Main_Coord(1:2,:) + Contact_message(i3).Main_CoordNEW(1:2,:))/2;
        L_mean = norm(L_mean(1,:)-L_mean(2,:));
        L_Dis = L_mean*(dis_new-dis_old);

        % Force_NOR = Force_TANLS*Contact_NOR_U/Contact_NOR_V*...
        %             (Contact_message(i3).Main_Nor1NEW+Contact_message(i3).Main_Nor1)/2*(L_Dis/dt);
        NOR_ALL = Force_TANLS*Contact_NOR_U*...
            B2*tanh(L_Dis/dt/Contact_NOR_V);

        % 主表面分配比例
        Main_BL = [max(1-(dis_old+dis_new)/2,0),min((dis_old+dis_new)/2),1];
        Main_ID = Contact_message(i3).Main_ID(1:2);

    elseif Contact_message(i3).dist2NEW<-1e-6 || Contact_message(i3).Contact_YES == 2
        Contact_message(i3).Contact_YES = 2;

        if Contact_message(i3).dist2_OLD*Contact_message(i3).dist2NEW<0
            Contact_message(i3).AHA = max(0.01,Contact_message(i3).AHA * 0.50);
        end

        Contact_SL = Contact_SL + 1e-6*Contact_TAN_E*Contact_message(i3).AHA;

        % 求解法向侵入（方向以作用于主面为正）
        B1 = (Contact_message(i3).Main_Tan2NEW+Contact_message(i3).Main_Tan2)/2;
        B2 = (Contact_message(i3).Main_Nor2NEW+Contact_message(i3).Main_Nor2)/2;
        Force_TAN = Contact_TAN_E*Contact_message(i3).AHA*B1*...
                Contact_message(i3).dist2NEW;

        % 目前总的法向力计算(总力+该步法向力,在垂直方向上的投影)
        A = (-Contact_Force_Up(Contact_message(i3).Cong_ID,:)+Force_TAN);

        Force_TANLS = abs(dot(A, B1) / norm(B1));
        TAN_ALL = (dot(A,B1) / dot(B1,B1)) * B1;

        % 计算初始时刻对于节点1和2的距离
        dis_old = pointProjection(Contact_message(i3).Cong_Coord,Contact_message(i3).Main_Coord(2:3,:));
        dis_new = pointProjection(Contact_message(i3).Cong_CoordNEW,Contact_message(i3).Main_CoordNEW(2:3,:));
        L_mean = (Contact_message(i3).Main_Coord(2:3,:) + Contact_message(i3).Main_CoordNEW(2:3,:))/2;
        L_mean = norm(L_mean(1,:)-L_mean(2,:));
        L_Dis = L_mean*(dis_new-dis_old);

        % Force_NOR = Force_TANLS*Contact_NOR_U/Contact_NOR_V*...
        %             (Contact_message(i3).Main_Nor2NEW+Contact_message(i3).Main_Nor2)/2*(L_Dis/dt);
        NOR_ALL = Force_TANLS*Contact_NOR_U*...
            B2*tanh(L_Dis/dt/Contact_NOR_V);
        
        % 主表面分配比例
        Main_BL = [max(1-(dis_old+dis_new)/2,0),min((dis_old+dis_new)/2,1)];
        Main_ID = Contact_message(i3).Main_ID(2:3);
    else
        Contact_message(i3).Contact_YES = 0;
        Force_TAN = zeros(1,nD);
        TAN_ALL = zeros(1,nD);
        NOR_ALL = zeros(1,nD);
        Main_BL = zeros(1,2);
        Main_ID = Contact_message(i3).Main_ID(1:2);
    end

    % 获取物质点附加应力
    if Contact_message(i3).Force_NOR_OLD~=0.0
        NOR_ALLCal = (NOR_ALL*0.99+Contact_message(i3).Force_NOR_OLD*0.01);
    else
        NOR_ALLCal = NOR_ALL;
    end
    Contact_message(i3).Force_NOR_OLD = NOR_ALLCal;
    Contact_message(i3).Force_TAN_OLD = Force_TAN;
    

    Contact_message(i3).dist1_OLD = Contact_message(i3).dist1NEW;
    Contact_message(i3).dist2_OLD = Contact_message(i3).dist2NEW;

    Contact_message(i3).Cong_Force = [Contact_message(i3).Cong_ID,-NOR_ALLCal-TAN_ALL];
    Contact_message(i3).Main_Force = [Main_ID,[(NOR_ALL+TAN_ALL)*Main_BL(1);(NOR_ALLCal+TAN_ALL)*Main_BL(2)]];
    Force_Cong(Contact_message(i3).Cong_ID,:) = Force_Cong(Contact_message(i3).Cong_ID,:) + Contact_message(i3).Cong_Force(:,2:end);
    Force_Main(Main_ID,:) = Force_Main(Main_ID,:) + Contact_message(i3).Main_Force(:,2:end);
    
    Contact_message(i3).JS = Contact_message(i3).JS + 1;
    % 更新接触坐标，用于计算下一次的不平衡力
    % Contact_message(i3).Cong_Coord = Contact_message(i3).Cong_CoordNEW;
    % Contact_message(i3).Main_Coord = Contact_message(i3).Main_CoordNEW;
    % Contact_message(i3).Main_MPCoord = Contact_message(i3).Main_MPCoordNEW;
    % Contact_message(i3).Main_Tan1 = Contact_message(i3).Main_Tan1;
    % Contact_message(i3).Main_Tan2 = Contact_message(i3).Main_Tan2;
    % Contact_message(i3).Main_Nor1 = Contact_message(i3).Main_Nor1;
    % Contact_message(i3).Main_Nor2 = Contact_message(i3).Main_Nor2;
    % 
    % Contact_message(i3).dist1 = Contact_message(i3).dist1NEW;
    % Contact_message(i3).dist2 = Contact_message(i3).dist2NEW;
end
Contact_Force_Up(Force_Cong(:,1)==0) = 0;
Contact_Force_Down(Force_Main(:,1)==0) = 0;

Force_Cong = Force_Cong - Contact_Force_Up;
Force_Main = Force_Main - Contact_Force_Down;


Force_TAN = sum([reshape([Contact_message.Force_TAN_OLD]',2,size(Contact_message,2))],2);
NOR_ALL = sum([reshape([Contact_message.Force_NOR_OLD]',2,size(Contact_message,2))],2);
Force_ADD = [Force_TAN,NOR_ALL]';
else
    Force_Cong = zeros(size(mpDataA,2),nD);
    Force_Main = zeros(size(mpDataB,2),nD);
    Force_ADD = zeros(2,2);
    Contact_SL = 0.0;
end
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

function t = pointProjection(A, BC)
% POINTPROJECTION 计算点A到线段BC的投影及归一化参数
%   A: [1,2] 点坐标 [x, y]
%   BC: [2,2] 线段端点，第一行是B，第二行是C
%   
%   输出:
%   proj_point: [1,2] 投影点坐标（可能在线段延长线上）
%   t: 标量，归一化位置参数
%      t=0 投影在B点，t=1 投影在C点
%      t<0 投影在B外侧，t>1 投影在C外侧

    B = BC(1,:);
    C = BC(2,:);
    
    % 线段向量
    v = C - B;              % [1,2] 从B指向C
    
    % 从B指向A的向量
    w = A - B;              % [1,2]
    
    % 计算投影参数 t = (w·v) / (v·v)
    % 这是标量投影长度与线段长度的比值
    v_dot_v = dot(v, v);    % 线段长度平方
    
    if v_dot_v < eps        % 处理B=C的退化情况
        t = 0;
        warning('B和C重合，无法定义线段');
        return;
    end
    
    t = dot(w, v) / v_dot_v;
end