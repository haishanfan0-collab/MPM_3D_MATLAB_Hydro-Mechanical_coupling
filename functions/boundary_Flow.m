function mpData = boundary_Flow(mpData,Flow_boundary,Flow_P,TZ_Position)
% 【核心-渗流边界条件配置与边界物质点几何修正】
% 基于边界虚拟物质点的渗流边界自动识别、流量分配与边坡几何修正，
% 实现降雨入渗、底面排水、坡面渗流等多类型边界条件的物质点法配置
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 29/04/2026
% 描述:
% 进行渗流边界条件的自动识别与配置，核心功能包括：
%   1. 边界类型识别：从虚拟物质点（mpType=3）中自动分类顶面、底面、
%      垂直面、边坡面四类几何边界
%   2. 边坡几何修正：对边坡面虚拟点进行线性回归，修正为严格共线排列，
%      消除边界锯齿效应，提升渗流计算稳定性
%   3. 流量等效分配：将面流量强度 [m/s] 转化为边界物质点等效节点流量 [m³/s]，
%      采用中点法积分（端点半宽、内点全宽）
%   4. 边界标记写入：为边界物质点写入Flow_Position、Flow_P、Flow_SIZE字段
%
% 边界处理策略：
%   - 顶面（Up）：降雨入渗边界，Flow_P > 0 为入渗，Flow_P < 0 为蒸发
%   - 底面（Down）：排水边界，通常 Flow_P = 0（零水头）或固定流量
%   - 垂直面（Vertical）：侧向边界，TZ_Position指定'Left'/'Right'
%   - 边坡面（Slope）：倾斜自由面，需几何修正为严格直线
%
% 边坡修正算法：
%   提取边坡面虚拟点后，删除水平共线冗余点（dy=0），对剩余点进行
%   最小二乘线性拟合：y = kx + b，强制将所有边坡点投影至该直线，
%   确保边坡面严格共线，避免数值振荡
%--------------------------------------------------------------------------
% mpData = boundary_Flow(mpData, Flow_boundary, Flow_P, TZ_Position)
%--------------------------------------------------------------------------
% 输入:
% mpData        - 物质点信息结构体数组 [1×mpM]（含内部点+虚拟边界点）
% Flow_boundary - 边界类型字符串胞数组，如 {'Slope','Up'} 或 {'Up'}
%                 可选值: 'Up'（顶面）、'Down'（底面）、
%                        'Slope'（边坡面）、'Vertical'（垂直面）
% Flow_P        - 边界流量强度 [m/s]（单位面积入渗/渗出速率）
%                 物理意义: 降雨强度/渗透系数比，或达西流速
%                 等效节点流量: Flow_P × 控制宽度 × 0.5（端点）或 × 1.0（中点）
% TZ_Position   - 垂直面方位标识字符串: 'Right'（右侧垂直面）或 'Left'（左侧垂直面）
%                 用于确定垂直面法向方向，配合边坡面识别
%--------------------------------------------------------------------------
% 输出:
% mpData        - 更新后的物质点结构体数组，边界物质点新增/更新字段：
%   边界识别与修正字段:
%     Flow_Position  - 边界位置标记字符串: 'Up'/'Down'/'Slope'/'Vertical'
%     Flow_SIZE      - 边界流量强度 [m/s]（与输入Flow_P一致，用于诊断输出）
%
%   等效节点流量字段（核心）:
%     Flow_P         - 等效节点流量 [m³/s] 或 [m²/s]（二维简化）
%                      计算方式: Flow_P × 相邻点间距 × 权重系数
%                      端点权重: 0.5（半宽贡献）
%                      内点权重: |x_{j+1} - x_{j-1}| × 0.5（全宽贡献）
%
%   几何修正字段（边坡面）:
%     mpC            - 修正后的坐标 [x,y]，严格位于拟合直线 y = kx + b 上
%
% 删除操作:
%   水平共线冗余边坡点（LS索引）将被删除，避免边界重叠导致的数值奇异
%--------------------------------------------------------------------------
% 边界识别算法细节:
% 1. 提取虚拟点: 从mpData中筛选mpType=3的边界物质点，构建坐标矩阵
%    mpC格式: [全局索引, mpType, x, y]
%
% 2. 顶面识别: y坐标最大值点群（容差0.001），按x排序
%    mpC_Up = mpC(y ≥ max(y)-0.001, :)
%
% 3. 底面识别: y坐标最小值点群（容差0.001）
%    mpC_Down = mpC(y ≤ min(y)+0.001, :)
%
% 4. 垂直面识别: 
%    TZ_Position='Right': x坐标最小值点群（右边界，x ≤ min(x)+0.001）
%    TZ_Position='Left':  x坐标最大值点群（左边界，x ≥ max(x)-0.001）
%    按y排序后提取
%
% 5. 边坡面识别: 剩余点即为边坡面候选点，删除水平共线点（dy=0）后，
%    进行线性拟合: k = (y_max - y_min)/(x_max - x_min)
%                b = (y_max·x_min - y_min·x_max)/(x_min - x_max)
%    强制投影: x_new = (y - b)/k（保持y不变，修正x至拟合直线）
%
% 6. 冗余点删除: 水平共线点（dy=0）视为内部过渡点，从mpData中删除
%--------------------------------------------------------------------------
% 流量分配公式:
% 对于排序后的边界点列（按边界走向排序，如顶面按x排序）:
%   首点:  Flow_P_point = Flow_P × |x₂ - x₁| × 0.5
%   末点:  Flow_P_point = Flow_P × |xₙ - xₙ₋₁| × 0.5
%   内点:  Flow_P_point = Flow_P × |xⱼ₊₁ - xⱼ₋₁| × 0.5
% 物理意义: 中点积分法则，将面流量转化为节点等效流量
%
% 注意: 二维问题中流量单位为 [m²/s]（厚度方向取单位1），
%       实际三维需乘以边界面积
%--------------------------------------------------------------------------
% 调用关系:
% 被调用: pavement_structure_couple（路基建模主函数，用于渗流边界配置）
% 内部调用: 无（纯MATLAB基础函数：sortrows, ismember, strcmp等）
%--------------------------------------------------------------------------
% 注意事项:
%   - 边坡修正仅修正x坐标，保持y坐标不变（假设边坡为y关于x的函数）
%   - 水平共线点删除后，相邻点间距重新计算，确保流量分配连续
%   - 若边坡面存在垂直段（k→∞），需特殊处理（当前代码假设边坡为有限斜率）
%   - 多边界类型同时配置时（如{'Slope','Up'}），各边界独立处理，无耦合效应
%--------------------------------------------------------------------------

%% 第一步边界物质点位置识别
mpM = size(mpData,2);
mpC = [(1:mpM)',[mpData.mpType]',(reshape([mpData.mpC]',2,mpM))'];
mpC(mpC(:,2)~=3,:)=[];

% 顶面查找
mpC_Up = mpC(mpC(:,end) >= max(mpC(:,end)-0.001),:);
mpC_Up = sortrows(mpC_Up,3);
mpC(ismember(mpC(:,1),mpC_Up(:,1)),:) = [];

% 底部查找
mpC_Down = mpC(mpC(:,end) <= min(mpC(:,end)+0.001),:);
mpC(ismember(mpC(:,1),mpC_Down(:,1)),:) = [];
mpC_Down = sortrows(mpC_Down,3);

% 垂直面查找
if strcmp(TZ_Position, 'Right')
    mpC_Vertical = mpC(mpC(:,end-1) <= min(mpC(:,end-1)+0.001),:);
    mpC_Vertical = sortrows(mpC_Vertical,4);
    mpC(ismember(mpC(:,1),mpC_Vertical(:,1)),:) = [];
    mpC_Slope = mpC; clear mpC
elseif strcmp(TZ_Position, 'Left')
    mpC_Vertical = mpC(mpC(:,end-1) >= max(mpC(:,end-1)-0.001),:);
    mpC_Vertical = sortrows(mpC_Vertical,4);
    mpC(ismember(mpC(:,1),mpC_Vertical(:,1)),:) = [];
    mpC_Slope = mpC; clear mpC
else
    error('边坡方向指定错误')
end

% 边坡面修正
mpC_Slope = sortrows(mpC_Slope,4);
LS = mpC_Slope((mpC_Slope(2:end-1,end) - mpC_Slope(1:end-2,end)) == 0,1);
mpC_Slope((mpC_Slope(2:end-1,end) - mpC_Slope(1:end-2,end)) == 0,:) = [];

[x1,N] = max(mpC_Slope(:,end-1));
y1 = mpC_Slope(N,end);

[x2,N] = min(mpC_Slope(:,end-1));
y2 = mpC_Slope(N,end);

k = (y2-y1)/(x2-x1);
b = (y2*x1-y1*x2)/(x1-x2);

mpC_Slope(:,end-1) = ((mpC_Slope(:,end) - b) /k);

% 坐标修正
for i=1:size(mpC_Slope,1)
    mpData(mpC_Slope(i,1)).mpC = mpC_Slope(i,3:end);
end

%% 增加流量边界信息
for i = 1:length(Flow_boundary)
    if strcmp(Flow_boundary(i), 'Up')
        mpC_Select = mpC_Up;
    elseif strcmp(Flow_boundary(i), 'Down')
        mpC_Select = mpC_Down;
    elseif strcmp(Flow_boundary(i), 'Slope')
        mpC_Select = mpC_Slope;
    elseif strcmp(Flow_boundary(i), 'Vertical')
        mpC_Select = mpC_Vertical;
    else
        error(' Flow_boundary 输入错误');
    end
    for j = 1:size(mpC_Select,1)
        if j == 1
            mpData(mpC_Select(j,1)).Flow_P = Flow_P*abs(mpC_Select(j+1,3)-mpC_Select(j,3))*0.50;
        elseif j == size(mpC_Select,1)
            mpData(mpC_Select(j,1)).Flow_P = Flow_P*abs(mpC_Select(j,3)-mpC_Select(j-1,3))*0.50;
        else
            mpData(mpC_Select(j,1)).Flow_P = Flow_P*abs(mpC_Select(j+1,3)-mpC_Select(j-1,3))*0.50;
        end
        mpData(mpC_Select(j,1)).Flow_Position = Flow_boundary(i);
        mpData(mpC_Select(j,1)).Flow_SIZE = Flow_P;
    end
end

mpData(LS) = [];

end