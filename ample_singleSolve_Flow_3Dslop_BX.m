%% 【主控程序-湿-力耦合物质点法路基长期性能分析】
% 基于隐式/显式动力求解器的湿-力耦合物质点法主控程序，实现循环荷载-降雨耦合

% 作用下路基-地基体系的变形-渗流-吸力演化全过程模拟
%--------------------------------------------------------------------------
% 作者: FAN Haishan（原作者：William Coombs）
% 修改代码范围 90%，新增模块：碰撞模块、渗流模块、耦合模块
% 重点优化：全向量化，通过高维矩阵实现刚度矩阵一次性求解
% 日期: 29/04/2026
% 描述:
% 进行路基-地基耦合体系的湿-力耦合时程分析，采用算子分裂策略实现：
%   1. 力学步：隐式/显式动力求解，更新应力-应变-位移场（考虑吸力-应力耦合）
%   2. 渗流步：隐式求解 Richards 方程，更新饱和度-吸力-水头场
%   3. 耦合步：将力学变形引起的孔隙变化反馈至渗流参数（渗透系数、水头梯度）
%
% 求解策略：
%   - 时间步长：基础步长 0.1s（lstps=10），自适应调整（dt_min=0.001s）
%   - 力学求解：Newton-Raphson迭代（NRitMax=6），瑞利阻尼 C = 0.002K + 1.5M
%   - 渗流求解：独立迭代回路（NRit_Flow），1.0s后激活湿度模块
%   - 输出控制：每 0.5s 绘图输出，每 20 个输出步自动保存
%
% 荷载配置：
%   - 循环车辆荷载：200kN 四组双轮，0.5s 周期（通过 AMPLE_Solve_Main_BX 内部控制）
%   - 自重荷载：逐步施加（min(1.0, t) 线性增长，1.0s 后满值）
%   - 降雨荷载：通过边界流量 Flow_SIZE 持续入渗
%
% 总分析时长：7200s（2小时），涵盖多次干湿循环与车辆荷载作用
%--------------------------------------------------------------------------
% 程序结构:
% 1. 初始化/断点续算：检测 SAVE_NAME 文件，支持断点续算
% 2. 模型准备：调用 pavement_structure_couple 生成 mpData + mesh
% 3. 参数配置：时间步、阻尼、容差、输出频率等
% 4. 时步循环：while Calculate_time ≤ Time_total
%    4.1 物质点-单元关联更新：elemMPinfo
%    4.2 外力计算：detExtForce（自重 + 接触力 + 车辆荷载）
%    4.3 力学求解：AMPLE_Solve_Main_BX（隐式）或 AMPLE_Solve_Main_BX_explicit（显式）
%    4.4 变形-渗流耦合：Couple_Flow_Stress_BX（孔隙变化 → 水头/渗透系数更新）
%    4.5 渗流求解：AMPLE_Solve_Flow_Main_BX（Richards方程，t≥1.0s激活）
%    4.6 扫尾：时间步调整 SAOWei，物质点更新 updateMPs_BX
%    4.7 结果输出：postPro（Paraview格式），自动保存
%--------------------------------------------------------------------------
% 输入: 无（纯参数化驱动，所有配置内部定义）
%--------------------------------------------------------------------------
% 输出:
%   文件保存: SAVE_NAME.mat（断点续算数据）、SAVE_NAME-60/80/100.mat（关键步备份）
%   图形输出: postPro 生成 Paraview 可视化文件（vtk格式）
%   屏幕输出: 每步计算时间、迭代误差、渗流误差、输出提示
%--------------------------------------------------------------------------
% 关键变量说明:
% Cal_par       - 求解控制参数结构体，含字段：
%   method          - 'dynamic'（动力分析）
%   solve_method    - 'implicit'/'explicit'（隐式/显式求解器切换）
%   c_par           - 瑞利阻尼系数 [A, B] = [0.002, 1.50]，Ct = A*Kt + B*Mt
%   dt              - 当前时间步长 [s]
%   dt_min          - 最小时间步长 0.001s（输出控制精度）
%   dt_SAVE         - 保存时间戳
%   Calculate_time  - 当前计算时刻 [s]
%   Time_total      - 总分析时长 7200s
%   NRitMax         - 最大Newton迭代次数 6
%   NRit/NRit_Flow  - 当前力学/渗流迭代次数
%   fErr/fErr_Flow  - 力学/渗流不平衡力误差
%   tStart/tStart0  - 单步/累计计时器
%   lstp_plot       - 输出步计数器
%   Time_plot       - 上次输出时刻
%   TZ_num          - 预留字段（当前未使用）
%
% mesh          - 背景网格结构体，含字段：
%   coord           - 节点坐标
%   etpl            - 单元拓扑
%   bc/bc_Flow      - 位移/渗流边界条件
%   Stifiness       - 刚度矩阵结构体（Kt, Mt, Ct）
%   H/HVA           - 渗流场（水头/水头梯度/水头变化率）
%   Flow            - 渗流速度场
%   Contact_Force   - 接触力（物质点间相互作用）
%   bounday_imfo    - 边界信息（elemMPinfo更新）
%
% mpData        - 物质点结构体数组，详见 pavement_structure_couple 注释
%
% uvw           - 节点位移向量 [nDoF×1] [m]
% VA            - 速度-加速度结构体，含字段：
%   dudt            - 节点速度 [nDoF×1] [m/s]
%   du2dt2          - 节点加速度 [nDoF×1] [m/s²]
% VA_old        - 上一步速度-加速度结构体
%
% fext          - 节点外力向量 [nDoF×1] [N]
% oobf/oobf_Flow- 力学/渗流出界力（残余不平衡力）
% frct/frct_Flow- 力学/渗流接触力
% fd/fd_Flow    - 自由自由度索引（位移/水头）
% tol           - 收敛容差 1e-6
% lstp          - 总时步计数器
% lstps         - 绘图频率 10（每10步绘图一次，即每1s基础步长绘图）
% Output_time   - 输出时间间隔 0.5s
% CT            - 自动保存计数器（每20个输出步保存一次）
%--------------------------------------------------------------------------
% 求解器切换说明:
% 隐式求解（implicit）:
%   - 适用：准静态或低频动力问题，大时间步稳定
%   - 特点：Newton-Raphson迭代，一致切线模量，二次收敛
%   - 条件：fErr < tol 且 NRit < NRitMax，或 NRit < 2（强制最少2次迭代）
%
% 显式求解（explicit）:
%   - 适用：高频动力或冲击问题，小时间步
%   - 特点：无需求解线性方程组，条件稳定（CFL条件限制）
%   - 调用：AMPLE_Solve_Main_BX_explicit（需单独配置质量矩阵对角化）
%
% 算子分裂策略:
%   力学步 → 变形耦合 → 渗流步，两步解耦但顺序耦合
%   变形对渗流的影响：通过 Couple_Flow_Stress_BX 更新孔隙比e、渗透系数K、水头H
%   渗流对力学的影响：通过吸力 STATEV(:,7) 反馈至本构模型（边界面/Mohr-Coulomb）
%--------------------------------------------------------------------------
% 调用关系:
% 主控调用: pavement_structure_couple（模型生成）、elemMPinfo（MP-单元关联）、
%           detExtForce（外力组装）、detFDoFs/detFDoFs_Flow（自由度识别）、
%           Par_initial/Par_initial_Flow（变量初始化）、
%           AMPLE_Solve_Main_BX/AMPLE_Solve_Main_BX_explicit（力学求解）、
%           Couple_Flow_Stress_BX（变形-渗流耦合）、
%           AMPLE_Solve_Flow_Main_BX（渗流求解）、
%           SAOWei（时间步调整）、updateMPs_BX（物质点更新）、
%           postPro（结果输出）
%--------------------------------------------------------------------------
% 注意事项:
%   - 程序支持断点续算：首次运行生成 SAVE_NAME.mat，后续自动加载续算
%   - 湿度模块（渗流求解）在 t < 1.0s 时跳过，避免初始瞬态不稳定
%   - 自重荷载线性增长：t < 1.0s 时 fext × t，t ≥ 1.0s 时满值
%   - 输出控制：ceil((t-dt_min)/0.5)*0.5 判定，确保 0.5s 整数倍输出
%   - 关键步备份：lstp_plot = 60/80/100 时单独保存，便于后处理对比
%--------------------------------------------------------------------------
clear;clc;warning('off');
SAVE_NAME = 'Slop_3D';

%% 求解方式 及 动力阻尼参数
try
    load (SAVE_NAME)
catch   
    Cal_par.dt = 0.10;
    %% 计算模式
    Cal_par.method = 'dynamic';
    Cal_par.method0 = Cal_par.method;
    Cal_par.c_par = [0.002,1.50];                                             % 阻尼：瑞利阻尼： C=A*K+B*M
    Cal_par.solve_method = 'implicit';
    Cal_par.Explicit = 0;
    Cal_par.Error = false;
    PLOT_FILE = [];
    %% 模型准备
    [mpData,mesh] = Slop_3D;                                                  % setupGrid2D/setupGrid3D/setupGrid_beam/setupGrid_collapse
    [nodes,nD] = size(mesh.coord);                                            % number of nodes and dimensions
    mesh.Contact_Force = zeros(size(mpData,2),nD);

    %% 参数初始化
    Cal_par.dt_SAVE = 0.0; Cal_par.Calculate_time = 0.0;
    Cal_par.dt_min=0.001; Cal_par.JS=0; Cal_par.Time_plot=-1;
    Cal_par.lstp_plot = 0; Cal_par.tStart = tic;                              % zero loadstep counter (for plotting function)
    tol = 1e-6;  CT=0;

    nDoF = nodes*nD;                                                          % total number of degrees of freedom
    lstp = 0;                                                                 % zero loadstep counter (for plotting function)
    % 初始化 & 绘图
    [uvw,VA,mesh.Stifiness,oobf,Cal_par.fErr,frct,Cal_par.NRit] = Par_initial(nDoF);
    [mesh.H,mesh.HVA,mesh.Flow,oobf_Flow,Cal_par.fErr_Flow,frct_Flow,Cal_par.NRit_Flow] = Par_initial_Flow(nodes);
    postPro(mpData,mesh,Cal_par,uvw,VA,SAVE_NAME);                            % plotting initial state & mesh

end

%% 加载参数
Cal_par.Gravity_Time = 1.0; Cal_par.Load_Time = 1.00; Cal_par.Stop_Time = 59.00;
Cal_par.lstps0 = 1/(Cal_par.Gravity_Time/10);Cal_par.lstps1 = 1/(Cal_par.Load_Time/20); Cal_par.lstps2 = 0.10;

Cal_par.tStart0 = tic; Cal_par.tStart = tic; aha=0;
                                       % 显式动力求解器
Time_total = 3600*24;
%% 开始迭代
while Cal_par.Calculate_time<=Time_total                                      % loadstep loop

    [uvw,VA,mesh.Stifiness,oobf,Cal_par.fErr,frct,Cal_par.NRit] = Par_initial(nDoF);
    [mesh.H,mesh.HVA,mesh.Flow,oobf_Flow,Cal_par.fErr_Flow,frct_Flow,Cal_par.NRit_Flow] = Par_initial_Flow(nodes);
    
    Time = rem(Cal_par.Calculate_time+1e-4-Cal_par.Gravity_Time,Cal_par.Load_Time + Cal_par.Stop_Time);
    
    % 分配不同的分析步步长
    if strcmp(Cal_par.solve_method, 'implicit')
        if Time<(Cal_par.dt+1e-4) && Time>=-1e-4
            Cal_par.dt = 0.20;
            lstps = Cal_par.lstps2;
            Cal_par.NRitMax = 5;
        elseif Time<-0.999
            Cal_par.dt = 0.10;
            lstps = Cal_par.lstps0;
            Cal_par.NRitMax = 5;
        end
        Cal_par.dt_min=0.001;
    else
        Cal_par.dt = 1e-4;
        lstps = 10000;
        Cal_par.NRitMax = 6;
        Cal_par.dt_min = 1e-4;
    end

    Time = Time + Cal_par.dt;
    % text output to screen (loadstep)
    Cal_par.Calculate_time = Cal_par.Calculate_time + Cal_par.dt;
    fprintf(1,'\n %s %s %4i %s %.4f %s %.4f %s %.4f %s\n',...
        Cal_par.solve_method, '- loadstep ',lstp,'Calculate_time ',Cal_par.Calculate_time,'Cost Time ',toc(Cal_par.tStart),' s/',toc(Cal_par.tStart0)/60,'min');
    Cal_par.tStart = tic;

    % material point - element information
    [mesh,mpData,VA_old,mesh.bounday_imfo] = elemMPinfo(mesh,mpData,Cal_par);

    fext = detExtForce(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force);        % external force calculation (total)
    fext = fext*min(1.0,Cal_par.Calculate_time);                              % current external force value
    fd   = detFDoFs(mesh);                                                    % free degrees of freedom

    fd_Flow   = detFDoFs_Flow(mesh);                                          % free degrees of freedom

    Cal_par.TZ_num = 0.0;                                                     % zero the iteration counter

    %% 进入循环_求解
    if strcmp(Cal_par.solve_method, 'implicit')
        while (Cal_par.fErr > tol) && (Cal_par.NRit < Cal_par.NRitMax) || (Cal_par.NRit < 2)   % global equilibrium loop
            [mesh,mpData,VA,VA_old,fext,oobf,fd,uvw,frct,Cal_par,~] = ...
                AMPLE_Solve_Main_BX(mesh,mpData,VA_old,...  % 物质点/网格参数
                fext,oobf,fd,...        % 外力数据
                uvw,frct,...            % 核心出装
                Cal_par,nDoF,lstp);     % 计算参数
        end
    elseif strcmp(Cal_par.solve_method, 'explicit')
        [mpData,VA,VA_old,fext,uvw,frct,Cal_par] = ...
            AMPLE_Solve_Main_BX_explicit(mesh,mpData,VA_old,...  % 物质点/网格参数
            fd,Cal_par);     % 计算参数
    else
        error(' Method must be ''explicit'' or ''explicit'' !!')
    end

    %% 隐式自动切换
    % 如隐式迭代误差>0.001或隐式算崩了，切显式
    if strcmp(Cal_par.solve_method, 'implicit') && ...
             (Cal_par.fErr>1e-2 || Cal_par.Error) && ...
              (Cal_par.Calculate_time+1e-5 - Cal_par.dt) >= Cal_par.Gravity_Time
        Cal_par.Calculate_time = Cal_par.Calculate_time - Cal_par.dt;
        Cal_par.solve_method = 'explicit';
        Cal_par.Explicit = 0;
        disp (['    === >>> error: ',num2str(Cal_par.fErr),', dt: ',num2str(Cal_par.dt),'切显式'])
        CT = 0;
        continue
    end
    % 如隐式迭代速度普遍很低，连续50步不出现高速物质点，切隐式
    if strcmp(Cal_par.solve_method, 'explicit') 
        if max(abs(VA_old.dudt))<0.1
            Cal_par.Explicit = Cal_par.Explicit + 1;
        else
            Cal_par.Explicit = 0;
        end
        if Cal_par.Explicit >=50
            Cal_par.Calculate_time = Cal_par.Calculate_time - Cal_par.dt;
            Cal_par.solve_method = 'implicit';
            Cal_par.Error = false;
            Cal_par.dt = 0.20;
            disp ('    === >>> 连续50代速度较低，切隐式')
            continue
        end
    end
    Time = rem(Cal_par.Calculate_time+1e-4-Cal_par.Gravity_Time,Cal_par.Load_Time + Cal_par.Stop_Time);
    disp (['    === >>> 变形计算已完成，error: ',num2str(Cal_par.fErr)])

    %% 根据变形更新H
    [mesh, mpData] = Couple_Flow_Stress_BX(mesh, mpData, Cal_par);
    disp ('    === >>> 已将变形导致的孔隙变化引起的水头影响更新至渗流参数（包括渗透系数、水头、水头梯度）')
    if Cal_par.Calculate_time>=1.0
        % 湿度模块
        while (Cal_par.fErr_Flow > 1e-3*tol) && (Cal_par.NRit_Flow < Cal_par.NRitMax) || (Cal_par.NRit_Flow < 3)   % global equilibrium loop
            [mesh,mpData,mesh.HVA,mesh.HVA_old,oobf_Flow,fd_Flow,mesh.H,frct_Flow,Cal_par] = ...
                AMPLE_Solve_Flow_Main_BX(mesh,mpData,mesh.HVA_old,...         % 物质点/网格参数
                oobf_Flow,fd_Flow,...                                         % 外力数据
                mesh.H,frct_Flow,Cal_par);                                    % 计算参数
        end
        disp (['    === >>> 渗流计算已完成，error: ',num2str(Cal_par.fErr_Flow)])
    end

    %% 扫尾工作
    Cal_par = SAOWei(Cal_par,lstps);                                          % 时间步调整，特定时间步输出
    lstp = lstp+1;
    mpData = updateMPs_BX(uvw,VA.dudt,VA.du2dt2,mpData,mesh);                 % update material points

    %% 结果输出
    % 分配不同的分析步步长
    if strcmp(Cal_par.solve_method, 'implicit')
        if (abs(Cal_par.Calculate_time-Cal_par.Gravity_Time) <= 1.1*Cal_par.dt_min ||...
                abs(rem(Time+0.1*Cal_par.dt_min, 10.0)) <= 1.1*Cal_par.dt_min ||...
                abs(Time) <= 1.1*Cal_par.dt_min) && ...
                (Cal_par.Calculate_time - Cal_par.Time_plot)>100*Cal_par.dt_min
            Cal_par.lstp_plot = Cal_par.lstp_plot+1;
            Cal_par.Time_plot = Cal_par.Calculate_time;
            postPro(mpData,mesh,Cal_par,uvw,VA,SAVE_NAME);                        % Plotting and post processing
            disp(['  ==>>> 已输出计算结果文件   Time = ',num2str(Cal_par.Calculate_time),'  lstp_plot=',num2str(Cal_par.lstp_plot)])
            PLOT_FILE(end+1,1) = Cal_par.Calculate_time;
            CT = CT+1;
            if CT>=6
                CT = 0;
                save(SAVE_NAME,'-V7')
                disp('  ==>>> 已保存MAT文件！！！')
            end
        end
    else
        CT = CT+1;
        if CT>=10
            CT = 0;
            save(SAVE_NAME,'-V7')
            Cal_par.lstp_plot = Cal_par.lstp_plot+1;
            Cal_par.Time_plot = Cal_par.Calculate_time;
            postPro(mpData,mesh,Cal_par,uvw,VA,SAVE_NAME);                        % Plotting and post processing
            PLOT_FILE(end+1,1) = Cal_par.Calculate_time;
            disp(['  ==>>> 已输出计算结果文件   Time = ',num2str(Cal_par.Calculate_time),'  lstp_plot=',num2str(Cal_par.lstp_plot)])
            disp('  ==>>> 已保存MAT文件！！！')
        end
    end
end