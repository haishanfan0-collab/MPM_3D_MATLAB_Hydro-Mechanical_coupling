function [mesh,mpData,VA,VA_old,fext,oobf,fd,uvw,frct,Cal_par,duvw] = ...
    AMPLE_Solve_Main_BX(mesh,mpData,VA_old,...  % 物质点/网格参数
                     fext,oobf,fd,...        % 外力数据
                     uvw,frct,...            % 核心出装
                     Cal_par,nDoF,lstp)
% 【核心-力学求解代码】，内置隐式动力/静力求解器模块、刚度矩阵拼接、本构计算 该计算代码在
% 源代码基础上重新梳理整理获得，基本全部重新编译
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行力学求解，涵盖物质点本构计算、整体刚度、质量、阻尼矩阵组装、有限元背景网格求解、不
% 平衡力迭代求解全过程，为物质点力学求解核心子函数
%
%--------------------------------------------------------------------------
% [mesh,mpData,VA,VA_old,fext,oobf,fd,uvw,frct,Cal_par,duvw] = 
%    AMPLE_Solve_Main_BX(mesh,mpData,VA_old, 
%                     fext,oobf,fd,        
%                     uvw,frct,           
%                     Cal_par,nDoF,lstp)
%--------------------------------------------------------------------------
% 输入:
% mesh     - 网格信息结构体参数
% mpData   - 物质点信息结构体参数
% VA_old   - 节点速度、加速度信息
% fext     - 节点外力矩阵
% oobf     - 不平衡力矩阵
% fd       - 需要求解的节点序列
% uvw      - 背景网格位移序列
% frct     - 节点内力矩阵
% Cal_par  - 求解参数结构体参数
% nDoF     - 求解维度，节点维度*节点数量
% lstp     - 求解步
%--------------------------------------------------------------------------
% 输出:
% mesh     - 网格信息结构体参数
% mpData   - 物质点信息结构体参数
% VA       - 节点速度、加速度（新）
% VA_old   - 节点速度、加速度（旧）
% fext     - 节点外力矩阵
% oobf     - 不平衡力矩阵
% fd       - 需要求解的节点序列
% uvw      - 背景网格位移序列
% frct     - 节点内力矩阵
% Cal_par  - 求解参数结构体参数
% duvw     - 该步骤计算的背景网格位移增量
%--------------------------------------------------------------------------
% 此子函数包含子函数
% detExtForce - 节点内力计算
% model_Solve - 求解器，包含
%                Static_Solve  - 隐式静力求解模块
%                Dynamic_Solve - 隐式动力求解模块
% detMPs_BX - 物质点本构模型计算及刚度、阻尼、质量矩阵拼接（向量版），包含
%                所有本构模型（线弹性、弹塑性） 同样为【核心】子函数
% updataLS_mp - 兜底函数，强行更新物质点位置，避免计算奇异
%--------------------------------------------------------------------------
[nodes,nD] = size(mesh.coord);
%% 荷载计算
if Cal_par.NRit == 0
    Time = rem(Cal_par.Calculate_time-Cal_par.Gravity_Time,Cal_par.Load_Time + Cal_par.Stop_Time);
    if Time<=0.0
        fext = detExtForce(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force,0.0);            % external force calculation (total)                                               % current external force value
        fext = min(1.0,Time+Cal_par.Gravity_Time)*fext;
    else
        fext_BS = zeros(nD,1);
        if Time<=Cal_par.Load_Time
            fext_BS(2) = 0.5*(1+cos(2*pi*Time/Cal_par.Load_Time+pi));
        else
            fext_BS(2) = 0;
        end
        fext = detExtForce(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force,fext_BS);            % external force calculation (total)                                               % current external force value
    end
    % % 每个循环加载开始将加速度清零
    % if rem(Cal_par.Calculate_time+0.0001-Cal_par.dt,1.0)<0.001
    %     VA_old.dudt(:) = 0.0;
    %     VA_old.du2dt2(:) = 0.0;
    % end
end

if Cal_par.Calculate_time > Cal_par.Gravity_Time
    Cal_par.method = Cal_par.method0;
else
    Cal_par.method = 'dynamic';
end

[duvw,VA.dudt,VA.du2dt2,drct] = model_Solve(mesh.bc, mesh.Stifiness.Kt, mesh.Stifiness.Mt, mesh.Stifiness.Ct, oobf,...
    Cal_par.NRit, fd, Cal_par.dt,...
    uvw,VA_old.dudt,VA_old.du2dt2,Cal_par.method);    % linear solver
uvw  = uvw+duvw;                                                                        % update displacements
frct = frct+drct;                                                                       % update reaction forces
[fint,mpData,mesh.Stifiness] = detMPs_BX(uvw,mpData,Cal_par.c_par,Cal_par);                % global stiffness & internal force

Cal_par.NRit = Cal_par.NRit+1;                                                          % increment the NR counter
FORCE = fext+frct;
fint_all = (fint+mesh.Stifiness.Mt*VA.du2dt2+mesh.Stifiness.Ct*VA.dudt);
Cal_par.fErr = norm(FORCE-fint_all)/norm(FORCE+eps);                            % normalised oobf error
%% 动态弹性求解器
if (Cal_par.fErr<0.01 || Cal_par.NRit<5) && (Cal_par.fErr<0.10 || Cal_par.NRit<4) && max(abs(VA_old.dudt))<0.1
    if Cal_par.NRit>6
        oobf = oobf*0.05+(FORCE-fint_all)*0.95;
    else
        oobf = (FORCE-fint_all);
    end
else
    % 缩小时间步
    if Cal_par.dt>Cal_par.dt_min && max(abs(VA_old.dudt))<0.1
        Cal_par.Calculate_time = Cal_par.Calculate_time - Cal_par.dt;
        Cal_par.dt=max(Cal_par.dt_min,Cal_par.dt/4);
        Cal_par.Calculate_time = Cal_par.Calculate_time + Cal_par.dt;
        fext = detExtForce(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force);
        fext = fext*min(1.0,Cal_par.Calculate_time);
        [uvw,VA,mesh.Stifiness,oobf,Cal_par.fErr,frct,Cal_par.NRit] = Par_initial(nDoF);
        disp(['  ==>>> 部分物质点速度/加速度/变形过大，已缩小时间步   dt = ',num2str(Cal_par.dt)])
        fprintf(1,'\n%s %4i %s %.4f %s %.4f %s \n','loadstep ',lstp,'Calculate_time ',Cal_par.Calculate_time,'Cost Time ',toc(Cal_par.tStart)/60,' min');             % text output to screen (loadstep)
        Cal_par.Error = false;
        % 极端物质点位置调整/删除操作
    else
        % Cal_par.TZ_num = Cal_par.TZ_num+1;
        % % 保障函数，删除速度过快物质点，认为已经飞出模型
        % if Cal_par.TZ_num<3
        %     [mesh,mpData,VA_old,Cal_par.NRit,uvw,frct,fd,mesh.bounday_imfo] = updataLS_mp(VA.dudt,VA.du2dt2,mpData,mesh,Cal_par,uvw);
        %     fext = detExtForce(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force);
        %     fext = fext*min(1.0,Cal_par.Calculate_time);
        % else
        Cal_par.NRit=1000;
        Cal_par.Error = true;
        % end
    end
end
fprintf(1,'%s %2i %s: %8.3e, %s: %8.3e, %s: %8.3e \n','  iteration',Cal_par.NRit,' NR error',Cal_par.fErr,'max(oobf)',norm(abs(FORCE-fint_all)), 'max(uvw)', norm(abs(uvw)));   % text output to screen (NR error)
end



%% 求解主程序，自动接入"动力"和"静力"模块
function [duvw,dudt,du2dt2,drct] = model_Solve(bc, Kt, Mt, Ct, oobf, NRit, fd, dt, uvw, dudt_old,du2dt2_old,method)
% 【核心-力学求解分发代码】，内置静力/动力求解自动路由模块，根据分析类型分发至对应求解器
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行物质点法力学求解路由控制，根据 method 参数自动分发至隐式静力求解或隐式动力求解模块。
% 静力分析调用 Static_Solve 直接求解位移增量与反力；动力分析调用 Dynamic_Solve 基于
% Newmark-β 平均加速度法求解位移增量、速度、加速度及反力。为 AMPLE_Solve_Main_BX 核心
% 子函数的统一求解入口。
%
%--------------------------------------------------------------------------
% [duvw, dudt, du2dt2, drct] = 
%    model_Solve(bc, Kt, Mt, Ct, oobf, NRit, fd, dt, uvw, dudt_old, du2dt2_old, method)
%--------------------------------------------------------------------------
% 输入:
% bc         - 位移边界条件矩阵 [nbc×2]，第1列为约束自由度编号，第2列为边界值
% Kt         - 整体切线刚度矩阵 [nDoF×nDoF]（稀疏）
% Mt         - 整体质量矩阵 [nDoF×nDoF]（稀疏，动力分析时传入）
% Ct         - 整体阻尼矩阵 [nDoF×nDoF]（稀疏，动力分析时传入）
% oobf       - 节点不平衡力向量 [nDoF×1]
% NRit       - Newton-Raphson 迭代计数器（从 0 开始）
% fd         - 自由自由度索引向量 [nfree×1]
% dt         - 时间步长 [标量]（动力分析时传入）
% uvw        - 当前时间步节点位移向量 [nDoF×1]（动力分析时传入）
% dudt_old   - 上一时间步节点速度向量 [nDoF×1]（动力分析时传入）
% du2dt2_old - 上一时间步节点加速度向量 [nDoF×1]（动力分析时传入）
% method     - 求解方法字符串，'static'（静力）或 'dynamic'（动力）
%--------------------------------------------------------------------------
% 输出:
% duvw       - 节点位移增量向量 [nDoF×1]
% dudt       - 更新后的节点速度向量 [nDoF×1]（静力时为零向量）
% du2dt2     - 更新后的节点加速度向量 [nDoF×1]（静力时为零向量）
% drct       - 节点约束反力向量 [nDoF×1]（仅边界自由度非零）
%--------------------------------------------------------------------------
% 此子函数包含子函数
% Static_Solve  - 隐式静力线性方程组求解模块，直接法求解位移增量与反力
% Dynamic_Solve - 隐式动力线性方程组求解模块，Newmark-β 平均加速度法求解
%                 位移增量、速度、加速度及反力
%--------------------------------------------------------------------------
if strcmp(method, 'static')
    [duvw,drct] = Static_Solve(bc,Kt,oobf,NRit,fd);
    nDoF = length(oobf);
    dudt = zeros(nDoF,1);                                                       % zero displacement increment
    du2dt2 = zeros(nDoF,1);                                                       % zero displacement increment
elseif strcmp(method, 'dynamic')
    [duvw,dudt,du2dt2,drct] = Dynamic_Solve(bc,Kt,Mt,Ct,oobf,NRit,fd,uvw,dudt_old,du2dt2_old,dt);
else
    error('method only can be ''dynamic'' and ''static'' !!!')
end

end