function [mpData,VA,VA_old,fext,uvw,frct,Cal_par] = ...
    AMPLE_Solve_Main_BX_explicit(mesh,mpData,VA_old,...  % 物质点/网格参数
                     fd,Cal_par)

% 【核心-力学求解代码】，内置显式动力求解器模块、刚度矩阵拼接、本构计算 该计算代码在
% 源代码基础上重新梳理整理获得，基本全部重新编译
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行力学求解，涵盖物质点本构计算、整体刚度、质量、阻尼矩阵组装、有限元背景网格求解、不
% 平衡力迭代求解全过程，为物质点力学求解核心子函数
%
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

[fint,mpData,mesh.Stifiness] = detMPs_BX_explicit(VA_old.dudt, Cal_par.dt, mpData,Cal_par.c_par,Cal_par);
[mpData, uvw, VA.dudt, VA.du2dt2, frct] = ExplicitSolve(mpData, fint, fext, mesh.Stifiness.Mt, mesh.Stifiness.Ct, Cal_par.dt, VA_old.dudt, mesh.bc, Cal_par.Calculate_time==0, fd);

FORCE = fext+frct;
fint_all = (fint+mesh.Stifiness.Mt*VA.du2dt2+mesh.Stifiness.Ct*VA.dudt);
Cal_par.fErr = norm(FORCE-fint_all)/norm(FORCE+eps);                            % normalised oobf error

% fprintf(1,'%s: %8.3e, %s: %8.3e, %s: %8.3e \n','NR error',Cal_par.fErr,'max(oobf)',norm(abs(FORCE-fint_all)), 'max(uvw)', norm(abs(uvw)));   % text output to screen (NR error)

end
