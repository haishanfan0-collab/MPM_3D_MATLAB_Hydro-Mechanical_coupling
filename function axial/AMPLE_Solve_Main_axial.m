%% 核心求解代码
function [mesh,mpData,VA,VA_old,oobf,fd,uvw,frct,Cal_par,fext] = ...
    AMPLE_Solve_Main_axial(mesh,mpData,VA_old,...  % 物质点/网格参数
    oobf,fd,...        % 外力数据
    uvw,frct,...            % 核心出装
    Cal_par,nDoF,lstp,fext)
if Cal_par.Calculate_time > 1.0
    Cal_par.method = 'dynamic';
else
    Cal_par.method = 'static';
end
[nodes,nD] = size(mesh.coord);
[duvw,VA.dudt,VA.du2dt2,drct] = model_Solve(mesh.bc, mesh.Stifiness.Kt, mesh.Stifiness.Mt, mesh.Stifiness.Ct, oobf,...
    Cal_par.NRit, fd, Cal_par.dt,...
    uvw,VA_old.dudt,VA_old.du2dt2,Cal_par.method);    % linear solver
uvw  = uvw+duvw;                                                                        % update displacements
frct = frct+drct;                                                                       % update reaction forces
[fint,mpData,mesh.Stifiness] = detMPs_axial(uvw,mpData,Cal_par.c_par,Cal_par);          % global stiffness & internal force
if Cal_par.NRit == 0
    if Cal_par.Calculate_time<=1.0
        fext = detExtForce_axial(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force,0.0);            % external force calculation (total)                                               % current external force value
        fext = min(1.0,Cal_par.Calculate_time)*fext;
    else
        fext_BS = 0.5*(1+cos(2*pi*Cal_par.Calculate_time/1.0+pi));
        fext = detExtForce_axial(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force,fext_BS);            % external force calculation (total)                                               % current external force value

    end
end
Cal_par.NRit = Cal_par.NRit+1;                                                          % increment the NR counter
FORCE = fext+frct;
fint_all = (fint+mesh.Stifiness.Mt*VA.du2dt2+mesh.Stifiness.Ct*VA.dudt);
Cal_par.fErr = norm(FORCE-fint_all)/norm(FORCE+eps);                                    % normalised oobf error
%% 动态弹性求解器
if (Cal_par.fErr<0.50 || Cal_par.NRit<10) && (Cal_par.fErr<10 || Cal_par.NRit<3) && (norm(abs(uvw))<1e6)
    if Cal_par.NRit>5
        oobf = oobf*0.5+(FORCE-fint_all)*0.5;
    else
        oobf = (FORCE-fint_all);
    end
else
    % 缩小时间步
    if Cal_par.dt>Cal_par.dt_min
        Cal_par.Calculate_time = Cal_par.Calculate_time - Cal_par.dt;
        Cal_par.dt=max(Cal_par.dt_min,Cal_par.dt/4);
        Cal_par.Calculate_time = Cal_par.Calculate_time + Cal_par.dt;
        fext = detExtForce_axial(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force);
        fext = fext*min(1.0,Cal_par.Calculate_time);
        [uvw,VA,mesh.Stifiness,oobf,Cal_par.fErr,frct,Cal_par.NRit] = Par_initial(nDoF);
        disp(['  ==>>> 部分物质点速度/加速度/变形过大，已缩小时间步   dt = ',num2str(Cal_par.dt)])
        fprintf(1,'\n%s %4i %s %.4f %s %.4f %s \n','loadstep ',lstp,'Calculate_time ',Cal_par.Calculate_time,'Cost Time ',toc(Cal_par.tStart)/60,' min');             % text output to screen (loadstep)
        % 极端物质点位置调整/删除操作
    else
        Cal_par.TZ_num = Cal_par.TZ_num+1;
        % 保障函数，删除速度过快物质点，认为已经飞出模型
        if Cal_par.TZ_num<3
            [mesh,mpData,VA_old,Cal_par.NRit,uvw,frct,fd,mesh.bounday_imfo] = updataLS_mp_axial(VA.dudt,VA.du2dt2,mpData,mesh,Cal_par.method,uvw);
            fext = detExtForce_axial(nodes,nD,mesh.g,mpData,'y',mesh.Contact_Force);
            fext = fext*min(1.0,Cal_par.Calculate_time);
        else
            Cal_par.NRit=1000;
        end
    end
end
fprintf(1,'%s %2i %s: %8.3e, %s: %8.3e, %s: %8.3e \n','  iteration',Cal_par.NRit,' NR error',Cal_par.fErr,'max(oobf)',norm(abs(FORCE-fint_all)), 'max(uvw)', norm(abs(uvw)));   % text output to screen (NR error)
end



%% 求解主程序，自动接入"动力"和"静力"模块
function [duvw,dudt,du2dt2,drct] = model_Solve(bc, Kt, Mt, Ct, oobf, NRit, fd, dt, uvw, dudt_old,du2dt2_old,method)

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