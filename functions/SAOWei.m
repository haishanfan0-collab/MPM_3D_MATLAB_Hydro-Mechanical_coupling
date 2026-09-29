%% 扫尾子函数
function Cal_par = SAOWei(Cal_par,lstps)

if Cal_par.fErr<1e-2 && Cal_par.dt<1/lstps
    Cal_par.JS = Cal_par.JS+1;
    if Cal_par.JS >1
        Cal_par.JS=0;
        Cal_par.dt=min(1/lstps,max(Cal_par.dt_SAVE*1.5,Cal_par.dt*1.5));
        % disp(['  ==>>> 连续两个分析步收敛不错，已适当放大时间步   dt = ',num2str(Cal_par.dt)])
    end
end
Cal_par.dt = max(Cal_par.dt, Cal_par.dt_SAVE);
Cal_par.dt_SAVE = 0.0;

if (Cal_par.Calculate_time+Cal_par.dt-(1/lstps*ceil(1e-5+(Cal_par.Calculate_time-Cal_par.Gravity_Time)*lstps)+Cal_par.Gravity_Time))>Cal_par.dt_min 
    Cal_par.dt_SAVE = Cal_par.dt;
    Cal_par.dt=max(Cal_par.dt_min,    (1/lstps*ceil(1e-5+(Cal_par.Calculate_time-Cal_par.Gravity_Time)*lstps)+Cal_par.Gravity_Time)-Cal_par.Calculate_time);
    % disp(['  ==>>> 为了整数输出   dt = ',num2str(Cal_par.dt)])
end

end