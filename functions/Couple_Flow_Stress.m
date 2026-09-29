function [mesh, mpData] = Couple_Flow_Stress(mesh, mpData, Cal_par)

%% 该程序用于将应变导致的水头变化反馈到吸力计算上
nmp = length(mpData);                                                       % number of material points
nD  = length(mpData(1).mpC);
t = [1 5 9];
for mp=1:nmp
    if mpData(mp).modulus>0.0 && mpData(mp).Gravity_Water_Content>0.0
        if mpData(mp).STATEV(7)~=0.0
            % 更新饱和度
            F   = mpData(mp).F;                                                     % deformation gradient
            mpData(mp).vp    = det(F)*mpData(mp).vp0;                               % update material point volumes

            Fa = mpData(mp).Flow_CP(1)/1000;
            Fm = mpData(mp).Flow_CP(2);
            Fn = mpData(mp).Flow_CP(3);
            % mpData(mp).e = (1+mpData(mp).e0)*mpData(mp).vp/mpData(mp).vp0-1;        % 孔隙比
            % 方法1：使用中间变量（推荐，清晰易读）
            mpData(mp).e = (1+mpData(mp).e0)*mpData(mp).vp/mpData(mp).vp0-1;        % 孔隙比
            e = mpData(mp).e;
            coeff = Fa / (e^(1/(Fm*Fn)));  % 对应公式中的 m3/e^{1/(m1m2)}
            mpData(mp).Sr = (1 + (mpData(mp).Scution_Flow ./ coeff).^Fn).^(-Fm);

            if mpData(mp).Sr>=1.00
                mpData(mp).Sr = 1.00;
            end

            % 更新吸力
            dH = (mpData(mp).STATEV(7) - mpData(mp).STATEVn(7))/-9.81;
            
        else
            % 更新饱和度
            F   = mpData(mp).F;                                                     % deformation gradient
            mpData(mp).vp    = det(F)*mpData(mp).vp0;                               % update material point volumes

            e_old = mpData(mp).e;
            Fa = mpData(mp).Flow_CP(1)/1000;
            Fm = mpData(mp).Flow_CP(2);
            Fn = mpData(mp).Flow_CP(3);
            mpData(mp).e = (1+mpData(mp).e0)*mpData(mp).vp/mpData(mp).vp0-1;        % 孔隙比
            e = mpData(mp).e;
            mpData(mp).Sr =   mpData(mp).Sr*(e_old/mpData(mp).e);                   % 由于降雨导致水头变化产生的吸力

            if mpData(mp).Sr>=1.00
                mpData(mp).Sr = 1.00;
            end

            % 更新吸力
            Scution_Flow_old = mpData(mp).Scution_Flow;
            mpData(mp).Scution_Flow = ((mpData(mp).Sr^(-1/Fm) - 1)^(1/Fn))*Fa/mpData(mp).e^(1/(Fm*Fn));

            dH = (mpData(mp).Scution_Flow - Scution_Flow_old)/-9.81;
        end

        % 更新吸力变化 (由于孔隙率变化导致的）
        mpData(mp).dH_sum_strain =  mpData(mp).dH_sum_strain + dH;

        % 更新物质点 dH (由于孔隙率变化导致的）
        mpData(mp).dH = mpData(mp).dH + dH/Cal_par.dt;
        
        if mpData(mp).Flow_Type == 99
            mpData(mp).dH_sum_strain = 0;
            mpData(mp).dH = 0;
            mpData(mp).Sr = 1.00;
            mpData(mp).Scution_Flow = 0.0;
        end

        % 更新干密度
        mpData(mp).Roud = mpData(mp).Gs/(1+e);

        % 计算质量含水率
        mpData(mp).Gravity_Water_Content = (mpData(mp).Sr .* e) ./ mpData(mp).Gs*1000;
        mpData(mp).Gravity_Water_Content_old = mpData(mp).Gravity_Water_Content;

        % 更新物质点质量
        mpData(mp).mpM = mpData(mp).vp*(1+mpData(mp).Gravity_Water_Content)*mpData(mp).Roud;

        % 更新 lp
        if mpData(mp).mpType == 2                                               % GIMPM only (update domain lengths)
            [V,D] = eig(F.'*F);                                                 % eigen values and vectors F'F
            U     = V*sqrt(D)*V.';                                              % material stretch matrix
            mpData(mp).lp = (mpData(mp).lp0).*U(t(1:nD));                       % update domain lengths
        end

        % 更新渗透系数（颗粒拉伸压缩）
        mpData(mp).Ksat(1) = mpData(mp).lp(2)/mpData(mp).lp0(2)*mpData(mp).Ksat0(1);
        mpData(mp).Ksat(2) = mpData(mp).lp(1)/mpData(mp).lp0(1)*mpData(mp).Ksat0(2);
    end

end

% 插值更新背景网格吸力水头
if Cal_par.Calculate_time>=1.0
    [mesh,mpData,~,~] = elemMPinfo(mesh,mpData,Cal_par.method);
end