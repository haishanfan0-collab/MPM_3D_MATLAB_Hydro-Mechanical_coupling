function [Dalg,sig,epsE,epsP,Out1] = Boundary_module_BX(eps,mCst,Cal_par,sign,epsn,epsEn,STATEVn,Mc,dMc,STATEV,Pore_Pressure)
% 【核心-非饱和土大周期弹塑性边界面本构模型】
% 基于临界状态土力学框架，考虑吸力效应、状态相关剪胀与模量循环记忆特性
% 适用于路基土在干湿循环-交通荷载耦合作用下的长期变形预测
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 29/04/2026
% 描述:
% 进行非饱和土应力-应变状态更新，涵盖吸力-应力耦合计算、边界面塑性映射、
% 状态相关硬化规律、模量循环记忆效应、干湿循环滞回特性等核心力学过程。
% 本构框架：边界面塑性理论 + 临界状态理论 + 非饱和土力学（吸力-饱和度耦合）
% 特殊功能：1）FHS模量循环记忆（0.5s周期内模量衰减-恢复）；2）干湿状态切换判别；
%           3）吸力相关边界面收缩/扩张；4）塑性流动方向状态依赖性
%
%--------------------------------------------------------------------------
% [Dalg,sig,epsE,epsP,Out1] = Boundary_module_BX(eps,mCst,Cal_par,sign,...
%                                   epsn,epsEn,STATEVn,Mc,dMc,STATEV,Pore_Pressure)
%--------------------------------------------------------------------------
% 输入:
% eps           - 当前步总应变向量 [N×6]，Voigt记法 [ε11,ε22,ε33,γ12,γ13,γ23]
% mCst          - 材料参数矩阵 [N×29]，每行对应一个物质点，列顺序如下：
%                 1:NCL_N, 2:NCL_λ, 3:CSL_Γ, 4:CSL_ω, 5:UnLoad_k, 6:nd, 7:M0
%                 8:Bounding_n, 9:Fabric_nc, 10:Fabric_PAR0, 11:Fabric_PAR1, 12:Fabric_PAR2
%                 13:Bounding_r, 14:pdmax, 15:Gs, 16:K0, 17:ν, 18:FHS_k0
%                 19:FHS_k1, 20:FHS_k2, 21:FHS_k3, 22:FHS_k4
%                 23:UnSAT_a, 24:UnSAT_b, 25:UnSAT_m1, 26:UnSAT_m2, 27:UnSAT_m3D
%                 28:UnSAT_m3W, 29:UnSAT_R
% Cal_par       - 求解控制参数结构体，含字段：
%                 dt        - 时间步长
%                 Calculate_time - 当前计算时刻
%                 NRit      - Newton-Raphson迭代标志（0=仅返回弹性矩阵）
% sign          - 上一步应力向量 [N×6]，Voigt记法（拉为正，压为负）
% epsn          - 上一步总应变 [N×6]
% epsEn         - 上一步弹性应变 [N×6]
% STATEVn       - 上一步状态变量 [N×19]，详见输出说明
% Mc            - 当前步含水质量 [N×1]
% dMc           - 含水质量增量 [N×1]
% STATEV        - 当前步状态变量存储数组 [N×19]（输出容器，部分值继承自STATEVn）
% Pore_Pressure - 孔隙水压力 [N×1]（kPa）
%--------------------------------------------------------------------------
% 输出:
% Dalg          - 一致切线模量张量 [6×6×N]（Voigt矩阵形式）
% sig           - 更新后应力向量 [N×6]
% epsE          - 更新后弹性应变 [N×6]
% epsP          - 更新后塑性应变 [N×6]
% Out1          - 单元胞数组，Out1{1}=EMOD_OLD/1e6（当前模量MPa），Out1{2}=STATEV
%--------------------------------------------------------------------------
% 状态变量 STATEV (19维，单位注意：应力/吸力类存储时除以1000):
%  1: plmax/1000    - 记忆预固结应力（边界面尺寸控制）
%  2: qr/1000       - 参考偏应力（模量计算基准）
%  3: S_r/1000      - 参考吸力（干湿循环滞回基准）
%  4: UnSAT_LS      - 吸力影响系数（NCL/CSL缩放因子）
%  5: UnLoad_k      - 卸载回弹模量指数
%  6: d_s           - 吸力增量
%  7: UnSAT_s/1000  - 当前吸力（kPa）
%  8: Ssa/1000      - 当前干湿状态判别吸力
%  9: Vs_all1       - 累积塑性应变（第一不变量）
% 10: Vs_all2       - 累积塑性应变×剪胀（第二不变量）
% 11: FHS_TIME      - 模量记忆时间戳（用于周期判别）
% 12: FHS_SR/1000   - 吸力-饱和度乘积（模量计算用）
% 13~17: modulus_dw - 0.0~0.5s周期内5个时刻的模量记忆（用于循环恢复）
% 18: σ22_history   - 历史竖向应力（用于后处理-附加竖向应力）
% 19: Stress_Ratio  - 应力比记忆（状态相关硬化参数）
%--------------------------------------------------------------------------
% 关键本构假设与理论框架:
% 1. 边界面方程: F = ln(p/pl) - (Q/M0/p)^Bounding_n·ln(Bounding_r) = 0
% 2. 临界状态线: e = CSL_Γ - CSL_ω·ln(p/PA) 
% 3. 状态相关剪胀: Ds = (Md² - η²)/(2η)，Md = M0·exp(nd·Ψ)
% 4. 吸力效应: 非饱和参数通过UnSAT_LS缩放NCL/CSL位置，影响边界面尺寸
% 5. SWCC模型: 基于Van Genuchten型，考虑主吸湿/脱湿曲线滞回（SsaD/SsaW）
% 6. 模量记忆: FHS模型，0.5s周期内记录5个时刻模量，实现循环荷载下模量衰减-恢复
% 7. 塑性流动: 关联流动法则，边界面内点通过映射距离控制塑性模量Kp
%--------------------------------------------------------------------------
% 此子函数包含/调用子函数
% Calculate_principal  - 批量3×3矩阵特征值分解（主应力计算）
% DDSDDE_Cal           - 弹性刚度矩阵组装（各向同性线弹性）
%--------------------------------------------------------------------------
%% 设定计算常数
N = size(sign,1);
PA = 101325;
FHS_JS = 1;
TEMP1 = 0.80;
TEMP2 = 0.40;TEMP_LS = TEMP2;
% FHS_TPJX = zeros(N,1);
% nc_PAR = zeros(N,1);
DTIME = Cal_par.dt;
% ================== Initial parameters
%       SET ZERO

TIME = Cal_par.Calculate_time - Cal_par.dt;
STRESS_OLD = sign;
STRAN_OLD = epsn;
STRESS = sign;
STRESS0 = sign;
DSTRAN_OLD = eps - epsn;
DSTRAN = DSTRAN_OLD;
HH = zeros(N , 6);

Vs = zeros(N,1);
Ds = zeros(N,1);
%Kp = zeros(N,1);

FF1 = zeros(6,6,N);
ONE = zeros(6,6,N);
DDSDDE_NEW = zeros(6,6,N);
for i = 1:6
    ONE(i,i,:) = 1;
end

% CALCULATE sij
SS1 = (STRESS0(:,1)+STRESS0(:,2)+STRESS0(:,3))/3.0;
STRESS_S(:,1:3) = -(STRESS0(:,1:3) - SS1);  % 或写成：-STRESS0(1:3) + SS1
STRESS_S(:,4:6) = -STRESS0(:,4:6);  % 修正索引！

% ================== READ PROPS PRAMTS ==================
NCL_N      = mCst(:,1);	    % Critical state soil mechanics
NCL_lapta  = mCst(:,2);	    % Critical state soil mechanics
CSL_T      = mCst(:,3);	    % Critical state soil mechanics
CSL_w      = mCst(:,4);	    % Critical state soil mechanics
UnLoad_k   = mCst(:,5);	    % Critical state soil mechanics
nd         = mCst(:,6);     % Critical state soil mechanics
M0         = mCst(:,7);     % Critical state soil mechanics
Bounding_n = mCst(:,8);	    % Parameters for bounding surface
Fabric_nc  = mCst(:,9);	    % Hardening parameter
Fabric_PAR0= mCst(:,10);	% Hardening parameter
Fabric_PAR1= mCst(:,11);	% Hardening parameter
Fabric_PAR2= mCst(:,12);	% Hardening parameter
Bounding_r = mCst(:,13);    % Boundary shape parameter
SUB_pdmax  = mCst(:,14);	% Subgrade soil parameters
SUB_Gs     = mCst(:,15);	% Subgrade soil parameters
SUB_K0     = mCst(:,16);	% Subgrade soil parameters
SUB_u      = mCst(:,17);    % Subgrade soil parameters
FHS_k0     = mCst(:,18);	% Elastic Parameters
FHS_k1     = mCst(:,19);	% Elastic Parameters
FHS_k2     = mCst(:,20);	% Elastic Parameters
FHS_k3     = mCst(:,21);	% Elastic Parameters
FHS_k4     = mCst(:,22);	% Elastic Parameters

% ------- Unsaturation parameters - (NEW PROPS)
UnSAT_a    = mCst(:,23);    % state dependent parameter of NCL and CSL
UnSAT_b    = mCst(:,24);    % state dependent parameter of NCL and CSL
UnSAT_m1   = mCst(:,25);    % SWCC Parameters
UnSAT_m2   = mCst(:,26);    % SWCC Parameters
UnSAT_m3D  = mCst(:,27);    % SWCC Parameters
UnSAT_m3W  = mCst(:,28);    % SWCC Parameters
UnSAT_R    = mCst(:,29);    % f(s) Parameters

% ------- Unsaturation parameters - (NEW PROPS)
UnSAT_Ts = 72.8/1000.0;

% ****************** READ STATEV
plmax = STATEVn(:,1)*1000.0;   % Memory stress: plmax
qr = STATEVn(:,2)*1000.0;      % Reference Stress: qr
S_r = STATEVn(:,3)*1000.0;     % Reference suction: sr
UnSAT_LS_OLD = STATEVn(:,4);
UnLoad_k_OLD = STATEVn(:,5);
ds_OLD = STATEVn(:,6);
UnSAT_s = max(1.0,STATEVn(:,7)*1000.0);
Ssa = STATEVn(:,8)*1000.0;
Vs_all1 = STATEVn(:,9);
Vs_all2 = STATEVn(:,10);
FHS_TIME = STATEVn(:,11);
FHS_SR = STATEVn(:,12)*1000.0;

Stress_Ratio0 = STATEVn(:,19);
FHS_TIMEOLD = FHS_TIME;

P_reverse0 = STATEVn(:,14);
Q_reverse0 = STATEVn(:,15);


% 历史pl
TIME_REVERSE0 = STATEVn(:,13);
dpl_old = STATEVn(:,16);
dpl_old2 = STATEVn(:,17);
dpl_old3 = STATEVn(:,20);

% 参考模量,保证一周期内模量能完全循环(0.5-0.6用0.4-0.5模量）
Time_CK = rem(TIME+1e-6-Cal_par.Gravity_Time,Cal_par.Load_Time + Cal_par.Stop_Time);
modulus_CK = [];
% if Time_CK<0.01
%     STATEV(:,13:17) = 0.0;
% end

% modulus_dw = STATEV(:,13:17);
% T0 = Cal_par.Load_Time/20;
% if ((T0-Cal_par.dt)<1e-4) && (Time_CK+1e-3 >= Cal_par.Load_Time/2.0) && (sum(sum(modulus_dw)>0)==size(modulus_dw,2)) && (Time_CK<=Cal_par.Load_Time-1e-3)
%     modulus_t = [...
%         Cal_par.Load_Time/2.0+4*T0,...
%         Cal_par.Load_Time/2.0+3*T0,...
%         Cal_par.Load_Time/2.0+2*T0,...
%         Cal_par.Load_Time/2.0+T0,...
%         Cal_par.Load_Time/2.0];
% 
%     for i = 1:length(modulus_t)-1
%         if Time_CK>=modulus_t(i+1) && Time_CK<modulus_t(i)+0.001
%             BS_LS = (modulus_t(i) - Time_CK) / (modulus_t(i) - modulus_t(i+1));
%             modulus_CK = modulus_dw(:,i) + BS_LS*(modulus_dw(:,i+1)-modulus_dw(:,i));
%             break;
%         end
%     end 
% end

S_rA = S_r;
% UnSAT_sA = UnSAT_s;
SsaA = Ssa;
d_s = zeros(N,1);

% ==================   2.0   Original state ==================
STRAN_V = STRAN_OLD(:,1)+STRAN_OLD(:,2)+STRAN_OLD(:,3);                 % Volum strain
STRAN_VEND = STRAN_V+DSTRAN_OLD(:,1)+DSTRAN_OLD(:,2)+DSTRAN_OLD(:,3);
e0 = SUB_Gs./(SUB_K0.*SUB_pdmax)-1.0;                               % Calculate e0
e = (1.0+e0).*(1+STRAN_V)-1.0;                                     % Void ratio
e_END = (1.0+e0).*(1+STRAN_VEND)-1.0;                              % Void ratio
emin = SUB_Gs./SUB_pdmax - 1.0;                                    % emin (KK = 100%)
FHS_KK = (emin+1.0)./((e_END+e)/2+1.0);                            % Degree of compaction

% ==================   2.1   Calculate Parameter considering unsaturation ==================
% Calculate Saturation
Sr = Mc.*SUB_Gs./e./1000;
dSr = (Mc+dMc).*SUB_Gs./e_END/1000 - Sr;

% Calculate suction
if TIME == 0.0
    % ((mpData(mp).Sr^(-1/Fm) - 1)^(1/Fn))*Fa/mpData(mp).e^(1/(Fm*Fn))
    UnSAT_s = (10.^(log10(UnSAT_m3D).*0.5+log10(UnSAT_m3W)*0.5))./...
        e.^(1./UnSAT_m1./UnSAT_m2).*(Sr.^(-1./UnSAT_m1)-1.0).^(1./UnSAT_m2);
    S_r = UnSAT_s;
    Ssa = UnSAT_m3W./e.^(1./UnSAT_m1./UnSAT_m2).*...
        (Sr.^(-1./UnSAT_m1)-1.0).^(1./UnSAT_m2);

    S_rA = S_r;
    SsaA = Ssa;
end

% Calculate EPSL
F_S=3.0*UnSAT_Ts./UnSAT_R./UnSAT_s.*...
    (sqrt(9.+8.*UnSAT_R.*UnSAT_s./UnSAT_Ts)-3.).*...
    (sqrt(9.+8.*UnSAT_R.*UnSAT_s./UnSAT_Ts)+1.)/16.0;
EPSL = F_S .* (1 - Sr);

% Update NCL and CSL
UnSAT_LS = 1.0 - UnSAT_a .* (1.0-exp(UnSAT_b.*EPSL));
NCL_N      = UnSAT_LS .* NCL_N;
NCL_lapta  = UnSAT_LS .* NCL_lapta;
CSL_T      = UnSAT_LS .* CSL_T;
CSL_w      = UnSAT_LS .* CSL_w;

% Calculate suction
if TIME == 0.0
    UnSAT_LS_OLD = UnSAT_LS;
    UnLoad_k_OLD = UnLoad_k;
end

%      UnLoad_k = UnLoad_k*(1-0.35*(1-Sr))

% ==================   2.2   Calculate p, q ==================
% Calculate principal stress
PS = Calculate_principal(STRESS_OLD);
% Calculate P and Q
P1=  max(1000,-(PS(:,1)+PS(:,2)+PS(:,3))/3.0 + Sr .* UnSAT_s - Pore_Pressure);
P1_FHS = -(PS(:,1)+PS(:,2)+PS(:,3))/3.0 - Pore_Pressure;
Q1 = sqrt(((PS(:,1)-PS(:,2)).^2+(PS(:,2)-PS(:,3)).^2+(PS(:,1)-PS(:,3)).^2)/2.);

% ==================   2.3   Original state ==================
if plmax == 0.0
    plmax = PA.*exp((NCL_N-(1+e0)-UnLoad_k.*log(P1/PA+1.0))./...
        (NCL_lapta-UnLoad_k))-PA;                     % INITIAL plmax
end

% ==================   2.3   Calculate pl and Update plmax ==================
if Cal_par.Calculate_time <=Cal_par.Gravity_Time ...
        || (Time_CK+0.1*Cal_par.dt) < Cal_par.dt
    FHS_SR = UnSAT_s.*Sr;
    qr = Q1;
end

% if Time_CK>=0 && Time_CK < Cal_par.dt
%     P_reverse = P1;
%     Q_reverse = Q1;
% end

% Bounding_r = 2.0;
Fi = e-(CSL_T-CSL_w.*log(P1./PA+1.0)-1.0);
Md = M0.*exp(nd.*Fi);

% ==================   2.4   Calculate plastic modulus K, G ==================

% Calculate modulus
if isempty(modulus_CK)
    EMOD_OLD = FHS_k0.*PA.*(max(1000,P1_FHS+FHS_SR)./PA).^FHS_k1 ...
        .*abs(qr./PA+1.0).^FHS_k2 ...
        .*abs(max(Q1-qr,0.)./PA+1.).^FHS_k3.*FHS_KK.^FHS_k4.*1000.0;

    EMOD_OLD = min(EMOD_OLD,3.*FHS_k0.*PA*1000.0);
    EMOD_OLD = max(EMOD_OLD,0.01.*FHS_k0.*PA*1000.0);
else
    EMOD_OLD = modulus_CK;
end

% EMOD_OLD = 100000000*ones(size(EMOD_OLD));

Modulus_K = EMOD_OLD./3.0./(1.0-2.0.*SUB_u);
Modulus_G = EMOD_OLD./2.0./(1.0 + SUB_u);

EG2 = 2.0*Modulus_G;
ELAM = Modulus_K - 2.0/3.0*Modulus_G;

% INITIAL DDSDDE
DDSDDE = DDSDDE_Cal(ELAM,EG2);

% ==================   2.5   Calculate plastic modulus Kp and flow rule Ds ==================
P0=PA.*exp((NCL_N-(1.0+e)-UnLoad_k.*log(P1./PA+1.0))./...
    (NCL_lapta-UnLoad_k))-PA;

P0_UnSAT=PA.*exp((NCL_N./UnSAT_LS.*UnSAT_LS_OLD-(1.0+e)...
    -UnLoad_k_OLD.*log(P1./PA+1.0))./...
    (NCL_lapta./UnSAT_LS.*UnSAT_LS_OLD-UnLoad_k_OLD))-PA;

if TIME > 0.0
    plmax = STATEVn(:,1) .* 1000.0 .* P0 ./ P0_UnSAT;
end

% 返回弹性矩阵
if Cal_par.NRit == 0 && strcmp(Cal_par.solve_method, 'implicit')
    Dalg = DDSDDE;
    sig = sign;
    epsE = epsEn;
    epsP = zeros(N,6);
    Out1 = [];
    return
end
ERR2 = ones(N,1);
DDSDDE_old = DDSDDE;
%
%% ================== Iteration position ==================
% ****************** Calculate DDSDDE
while (FHS_JS < 2 || (any(ERR2 > 5e-4) && FHS_JS <= 25))

    % ****************** Calculate DDSDDE
    % Calculate dp, dq
    Fresh_YES = (ERR2 > 5e-4);

    dStran = permute(DSTRAN, [2, 3, 1]);        % 变为 6×1×N
    dStress = pagemtimes(DDSDDE, dStran);       % 结果 6×1×N
    STRESS0 = STRESS + permute(dStress, [3, 1, 2]);

    % Calculate principal stress
    PS2 = Calculate_principal(STRESS0);

    % Calculate P and Q
    P2=  max(1000,-(PS2(:,1)+PS2(:,2)+PS2(:,3))/3.0 + (Sr+dSr).*(UnSAT_s+d_s) - Pore_Pressure);
    P2_FHS = -(PS2(:,1)+PS2(:,2)+PS2(:,3))/3.0 - Pore_Pressure;
    Q2 = sqrt(((PS2(:,1)-PS2(:,2)).^2+(PS2(:,2)-PS2(:,3)).^2.+(PS2(:,1)-PS2(:,3)).^2)/2.0);

    % Calculate dP and dQ
    dP = P2_FHS-P1_FHS;
    dQ = Q2-Q1;
    dPP = dP;

    % P3 = P2 - 0.5*dP;
    Q3 = Q2 - 0.5*dQ;

    pl = P1.*exp(log(Bounding_r).*(Q1./M0./P1).^Bounding_n);
    pl = min(0.99*P0,pl);
    plmax = max(plmax,pl);
    plmax = min(plmax,P0);
    
    % pl 计算
    dpl = min((P2.*exp(log(Bounding_r).*(Q2./M0./P2).^Bounding_n)),0.99*P0) - ...
          min((P1.*exp(log(Bounding_r).*(Q1./M0./P1).^Bounding_n)),0.99*P0) ;
     
    % 1. FHS_TPJX 触发（逐点独立） 
    FHS_TPJX = ((Q3-qr) >= 1000) & (dpl > 0) & (STRESS(:,2) <= -1000);
    
    % 应力反转点判定 !!!
    Reverse_YES = (dpl.*dpl_old<0 & dpl.*dpl_old2<0 & dpl.*dpl_old3<0) | ...
        (Time_CK>=0 & (Time_CK+0.1*Cal_par.dt) < Cal_par.dt);

    P_reverse = P_reverse0;
    Q_reverse = Q_reverse0;
    TIME_REVERSE = TIME_REVERSE0;
    Reverse_YES = Reverse_YES & ((TIME - TIME_REVERSE)>=Cal_par.Load_Time*0.25);

    if sum(Reverse_YES)>0.0
        P_reverse(Reverse_YES) = P1(Reverse_YES);
        Q_reverse(Reverse_YES) = Q1(Reverse_YES);
        TIME_REVERSE(Reverse_YES) = TIME;
    end

    plR = P_reverse.*exp(log(Bounding_r).*(Q_reverse./M0./(1+P_reverse)).^Bounding_n);

    % Calculate stress ratio
    Stress_Ratio = Stress_Ratio0;
    stress0 =      min(3,Q1./max(P1_FHS+FHS_SR,1000));
    update_condition = (Stress_Ratio0 <= stress0);
    if sum(update_condition)>0
        if Cal_par.Calculate_time>(Cal_par.Gravity_Time+Cal_par.Load_Time+Cal_par.Stop_Time)
            Stress_Ratio(update_condition) = Stress_Ratio(update_condition)*0.999+...
                                             stress0(update_condition)*0.001;
            Stress_Ratio(~update_condition)= Stress_Ratio(~update_condition)*0.99+...
                                             stress0(~update_condition)*0.01;
        else
            Stress_Ratio(update_condition) = stress0(update_condition);
            Stress_Ratio(~update_condition)= Stress_Ratio(~update_condition)/2+...
                                             stress0(~update_condition)/2;
        end
    end
    % if Cal_par.Calculate_time <=Cal_par.Gravity_Time || Time_CK < Cal_par.dt
    %     Stress_Ratio = min(3,Q1./max(P1,1000));
    % end


    stress0 = max(stress0,0.01);
    Stress_Ratio = max(Stress_Ratio,0.01);

    % Calculate modulus
    if isempty(modulus_CK)
        EMOD_NEW = FHS_k0.*PA...
            .*abs(max(1000,(P1_FHS+P2_FHS)/2.0+FHS_SR)./PA).^FHS_k1...
            .*abs(qr./PA+1.0).^FHS_k2...
            .*abs(max(Q3-qr,0.)./PA+1.).^FHS_k3.*FHS_KK.^FHS_k4*1000.;

        EMOD_NEW = min(EMOD_NEW,3*FHS_k0*PA*1000.0);
        EMOD_NEW = max(EMOD_NEW,0.01*FHS_k0*PA*1000.0);
    else
        EMOD_NEW = modulus_CK;
        EMOD_OLD = modulus_CK;
    end
    
    EMOD_OLD =  EMOD_NEW*TEMP1 + EMOD_OLD*(1-TEMP1);
    ERR2 = abs(EMOD_NEW-EMOD_OLD)./EMOD_NEW*10.0;

    % EMOD_OLD = 100000000*ones(size(EMOD_OLD));

    Modulus_K = EMOD_OLD/3.0./(1.0-2.0*SUB_u);
    Modulus_G = EMOD_OLD/2.0./(1.0 + SUB_u);

    EG2 = 2.0*Modulus_G;
    ELAM = Modulus_K - 2.0/3.0*Modulus_G;

    % Calculate DDSDDE_E
    DDSDDE_E = DDSDDE_Cal(ELAM,EG2);
    Vs_YES = (ones(N,1)>2);
    %!!!!!!  ========================  ================================================ UPDATA DDSDDE
    %!!!!!!  ========================  LOOK LOOK there   ========================
    %!!!!!!  ========================  ================================================ UPDATA DDSDDE
    if TIME > 0.0
        % Calculate Kp
        nc_PAR1 = Fabric_PAR1+(plmax/PA).^(Fabric_PAR2.*plmax/PA);
        % nc_PAR1 = min(1e4,nc_PAR1);
        
        nc_PAR = nc_PAR1.*(log((P0./plmax).^(Fabric_PAR0.*nc_PAR1))+1.0);
        FHS_FF = (P0./plmax).^nc_PAR;
        FHS_FF = min(1e18,FHS_FF);
        
        Kp=(1.0+(e+e_END)/2)./(log(Bounding_r).*...
            (NCL_lapta-UnLoad_k)).*...
            (Md.^2.*max(1.0,min(100,(P0-plR)./(pl-plR))).^(Fabric_nc).*FHS_FF-...
            (stress0).^2)./(2.*stress0);
        
        FHS_TPJX(( ((P0-plR)<0) | ((pl-plR)<0) | Kp>1e15) & FHS_TPJX>0) = 0;

        if sum(~isreal(Kp),'all')> 0
            error('  塑性模量计算出现异常虚数，请检查！！！')
        end

        SsaD = UnSAT_m3D./e.^(1./UnSAT_m1./UnSAT_m2).*...
            (Sr.^(-1./UnSAT_m1)-1.0).^(1./UnSAT_m2);
        SsaW = UnSAT_m3W./e.^(1./UnSAT_m1./UnSAT_m2).*...
            (Sr.^(-1./UnSAT_m1)-1.0).^(1./UnSAT_m2);


        if FHS_JS > 1
            
            % Calculate Ssa
            Yes_LS = (d_s/DTIME > 10.0) & (ds_OLD/DTIME > 10.0);
            Ssa = SsaW;
            Ssa(Yes_LS) = SsaD(Yes_LS);

            % Update S_r
            Yes_LS = abs((Ssa-SsaA)./(UnSAT_m3W-UnSAT_m3D)) > 0.20;
            S_r = S_rA;
            FHS_TIME = FHS_TIMEOLD;

            S_r(Yes_LS) = UnSAT_s(Yes_LS);
            FHS_TIME(Yes_LS) = TIME;

            % 异常反复横跳
            Yes_LS = abs(FHS_TIMEOLD-FHS_TIME) < 0.05;
            FHS_TIME(Yes_LS) = FHS_TIMEOLD(Yes_LS);
            S_r(Yes_LS) = S_rA(Yes_LS);
            Yes_LS1 = Yes_LS & SsaA>S_r;
            Yes_LS2 = Yes_LS & SsaA<=S_r;
            Ssa(Yes_LS1) = SsaD(Yes_LS1);
            Ssa(Yes_LS2) = SsaW(Yes_LS2);
        end

        % Calculate Rs
        if UnSAT_m3W==UnSAT_m3D
            Rs = 1.0;
        else
            Rs=(log10(UnSAT_s/1000.)-log10(S_r/1000.))./...
                (log10(Ssa/1000.)-log10(S_r/1000.));
            Rs = max(0.01,Rs);
            Rs = min(1.000,Rs);
        end

        %% ！！！后续需要替换回来

        % !!!!!!  ======================== Calculate partial derivative /Fh
        dSsa_dSr = UnSAT_s./UnSAT_m1./UnSAT_m2./(Sr.^(1./UnSAT_m1)-1.)./Sr;
        dSsa_dsvp= UnSAT_s.*(1.0+e)./UnSAT_m1./UnSAT_m2./e;
        dSsa_dP  = UnSAT_s.*(1.0+e)./UnSAT_m1./UnSAT_m2./e./Modulus_K;
        dFh_ds   = ones(N,1);
        dFh_dSsa = -ones(N,1);

        % !!!!!!  ======================== Calculate partial derivative /Fs
        if any(FHS_TPJX)
            dFs_dp=-Bounding_n.*(Q1./M0./P1).^Bounding_n./P1+1./P1./log(Bounding_r);
            dFs_dq = Bounding_n.*(1.0./M0./P1).^Bounding_n.*Q1.^(Bounding_n-1.0);
            dFs_dP0 = -1./log(Bounding_r)./P0;
            % Ds = min((Md.^2.*max(0.00,min(1.00,pl./P0)).^(Fabric_nc./2)-(stress0).^2)./(2.*stress0),100.0);
            Ds = min((Md.^2-(stress0).^2)./(2.*stress0),100.0);
            dP0_dEPSL = (P0+PA).*UnSAT_a.*UnSAT_b.*exp(UnSAT_b.*EPSL).*...
                (NCL_N./UnSAT_LS./(NCL_lapta-UnLoad_k)-...
                (NCL_N-(1.0+e)-UnLoad_k.*log(P1./PA+1.)).*NCL_lapta...
                ./UnSAT_LS./(NCL_lapta-UnLoad_k).^2);
            dFs_ds = 3.*UnSAT_Ts/16./UnSAT_R.*(...
                (sqrt(9./UnSAT_s+8.*UnSAT_R./UnSAT_Ts)-3.*UnSAT_s.^(-0.5)).*...
                (-4.5.*(9./UnSAT_s+8.*UnSAT_R./UnSAT_Ts).^(-0.50).*UnSAT_s.^(-2)-...
                0.50.*UnSAT_s.^(-1.50))+...
                (sqrt(9./UnSAT_s+8.*UnSAT_R./UnSAT_Ts)+UnSAT_s.^(-0.5)).*...
                (-4.5.*(9./UnSAT_s+8.*UnSAT_R./UnSAT_Ts).^(-0.50).*UnSAT_s.^(-2)+...
                1.50.*UnSAT_s.^(-1.50)) );

            dEPSL_ds = (1-Sr).*dFs_ds;
            dEPSL_dSr = -F_S;

            dP0_dP = zeros(N,1);

            % Calculate FHS_A
            FHS_A = zeros(2,2,N);
            FHS_A(1,1,:) = Kp;
            FHS_A(1,2,:) = -(dFs_dP0.*dP0_dEPSL.*dEPSL_ds);
            FHS_A(2,1,:) = dFh_dSsa.*dSsa_dsvp.*Ds;
            FHS_A(2,2,:) = dFh_ds;

            % Calculate FHS_Y
            FHS_Y = zeros(2,1,N);
            FHS_Y(1,1,:) = (dFs_dp+dFs_dP0.*dP0_dP).*dP+dFs_dq.*dQ+...
                dFs_dP0.*dP0_dEPSL.*dEPSL_dSr.*dSr;
            FHS_Y(2,1,:) = -dFh_dSsa.*(dSsa_dP.*dP + dSsa_dSr./Rs.*dSr);

            % FHS_A: 2×2×N,  FHS_Y: 2×1×N
            FHS_X = pagemldivide(FHS_A, FHS_Y);  % 结果: 2×1×N
            % 转为 N×2
            FHS_X = reshape(FHS_X, 2, N)';       % 先变 2×N，再转置为 N×2

            Vs  = FHS_X(:,1);
            d_s = FHS_X(:,2);
            
            Vs_YES = (Vs.*dpl < 0.0);
            
            if sum(Vs_YES)>0.0
                d_s(Vs_YES) = -dFh_dSsa(Vs_YES)...
                    .*(dSsa_dP(Vs_YES).*dP(Vs_YES) + dSsa_dSr(Vs_YES)...
                    ./Rs.*dSr(Vs_YES));

                Vs(FHS_TPJX & Vs_YES) = 0;
                Ds(FHS_TPJX & Vs_YES) = 0;
                Kp(FHS_TPJX & Vs_YES) = 0;
            end

            % Caculate FHS_NB
            FHS_NB = dFs_dP0.*dP0_dEPSL.*(dEPSL_dSr.*dSr+dEPSL_ds.*d_s);
            FHS_NB = ((dFs_dp+dFs_dP0.*dP0_dP).*dP+dFs_dq.*dQ + FHS_NB)./...
                ((dFs_dp+dFs_dP0.*dP0_dP).*dPP+dFs_dq.*dQ);
            FHS_NB = max(0.50,FHS_NB);
            FHS_NB = min(2.00,FHS_NB);

            % !!!!!!  ========================  ================================================ UPDATA DDSDDE
            % !!!!!!  ========================  LOOK LOOK there   ========================
            % !!!!!!  ========================  ================================================ UPDATA DDSDDE

            % Calculate Hij
            valid = FHS_TPJX & ~Vs_YES;
            for I=1:3
                HH(valid,I) = FHS_NB(valid)./Kp(valid)...
                    .*(3.0*STRESS_S(valid,I)/2.0./Q1(valid)+Ds(valid)/3.0);
            end
            for I=4:6
                HH(valid,I) = FHS_NB(valid)./Kp(valid)...
                    .*(3.0*STRESS_S(valid,I)/2.0./Q1(valid));
            end

            % Calculate FF1
            for I=1:6
                for J=1:3
                    FF1(I,J,valid)=HH(valid,I).*(1./3.*...
                        (dFs_dp(valid)+dFs_dP0(valid).*dP0_dP(valid)) + ...
                        3./2./Q1(valid).*dFs_dq(valid).*STRESS_S(valid,J));
                end
                for J=4:6
                    FF1(I,J,valid)=HH(valid,I).*(3.0./Q1(valid)...
                        .*dFs_dq(valid).*STRESS_S(valid,J));
                end

                if I>3
                    for J=1:6
                        FF1(I,J,valid) = FF1(I,J,valid)*2.0;
                    end
                end
            end

            % Multiply of Matrix
            DDSDDE_DE =  ONE(:,:,valid) + pagemtimes(DDSDDE_E(:,:,valid),FF1(:,:,valid));
            DDSDDE_NEW(:,:,valid) = pagemldivide(DDSDDE_DE, DDSDDE_E(:,:,valid));
        end
        
        d_s(~FHS_TPJX) = -dFh_dSsa(~FHS_TPJX)...
            .*(dSsa_dP(~FHS_TPJX).*dP(~FHS_TPJX) + ...
            dSsa_dSr(~FHS_TPJX)./Rs.*dSr(~FHS_TPJX));

        Vs(~FHS_TPJX) = 0;
        Ds(~FHS_TPJX) = 0;
        %Kp(~FHS_TPJX) = 0;

        % !!!!!! UPDATA DDSDDE
        if (FHS_JS > 1)
            TEMP_LS = min(TEMP_LS*1.20,0.80);
        end

        DDSDDE(:,:,Fresh_YES) = DDSDDE(:,:,Fresh_YES) * (1-TEMP_LS) + DDSDDE_NEW(:,:,Fresh_YES) * TEMP_LS;
        % DDSDDE(:,:,Fresh_YES) = max(DDSDDE(:,:,Fresh_YES), 0.01*DDSDDE_E(:,:,Fresh_YES));
        
        % Vs<0的点：直接设为弹性矩阵（等价于原版continue） 
        % ~FHS_TPJX的点：直接设为弹性矩阵
        DDSDDE(:,:,Fresh_YES & (Vs_YES | ~FHS_TPJX)) = DDSDDE_E(:,:,Fresh_YES & (Vs_YES | ~FHS_TPJX));
        
        % !!!!!! Calculate Err
        error_mat = abs((DDSDDE-DDSDDE_old)./DDSDDE);
        error_mat(isnan(error_mat)) = 0;
        ERR2 = squeeze(sum(error_mat, [1 2]));  % 结果 N×1
        DDSDDE_old = DDSDDE;
        % if strcmp(Cal_par.solve_method, 'explicit')
        %     Vs(Fresh_YES & (Vs_YES | ~FHS_TPJX)) = 0;
        %     Ds(Fresh_YES & (Vs_YES | ~FHS_TPJX)) = 0;
        %     disp(['       == >>> 总数 = ',num2str(sum(Fresh_YES)),' Vs_异常 = ',num2str(sum(Vs_YES)),' 塑性判定 = ',...
        %         num2str(sum(FHS_TPJX)), ' 塑性 = ',num2str(sum(Fresh_YES) - sum(Fresh_YES & (Vs_YES | ~FHS_TPJX)))])
        % end

    else
        Vs = zeros(N,1);
        Ds = zeros(N,1);
        %Kp = zeros(N,1);
        % !!!!!! UPDATA DDSDDE
        DDSDDE = DDSDDE_E;
    end

    FHS_JS = FHS_JS + 1;

end

% Update STRESS
DSTRESS_3d = pagemtimes(DDSDDE,dStran);
STRESS = STRESS + permute(DSTRESS_3d, [3, 1, 2]);

% Update Suction
UnSAT_s = max(1000,UnSAT_s + d_s);

% ================== SAVE STATEV ==================
STATEV(:,1) = plmax/1000.0;                         % Void ratio: e
STATEV(:,2) = qr/1000.0;                            % Initial void ratio: Vs
STATEV(:,3) = S_r/1000.0;
STATEV(:,4) = UnSAT_LS;
STATEV(:,5) = UnLoad_k;
STATEV(:,6) = d_s;
STATEV(:,7) = UnSAT_s/1000.0;
STATEV(:,8) = Ssa/1000.0;
if Cal_par.Calculate_time>1
    STATEV(:,9) = Vs_all1 + Vs;
    STATEV(:,10)= Vs_all2 + Vs.*Ds;
else
    STATEV(:,9) = 0.0;
    STATEV(:,10)= 0.0;
end
STATEV(:,11)= FHS_TIME;
STATEV(:,12)= FHS_SR/1000.0;


STATEV(:,14)= P_reverse;
STATEV(:,15)= Q_reverse;

STATEV(:,13) = TIME_REVERSE;
STATEV(:,16)= dpl;
STATEV(:,17)= dpl_old;
STATEV(:,20)= dpl_old2;


%% 储存应力状态，确保单周期能够回到原点
% if Cal_par.dt>T0*0.95
%     if abs(Time_CK) < 0.001               % 0.0~0.1模量
%         STATEV(:,13)= EMOD_OLD;
%     elseif abs(Time_CK-T0) < 0.001        % 0.1~0.2模量
%         STATEV(:,14)= EMOD_OLD;
%     elseif abs(Time_CK-2*T0) < 0.001      % 0.2~0.3模量
%         STATEV(:,15)= EMOD_OLD;
%     elseif abs(Time_CK-3*T0) < 0.001      % 0.3~0.4模量
%         STATEV(:,16)= EMOD_OLD;
%     elseif abs(Time_CK-4*T0) < 0.001      % 0.4~0.5模量
%         STATEV(:,17)= EMOD_OLD;
%     end
% end
if Time_CK < 1e-6
    STATEV(:,18)= STRESS(:,2);
end

STATEV(:,19)= Stress_Ratio;

% mpData.Scution_Flow = STATEV(7);
% mpData.modulus = STATEV(18);
% mpData.STATEV = STATEV;

Out1{1} = EMOD_OLD/1.0E6;
Out1{2} = STATEV;

%% 返回主程序
Dalg = DDSDDE;
sig = STRESS;

depsE = pagemldivide(DDSDDE_E, DSTRESS_3d);
epsE = epsEn + permute(depsE, [3, 1, 2]);
epsP = eps - epsE;

% if sum(Reverse_YES)>0.0
% disp(['  反转吧：',num2str(TIME)]);
% end
end



%% 批量主应力计算函数
function PS = Calculate_principal(STRESS)
% 输入: STRESS (N×6) [S11, S22, S33, S12, S13, S23]
% 输出: PS (N×3) 排序后的主应力 [大, 中, 小]
N = size(STRESS,1);
if N == 0
    PS = zeros(0,3);
    return;
end

% 构造3×3×N张量
S = zeros(3,3,N);
S(1,1,:) = STRESS(:,1);
S(2,2,:) = STRESS(:,2);
S(3,3,:) = STRESS(:,3);
S(1,2,:) = STRESS(:,4); S(2,1,:) = STRESS(:,4);
S(1,3,:) = STRESS(:,5); S(3,1,:) = STRESS(:,5);
S(2,3,:) = STRESS(:,6); S(3,2,:) = STRESS(:,6);

% 批量特征值分解
[~, D] = pageeig(S);  % D是3×3×N对角阵
eig_vals = [squeeze(D(1,1,:)), squeeze(D(2,2,:)), squeeze(D(3,3,:))];

% 排序 (每行降序)
PS = sort(eig_vals, 2, 'descend');
end

function DDSDDE = DDSDDE_Cal(ELAM,EG2)
N = size(ELAM,1);
DDSDDE = zeros(6,6,N);
for I=1:3
    for J=1:3
        DDSDDE(I,J,:)=ELAM;
    end
    DDSDDE(I,I,:)=EG2+ELAM;
end
for I=4:6
    DDSDDE(I,I,:)=EG2/2.0;
end
end