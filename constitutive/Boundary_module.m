function [Dalg,sig,epsE,epsP,mpData] = Boundary_module(eps,mpData,Cal_par)

%% 设定计算常数
PA = 101325;
FHS_JS = 1;
TEMP1 = 0.90;
TEMP2 = 0.50; TEMP_LS = TEMP2;
FHS_TPJX = 0.00;
nc_PAR = 0.0;
DTIME = Cal_par.dt;
% ================== Initial parameters
%       SET ZERO
TIME = Cal_par.Calculate_time - Cal_par.dt;
STRESS_OLD = mpData.sign;
STRAN_OLD = mpData.epsEn + mpData.epsPn;
STRESS = mpData.sign;
STRESS0 = mpData.sign;
DSTRAN_OLD = eps - mpData.epsEn - mpData.epsPn;
DSTRAN = eps - mpData.epsEn - mpData.epsPn;
HH = zeros(6,1);

FF1 = zeros(6,6);
ONE = eye(6);
FF2 = zeros(6,1);

% CALCULATE sij
SS1 = (STRESS0(1)+STRESS0(2)+STRESS0(3))/3.0;
STRESS_S(1:3) = -(STRESS0(1:3) - SS1);  % 或写成：-STRESS0(1:3) + SS1
STRESS_S(4:6) = -STRESS0(4:6);  % 修正索引！

PROPS = mpData.mCst;
% ================== READ PROPS PRAMTS ==================
NCL_N      = PROPS(1);	    % Critical state soil mechanics
NCL_lapta  = PROPS(2);	    % Critical state soil mechanics
CSL_T      = PROPS(3);	    % Critical state soil mechanics
CSL_w      = PROPS(4);	    % Critical state soil mechanics
UnLoad_k   = PROPS(5);	    % Critical state soil mechanics
nd         = PROPS(6);      % Critical state soil mechanics
M0         = PROPS(7);      % Critical state soil mechanics
Bounding_n = PROPS(8);	    % Parameters for bounding surface
Fabric_nc  = PROPS(9);	    % Hardening parameter
Fabric_PAR0= PROPS(10);	    % Hardening parameter
Fabric_PAR1= PROPS(11);	    % Hardening parameter
Fabric_PAR2= PROPS(12);	    % Hardening parameter
Bounding_r = PROPS(13);     % Boundary shape parameter
SUB_pdmax  = PROPS(14);	    % Subgrade soil parameters
SUB_Gs     = PROPS(15);	    % Subgrade soil parameters
SUB_K0     = PROPS(16);	    % Subgrade soil parameters
SUB_u      = PROPS(17);    	% Subgrade soil parameters
FHS_k0     = PROPS(18);	    % Elastic Parameters
FHS_k1     = PROPS(19);	    % Elastic Parameters
FHS_k2     = PROPS(20);	    % Elastic Parameters
FHS_k3     = PROPS(21);	    % Elastic Parameters
FHS_k4     = PROPS(22);	    % Elastic Parameters

% ------- Unsaturation parameters - (NEW PROPS)
UnSAT_a    = PROPS(23);       % state dependent parameter of NCL and CSL
UnSAT_b    = PROPS(24);       % state dependent parameter of NCL and CSL
UnSAT_m1   = PROPS(25);       % SWCC Parameters
UnSAT_m2   = PROPS(26);       % SWCC Parameters
UnSAT_m3D  = PROPS(27);       % SWCC Parameters
UnSAT_m3W  = PROPS(28);       % SWCC Parameters
UnSAT_R    = PROPS(29);       % f(s) Parameters

% ------- Unsaturation parameters - (NEW PROPS)
% Moisture_a = PROPS(30);       % moisture parameter
% Moisture_b = PROPS(31);       % moisture parameter
% Moisture_N = PROPS(32);       % moisture parameter

% ------- Unsaturation parameters - (NEW PROPS)
UnSAT_Ts = 72.8/1000.0;

% ****************** READ STATEV
STATEVn = mpData.STATEVn;
plmax = STATEVn(1)*1000.0;   % Memory stress: plmax
qr = STATEVn(2)*1000.0;      % Reference Stress: qr
S_r = STATEVn(3)*1000.0;     % Reference suction: sr
UnSAT_LS_OLD = STATEVn(4);
UnLoad_k_OLD = STATEVn(5);
ds_OLD = STATEVn(6);
UnSAT_s = STATEVn(7)*1000.0;
Ssa = STATEVn(8)*1000.0;
Vs_all1 = STATEVn(9);
Vs_all2 = STATEVn(10);
FHS_TIME = STATEVn(11);
FHS_SR = STATEVn(12)*1000.0;
FHS_TIMEOLD = FHS_TIME;

S_rA = S_r;
% UnSAT_sA = UnSAT_s;
SsaA = Ssa;
d_s = 0.0;

% ****************** Calculate Mc dMc
Mc = mpData.Gravity_Water_Content_old;
dMc = mpData.Gravity_Water_Content - mpData.Gravity_Water_Content_old;

% ==================   2.0   Original state ==================
STRAN_V = STRAN_OLD(1)+STRAN_OLD(2)+STRAN_OLD(3);                 % Volum strain
STRAN_VEND = STRAN_V+DSTRAN_OLD(1)+DSTRAN_OLD(2)+DSTRAN_OLD(3);
e0 = SUB_Gs/(SUB_K0*SUB_pdmax)-1.0;                               % Calculate e0
e = (1.0+e0)*(1+STRAN_V)-1.0;                                     % Void ratio
e_END = (1.0+e0)*(1+STRAN_VEND)-1.0;                              % Void ratio
% e_P0 =  (1.0+e0)*(1-Vs_all2)-1.0;
emin = SUB_Gs/SUB_pdmax - 1.0;                                    % emin (KK = 100%)
FHS_KK = (emin+1.0)/((e_END+e)/2+1.0);                            % Degree of compaction
% Calculate body load
mpData.mpM = FHS_KK*SUB_pdmax*(1+Mc)*mpData.vp;

% ==================   2.1   Calculate Parameter considering unsaturation ==================
% Calculate Saturation
Sr = Mc*SUB_Gs/e/1000;
dSr = (Mc+dMc)*SUB_Gs/e_END/1000 - Sr;

% Calculate suction
if TIME == 0.0
    % ((mpData(mp).Sr^(-1/Fm) - 1)^(1/Fn))*Fa/mpData(mp).e^(1/(Fm*Fn))
    UnSAT_s = (10^(log10(UnSAT_m3D)*0.5+log10(UnSAT_m3W)*0.5))/...
        e^(1./UnSAT_m1/UnSAT_m2)*(Sr^(-1./UnSAT_m1)-1.0)^(1./UnSAT_m2);
    S_r = UnSAT_s;
    Ssa = UnSAT_m3W/e^(1./UnSAT_m1/UnSAT_m2)*...
        (Sr^(-1./UnSAT_m1)-1.0)^(1./UnSAT_m2);

    S_rA = S_r;
    % UnSAT_sA = UnSAT_s;
    SsaA = Ssa;
end

% Calculate EPSL
F_S=3.0*UnSAT_Ts/UnSAT_R/UnSAT_s*...
    (sqrt(9.+8.*UnSAT_R*UnSAT_s/UnSAT_Ts)-3.)*...
    (sqrt(9.+8.*UnSAT_R*UnSAT_s/UnSAT_Ts)+1.)/16.0;
EPSL = F_S * (1 - Sr);

% Update NCL and CSL
UnSAT_LS = 1.0 - UnSAT_a * (1.0-exp(UnSAT_b*EPSL));
NCL_N      = UnSAT_LS * NCL_N;
NCL_lapta  = UnSAT_LS * NCL_lapta;
CSL_T      = UnSAT_LS * CSL_T;
CSL_w      = UnSAT_LS * CSL_w;

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
P=  -(PS(1)+PS(2)+PS(3))/3.0 + Sr * UnSAT_s;
P1_FHS = -(PS(1)+PS(2)+PS(3))/3.0;
Q = sqrt(((PS(1)-PS(2))^2+(PS(2)-PS(3))^2+(PS(1)-PS(3))^2)/2.);

% ==================   2.3   Original state ==================
if plmax == 0.0
    plmax = PA*exp((NCL_N-(1+e0)-UnLoad_k*log(P/PA+1.0))/...
        (NCL_lapta-UnLoad_k))-PA;                     % INITIAL plmax
end

% ==================   2.3   Calculate pl and Update plmax ==================
% Bounding_r = 2.0;
Fi = e-(CSL_T-CSL_w*log(P/PA+1.0)-1.0);
Md = M0*exp(nd*Fi);
if P == 0.0
    pl = 0.0;
else
    pl = P*exp(log(Bounding_r)*(Q/M0/P)^Bounding_n);
end

% ==================   2.4   Calculate plastic modulus K, G ==================
if (abs(TIME-round(TIME)) < 0.01) || (TIME<=1.01)
    FHS_SR = UnSAT_s*Sr;
    qr = Q;
    % PR = P;           % Initial time in every analysis step
end

% Calculate modulus
EMOD_OLD = FHS_k0*PA*abs((max(0.0,P1_FHS)+FHS_SR)/PA)^FHS_k1 ...
    *abs(qr/PA+1.0)^FHS_k2 ...
    *abs(max(Q-qr,0.)/PA+1.)^FHS_k3*FHS_KK^FHS_k4*1000.0;

EMOD_OLD = min(EMOD_OLD,3*FHS_k0*PA*1000.0);
EMOD_OLD = max(EMOD_OLD,0.33*FHS_k0*PA*1000.0);

Modulus_K = EMOD_OLD/3.0/(1.0-2.0*SUB_u);
Modulus_G = EMOD_OLD/2.0/(1.0 + SUB_u);

EG2 = 2.0*Modulus_G;
ELAM = Modulus_K - 2.0/3.0*Modulus_G;

% INITIAL DDSDDE
DDSDDE = DDSDDE_Cal(ELAM,EG2);

% ==================   2.5   Calculate plastic modulus Kp and flow rule Ds ==================
P0=PA*exp((NCL_N-(1.0+e)-UnLoad_k*log(P/PA+1.0))/...
    (NCL_lapta-UnLoad_k))-PA;

P0_UnSAT=PA*exp((NCL_N/UnSAT_LS*UnSAT_LS_OLD-(1.0+e)...
    -UnLoad_k_OLD*log(P/PA+1.0))/...
    (NCL_lapta/UnSAT_LS*UnSAT_LS_OLD-UnLoad_k_OLD))-PA;

if TIME > 0.0
    plmax = STATEVn(1) * 1000.0 * P0 / P0_UnSAT;
end

pl = min(0.99*P0,pl);
plmax = max(plmax,pl);
plmax = min(plmax,P0);

% 返回弹性矩阵
if Cal_par.NRit == 0
    Dalg = DDSDDE;
    sig = mpData.sign;
    epsE = mpData.epsEn;
    epsP = mpData.epsPn;
    return
end

%      IF((TIME(1).GT.0.0))THEN
%      OPEN(UNIT=10,
%     1      file='E:\ABAQUS Temp\20250617 memory surface\outputs.txt',
%     2      STATUS = 'NEW')
%      END IF
%% ================== Iteration position ==================
% ****************** Calculate DDSDDE
while (FHS_JS < 2 || (ERR2 > 0.005 && FHS_JS <= 15))
    % ****************** Calculate DDSDDE
    % Calculate dp, dq
    STRESS0 = STRESS + DDSDDE*DSTRAN;	% INITIAL STRESS STATE !
    % Calculate principal stress
    PS2 = Calculate_principal(STRESS0);
    % Calculate P and Q
    P2=  -(PS2(1)+PS2(2)+PS2(3))/3.0 + (Sr+dSr)*(UnSAT_s+d_s);
    P2_FHS = -(PS2(1)+PS2(2)+PS2(3))/3.0;
    Q2 = sqrt(((PS2(1)-PS2(2))^2+(PS2(2)-PS2(3))^2.+(PS2(1)-PS2(3))^2)/2.0);

    dPP = P2_FHS - P1_FHS;

    % Calculate dP and dQ
    dP = P2-P;
    dQ = Q2-Q;

    P3 = P2 - 0.5*dP;
    Q3 = Q2 - 0.5*dQ;

    if FHS_TPJX == 0.0 && (Q3-qr) > 500 && dQ > 1.0 && (abs(TIME-round(TIME)) > 0.01) && STRESS(2) < -1000
        FHS_TPJX = 1.0;
    end

    % Calculate stress ratio
    Stress_Ratio=max(Q3/(max(P2_FHS+P1_FHS,0.0)/2.0+FHS_SR),0.10);


    % Calculate modulus
    EMOD_NEW = FHS_k0*PA...
        *abs((max(P2_FHS+P1_FHS,0.0)/2.0+FHS_SR)/PA)^FHS_k1...
        *abs(qr/PA+1.0)^FHS_k2...
        *abs(max(Q3-qr,0.)/PA+1.)^FHS_k3*FHS_KK^FHS_k4*1000.;

    EMOD_NEW = min(EMOD_NEW,3*FHS_k0*PA*1000.0);
    EMOD_NEW = max(EMOD_NEW,0.33*FHS_k0*PA*1000.0);

    ERR2 = abs((EMOD_NEW-EMOD_OLD)/EMOD_NEW)*100.0;

    EMOD_OLD =  EMOD_NEW*TEMP1 + EMOD_OLD*(1-TEMP1);
    Modulus_K = EMOD_OLD/3.0/(1.0-2.0*SUB_u);
    Modulus_G = EMOD_OLD/2.0/(1.0 + SUB_u);

    EG2 = 2.0*Modulus_G;
    ELAM = Modulus_K - 2.0/3.0*Modulus_G;

    %!!!!!!  ========================  ================================================ UPDATA DDSDDE
    %!!!!!!  ========================  LOOK LOOK there   ========================
    %!!!!!!  ========================  ================================================ UPDATA DDSDDE
    if TIME > 0.0
        % Calculate Kp
        nc_PAR1 = Fabric_PAR1+(plmax/PA)^(Fabric_PAR2*plmax/PA);
        nc_PAR1 = nc_PAR1/Stress_Ratio;

        nc_PAR = nc_PAR1*(log((P0/plmax)^(Fabric_PAR0*nc_PAR1))+1.0);
        FHS_FF = (P0/plmax)^nc_PAR;

        Kp=Fabric_nc*(1.0+(e+e_END)/2)/(log(Bounding_r)*...
            (NCL_lapta-UnLoad_k))*...
            (Md^2*(P0/pl)^2*FHS_FF-...
            (Q3/P3)^2)/(2*Q3/P3);

        SsaD = UnSAT_m3D/e^(1./UnSAT_m1/UnSAT_m2)*...
            (Sr^(-1./UnSAT_m1)-1.0)^(1./UnSAT_m2);
        SsaW = UnSAT_m3W/e^(1./UnSAT_m1/UnSAT_m2)*...
            (Sr^(-1./UnSAT_m1)-1.0)^(1./UnSAT_m2);

        % Calculate Ssa
        if FHS_JS > 1
            if (d_s/DTIME > 10.0) && (ds_OLD/DTIME > 10.0)
                Ssa = SsaD;
                % OFF = 1.0;
            else
                Ssa = SsaW;
                % OFF = 0.0;
            end

            % Update S_r
            if abs((Ssa-SsaA)/(UnSAT_m3W-UnSAT_m3D)) > 0.20
                S_r = UnSAT_s;
                FHS_TIME = TIME;
            else
                S_r = S_rA;
                FHS_TIME = FHS_TIMEOLD;
            end

            if abs(FHS_TIMEOLD-FHS_TIME) < 0.05
                FHS_TIME = FHS_TIMEOLD;
                S_r = S_rA;
                if(SsaA>S_r)
                    Ssa = SsaD;
                else
                    Ssa = SsaW;
                end
            end
        end

        % Calculate Rs
        % Rs=(log10(UnSAT_s/1000.)-log10(S_r/1000.))/...
        %     (log10(Ssa/1000.)-log10(S_r/1000.));
        % Rs = max(0.01,Rs);
        % Rs = min(1.000,Rs);

        %% ！！！
        Rs = 1.0;

        % !!!!!!  ======================== Calculate partial derivative /Fh
        dSsa_dSr = UnSAT_s/UnSAT_m1/UnSAT_m2/(Sr^(1./UnSAT_m1)-1.)/Sr;
        dSsa_dsvp= UnSAT_s*(1.0+e)/UnSAT_m1/UnSAT_m2/e;
        dSsa_dP  = UnSAT_s*(1.0+e)/UnSAT_m1/UnSAT_m2/e/Modulus_K;
        dFh_ds   = 1.0;
        dFh_dSsa = -1.0;

        % !!!!!!  ======================== Calculate partial derivative /Fs
        if FHS_TPJX>0.1
            dFs_dp=-Bounding_n*(Q/M0/P)^Bounding_n/P+1./P/log(Bounding_r);
            dFs_dq = Bounding_n*(1.0/M0/P)^Bounding_n*Q^(Bounding_n-1.0);
            dFs_dP0 = -1./log(Bounding_r)/P0;
            Ds = (Md^2-(Q3/P3)^2)/(2*Q3/P3);
            dP0_dEPSL = (P0+PA)*UnSAT_a*UnSAT_b*exp(UnSAT_b*EPSL)*...
                (NCL_N/UnSAT_LS/(NCL_lapta-UnLoad_k)-...
                (NCL_N-(1.0+e)-UnLoad_k*log(P/PA+1.))*NCL_lapta...
                /UnSAT_LS/(NCL_lapta-UnLoad_k)^2);
            dFs_ds = 3.*UnSAT_Ts/16./UnSAT_R*(...
                (sqrt(9./UnSAT_s+8.*UnSAT_R/UnSAT_Ts)-3.*UnSAT_s^(-0.5))*...
                (-4.5*(9./UnSAT_s+8.*UnSAT_R/UnSAT_Ts)^(-0.50)*UnSAT_s^(-2)-...
                0.50*UnSAT_s^(-1.50))+...
                (sqrt(9./UnSAT_s+8.*UnSAT_R/UnSAT_Ts)+UnSAT_s^(-0.5))*...
                (-4.5*(9./UnSAT_s+8.*UnSAT_R/UnSAT_Ts)^(-0.50)*UnSAT_s^(-2)+...
                1.50*UnSAT_s^(-1.50)) );

            dEPSL_ds = (1-Sr)*dFs_ds;
            dEPSL_dSr = -F_S;

            dP0_dP = 0.0;

            % Calculate FHS_A
            FHS_A(1,1) = Kp;
            FHS_A(1,2) = -(dFs_dP0*dP0_dEPSL*dEPSL_ds);
            FHS_A(2,1) = dFh_dSsa*dSsa_dsvp*Ds;
            FHS_A(2,2) = dFh_ds;

            % Calculate FHS_Y
            FHS_Y(1,1) = (dFs_dp+dFs_dP0*dP0_dP)*dP+dFs_dq*dQ+...
                dFs_dP0*dP0_dEPSL*dEPSL_dSr*dSr;
            FHS_Y(2,1) = -dFh_dSsa*(dSsa_dP*dP + dSsa_dSr/Rs*dSr);

            FHS_X = FHS_A\FHS_Y;
            Vs  = FHS_X(1,1);

            if Vs < 0.0 || (nc_PAR >= (log(EMOD_NEW*1000.0) / log(P0/plmax)))
                d_s = -dFh_dSsa*(dSsa_dP*dP + dSsa_dSr/Rs*dSr);
                Vs = 0.0;
                Ds = 0.0;
                Kp = 0.0;

                if ((UnSAT_s + d_s) <= SsaW)
                    d_s = (SsaW - UnSAT_s)*0.90;
                end

                if ((UnSAT_s + d_s) >= SsaD)
                    d_s = (SsaD - UnSAT_s)*0.90;
                end

                % !!!!!! UPDATA DDSDDE
                DDSDDE = DDSDDE_Cal(ELAM,EG2);
                FHS_JS = FHS_JS + 1;
                continue
            end

            d_s = FHS_X(2,1);

            if ((UnSAT_s + d_s) <= SsaW)
                d_s = (SsaW - UnSAT_s)*0.90;
            end

            if ((UnSAT_s + d_s) >= SsaD)
                d_s = (SsaD - UnSAT_s)*0.90;
            end

            % Caculate FHS_NB
            FHS_NB = dFs_dP0*dP0_dEPSL*(dEPSL_dSr*dSr+dEPSL_ds*d_s);
            FHS_NB = ((dFs_dp+dFs_dP0*dP0_dP)*dP+dFs_dq*dQ + FHS_NB)/...
                ((dFs_dp+dFs_dP0*dP0_dP)*dPP+dFs_dq*dQ);
            FHS_NB = max(0.50,FHS_NB);
            FHS_NB = min(2.00,FHS_NB);

            % !!!!!!  ========================  ================================================ UPDATA DDSDDE
            % !!!!!!  ========================  LOOK LOOK there   ========================
            % !!!!!!  ========================  ================================================ UPDATA DDSDDE

            % Calculate Hij
            for I=1:3
                HH(I) = FHS_NB/Kp*(3.0*STRESS_S(I)/2.0/Q+Ds/3.0);
            end
            for I=4:6
                HH(I) = FHS_NB/Kp*(3.0*STRESS_S(I)/2.0/Q);
            end

            % Calculate FF1
            for I=1:6
                for J=1:3
                    FF1(I,J)=HH(I)*(1./3.*(dFs_dp+dFs_dP0*dP0_dP) + ...
                        3./2./Q*dFs_dq*STRESS_S(J));
                end
                for J=4:6
                    FF1(I,J)=HH(I)*(3.0/Q*dFs_dq*STRESS_S(J));
                end

                if I>3
                    for J=1:6
                        FF1(I,J) = FF1(I,J)*2.0;
                    end
                end
            end

            % Calculate FF2
            for I=1:6
                FF2(I,1) = FHS_NB*HH(I);
                if I>3
                    FF2(I,1) = FF2(I,1)*2.0;
                end
            end

            % Calculate DDSDDE_E
            DDSDDE_E = DDSDDE_Cal(ELAM,EG2);

            % Multiply of Matrix
            DDSDDE_DE =  ONE + DDSDDE_E*FF1;
            DDSDDE_NEW = DDSDDE_DE\DDSDDE_E;

            % !!!!!! Calculate Err
            ERR2=abs((DDSDDE-DDSDDE_NEW)./DDSDDE_NEW);
            ERR2(isnan(ERR2)) = 0;
            ERR2=sum(ERR2(:));

            % !!!!!! UPDATA DDSDDE
            if (FHS_JS > 1)
                TEMP_LS = min(TEMP_LS*1.20,0.80);
            end

            DDSDDE = DDSDDE * (1-TEMP_LS) + DDSDDE_NEW * TEMP_LS;
            DDSDDE = max(DDSDDE, 0.01*DDSDDE_E);

        else

            d_s = -dFh_dSsa*(dSsa_dP*dP + dSsa_dSr/Rs*dSr);
            Vs = 0.0;
            Ds = 0.0;
            Kp = 0.0;

            if (UnSAT_s + d_s) <= SsaW
                d_s = (SsaW - UnSAT_s)*0.90;
            end

            if (UnSAT_s + d_s) >= SsaD
                d_s = (SsaD - UnSAT_s)*0.90;
            end

            %!!!!!! UPDATA DDSDDE
            DDSDDE = DDSDDE_Cal(ELAM,EG2);
        end

    else
        Vs = 0.0;
        Ds = 0.0;
        Kp = 0.0;
        % !!!!!! UPDATA DDSDDE
        DDSDDE = DDSDDE_Cal(ELAM,EG2);
    end
    FHS_JS = FHS_JS + 1;
end

DDSDDE_E = DDSDDE_Cal(ELAM,EG2);

% Update STRESS
DSTRESS = DDSDDE*DSTRAN;
STRESS=STRESS+DSTRESS;

% Update Suction
UnSAT_s = UnSAT_s + d_s;

% ================== SAVE STATEV ==================
STATEV(1) = plmax/1000.0;                         % Void ratio: e
STATEV(2) = qr/1000.0;                            % Initial void ratio: Vs
STATEV(3) = S_r/1000.0;
STATEV(4) = UnSAT_LS;
STATEV(5) = UnLoad_k;
STATEV(6) = d_s;
STATEV(7) = UnSAT_s/1000.0;
STATEV(8) = Ssa/1000.0;
STATEV(9) = Vs_all1 + Vs;
STATEV(10)= Vs_all2 + Vs*Ds;
STATEV(11)= FHS_TIME;
STATEV(12)= FHS_SR/1000.0;
STATEV(13)= Sr*100;
STATEV(14)= nc_PAR;
STATEV(15)= P0/1.0E3;
STATEV(16)= Kp/1.0E6;
STATEV(17)= Md;
STATEV(18)= EMOD_NEW/1.0E6;
STATEV(19)= Stress_Ratio;

mpData.Scution_Flow = STATEV(7);
mpData.modulus = STATEV(18);
mpData.STATEV = STATEV;

%% 返回主程序
Dalg = DDSDDE;
sig = STRESS;
epsE = mpData.epsEn + DDSDDE_E\DSTRESS;
epsP = eps - epsE;

% 渗流参数
Name_mp = fieldnames(mpData);
if sum(strcmp(Name_mp,'Flow_Type'))>0
    mpData.Sr = STATEV(13)/100;
    mpData.Scution_Flow = STATEV(7);
    mpData.e = e;
end


end



function PS = Calculate_principal(STRESS)
%ABAQUS_PRINCIPAL_STRESS 计算ABAQUS应力分量的主应力
%   输入: STRESS_OLD - 6×1 向量，ABAQUS顺序 [S11, S22, S33, S12, S13, S23]
%   输出: PS - 3×1 向量，[大主应力, 中主应力, 小主应力]（降序排列）

% 将6个应力分量重构为3×3对称张量
% ABAQUS顺序: 11, 22, 33, 12, 13, 23
S = [STRESS(1), STRESS(4), STRESS(5);
    STRESS(4), STRESS(2), STRESS(6);
    STRESS(5), STRESS(6), STRESS(3)];

% 计算特征值（主应力）并降序排序
eigenvalues = eig(S);
PS = sort(eigenvalues, 'descend');

% 确保输出为列向量
if size(PS,1) == 1
    PS = PS';
end
end

function DDSDDE = DDSDDE_Cal(ELAM,EG2)
DDSDDE = zeros(6,6);
for I=1:3
    for J=1:3
        DDSDDE(I,J)=ELAM;
    end
    DDSDDE(I,I)=EG2+ELAM;
end
for I=4:6
    DDSDDE(I,I)=EG2/2.0;
end
end