function [Dalg, sig, epsE, epsP, STATEV] = MohrCoulomb3d_BX(...
    eps, mCst, STATEV, STATEVn, Cal_par, epsn, sign, epsEn, e0, Flow_CP, Mc, dMc, Pore_Pressure)
% 【核心-Mohr-Coulomb非饱和弹塑性本构模型】
% 基于主应力空间Return Mapping的批量向量化Mohr-Coulomb模型，考虑吸力-饱和度耦合
% 与孔隙水压力效应，适用于路基地基非饱和区/饱和区的广义有效应力分析
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 29/04/2026
% 描述:
% 进行Mohr-Coulomb弹塑性应力-应变关系计算，采用Voigt记法批量处理N个物质点。
% 本模型基于Drucker-Prager型广义有效应力框架，引入非饱和土力学中的吸力-饱和度
% 耦合效应（Bishop有效应力），并考虑孔隙水压力对骨架应力的贡献。在主应力空间
% 执行Return Mapping算法，自动判别光滑面/角点/顶点三种应力返回区域，采用向量化
% 批量处理实现大规模物质点高效计算。
%
% 理论框架：
%   屈服准则: f = (σ₁-σ₃) + (σ₁+σ₃)·sinφ - 2C·cosφ = 0  （Mohr-Coulomb六面锥）
%   流动法则: 非关联流动，剪胀角ψ控制塑性体积应变（ψ < φ 剪缩，ψ > φ 剪胀）
%   硬化规律: 理想塑性（C, φ, ψ = const），无各向同性/随动硬化
%   有效应力: σ' = σ - Sr·s·I - pw·I  （Bishop型，s为吸力，pw为孔隙水压力）
%
% 特殊处理：
%   1. 饱和区判定：通过吸力历史状态（STATEVn(:,7)）与当前吸力联合判别
%   2. 孔压效应：Pore_Pressure压为正，在有效应力中扣除（拉为正约定）
%   3. 吸力演化：基于SWCC（Van Genuchten型）由饱和度反算吸力，考虑干湿滞回
%   4. 角点处理：自动识别右角点（σ₁=σ₂）、左角点（σ₂=σ₃）、顶点（σ₁=σ₂=σ₃）
%   5. 塑性应变分解：严格区分弹性/塑性应变，避免虚假体应变（孔压不产塑性）
%
%--------------------------------------------------------------------------
% [Dalg, sig, epsE, epsP, STATEV] = MohrCoulomb3d_BX(...)
%--------------------------------------------------------------------------
% 输入:
% eps           - 当前步总应变向量 [N×6]，Voigt记法
%                 [ε₁₁, ε₂₂, ε₃₃, γ₁₂, γ₁₃, γ₂₃]（剪应变工程记法）
% mCst          - 材料参数矩阵 [N×6]，每行对应一个物质点：
%                 第1列: 弹性模量 E [Pa]
%                 第2列: 泊松比 ν [-]
%                 第3列: 粘聚力 C [Pa]
%                 第4列: 内摩擦角 φ [°]（输入为度，内部转弧度）
%                 第5列: 剪胀角 ψ [°]（输入为度，内部转弧度；ψ=0为无剪胀）
%                 第6列: 土粒比重 Gs [-]（用于饱和度-孔隙比计算）
% STATEV        - 当前步状态变量矩阵 [N×19]（输出容器，部分继承STATEVn）
% STATEVn       - 上一步状态变量矩阵 [N×19]，关键列：
%                 第7列: 上一步吸力 [kPa]（存储时除以1000，使用需乘回）
%                 第11~16列: 上一步塑性应变参考值（用于增量计算）
% Cal_par       - 求解控制参数结构体，含字段：
%                 Calculate_time - 当前计算时刻
%                 dt             - 时间步长
% epsn          - 上一步总应变 [N×6]（缺省为零）
% sign          - 上一步应力 [N×6]（缺省为零，压为负）
% epsEn         - 上一步弹性应变 [N×6]（缺省为零）
% e0            - 初始孔隙比 [N×1]（用于当前孔隙比计算）
% Flow_CP       - SWCC参数矩阵 [N×3]：
%                 第1列: Fa [kPa]（进气值相关参数）
%                 第2列: Fm [-]（VG模型参数m）
%                 第3列: Fn [-]（VG模型参数n）
% Mc            - 当前步含水质量 [N×1] [g]
% dMc           - 含水质量增量 [N×1] [g]
% Pore_Pressure - 孔隙水压力 [N×1] [Pa]（压为正，饱和区非零，非饱和区为零）
%--------------------------------------------------------------------------
% 输出:
% Dalg          - 算法一致切线模量 [6×6×N]，Voigt形式
%                 弹性点: Dalg = Dᵉ
%                 塑性点: Dalg = Dᵉ - Dᵉ·dg·(∂g/∂σ)ᵀ·Dᵉ / H（光滑面/角点/顶点区分）
% sig           - 更新后总应力 [N×6]，Voigt记法（含吸力项与孔压项贡献）
% epsE          - 更新后弹性应变 [N×6]（由净应力反推：εᵉ = Cᵉ:(σ - pw·I)）
% epsP          - 更新后塑性应变 [N×6]（εᵖ = εᵗᵒᵗᵃˡ - εᵉ，塑性点非零）
% STATEV        - 更新后状态变量 [N×19]，更新列：
%                 第9列: 等效塑性偏应变（√(2/3)·‖dev(εᵖ)‖）
%                 第10列: 塑性体应变（-tr(εᵖ)，压缩为正）
%                 第11~16列: 当前步塑性应变参考值（供下一步增量计算）
%                 第18列: σ₂₂历史值（周期时刻记录，用于K₀相关分析）
%--------------------------------------------------------------------------
% 关键算法步骤:
% 1. 吸力-饱和度计算: 基于当前含水质量Mc+dMc与孔隙比e，由SWCC反算吸力s
% 2. 弹性预测（有效应力空间）: σᵉᶠᶠ,ᵗʳ = σₙ - Srₙ·sₙ·I - pw·I + Dᵉ:Δε - Δ(Sr·s)·I
% 3. 主应力分解: 对σᵉᶠᶠ,ᵗʳ进行特征值分解，获得排序主应力σ₁≥σ₂≥σ₃及主方向Q
% 4. 屈服判别: f = (σ₁-σ₃) + (σ₁+σ₃)·sinφ - 2C·cosφ，f > tol进入Return Mapping
% 5. Return Mapping（主应力空间）:
%    - 光滑面: 单乘子γ，返回至单屈服面
%    - 右角点: σ₁=σ₂，双乘子[γ₁,γ₂]，同时满足f₁=0,f₂=0
%    - 左角点: σ₂=σ₃，双乘子[γ₁,γ₃]，同时满足f₁=0,f₃=0
%    - 顶点: σ₁=σ₂=σ₃，三向等压状态，以光滑面近似处理
% 6. 应力重构: σᵉᶠᶠ = Q·diag(σ₁,σ₂,σ₃)·Qᵀ，总应力σ = σᵉᶠᶠ + Sr·s·I + pw·I
% 7. 一致切线模量: 基于返回区域类型，构造Dalg = Dᵉ - Dᵉ·N·(Nᵀ·Dᵉ·N)⁻¹·Nᵀ·Dᵉ
% 8. 弹性应变反推: εᵉ = Cᵉ:(σ - pw·I)（确保孔压不产生虚假塑性体应变）
% 9. 塑性应变: εᵖ = ε - εᵉ
%
% 注意: 
%   - 应力符号约定: 压为负（土力学惯例），孔压Pore_Pressure输入压为正
%   - 剪应变分量: 采用工程记法γ=2ε，刚度矩阵剪切项系数为G=E/[2(1+ν)]
%   - 饱和度限制: Sr ∈ [1e-3, 1.0]，避免零饱和度数值奇异
%   - 角点容差: ct = 1e-6·max|σᵢ| + 1e-9，自动判别应力返回区域
%   - 饱和区线弹性: 若饱和区无塑性（本模型假设），需外部判定是否调用本模型
%--------------------------------------------------------------------------
% 调用关系:
% 被调用: detMPs_BX (物质点主程序，作为地基非饱和区MC本构选项)
% 内部调用: 无（主应力分解直接调用MATLAB eig，Return Mapping完全向量化）
%--------------------------------------------------------------------------
%% 0. 接口兼容
if nargin < 5 || isempty(Cal_par),  Cal_par.Calculate_time = 0; end
if nargin < 3 || isempty(STATEV),   STATEV = zeros(size(eps,1), 19); end
if nargin < 13 || isempty(Pore_Pressure), Pore_Pressure = zeros(size(eps,1), 1); end

N = size(eps,1);
if nargin < 6 || isempty(epsn),  epsn = zeros(N,6); end
if nargin < 7 || isempty(sign),   sign = zeros(N,6);  end
if nargin < 8 || isempty(epsEn), epsEn = zeros(N,6); end

tol = 1e-6;
E   = mCst(:,1);  nu = mCst(:,2);  C = mCst(:,3);
phi = mCst(:,4)*pi/180;  sf = sin(phi);  cf = cos(phi);
psi = mCst(:,5)*pi/180;  sp_ = sin(psi);
SUB_Gs = mCst(:,6);

lam = E.*nu./((1+nu).*(1-2*nu));
mu  = E./(2*(1+nu));

TIME = Cal_par.Calculate_time - Cal_par.dt;
%% 吸力与饱和度
STRAN_V = eps(:,1)+eps(:,2)+eps(:,3);
e_END = (1.0+e0).*(1+STRAN_V)-1.0;
Sr_old = min(1.0, max(1e-3, Mc.*SUB_Gs./e_END/1000));
Sr_new = min(1.0, max(1e-3, (Mc+dMc).*SUB_Gs./e_END/1000));

Fa = Flow_CP(:,1);
Fm = Flow_CP(:,2);
Fn = Flow_CP(:,3);
Suction_new = ((Sr_new.^(-1./Fm) - 1).^(1./Fn)).*Fa./e_END.^(1./(Fm.*Fn));
Suction_old = max(0.0, STATEVn(:,7)*1000.0);

%% 孔压项（拉为正：压应力取负值；Pore_Pressure压为正输入）
pw_term = (-Pore_Pressure) .* [1, 1, 1, 0, 0, 0];   % N×6，饱和区为负，非饱和区为零

%% 1. 弹性刚度 De 与柔度 Ce
factor = E./((1+nu).*(1-2*nu));
bm1 = [1;1;1;0;0;0];
De = zeros(6,6,N);
Ce = zeros(6,6,N);
for k = 1:N
    De(:,:,k) = factor(k)*( nu(k)*(bm1*bm1') + (1-2*nu(k))*diag([1,1,1,0.5,0.5,0.5]) );
    Ce_top = -nu(k)*ones(3) + (1+nu(k))*eye(3);
    Ce_bot = 2*(1+nu(k))*eye(3);
    Ce(:,:,k) = [Ce_top, zeros(3); zeros(3), Ce_bot]/E(k);
end

%% 2. 弹性预测（有效应力空间，扣吸力 + 扣孔压）
deps = eps - epsn;
deps_3 = reshape(deps', 6, 1, N);

suction_old_term = (Sr_old .* Suction_old) .* [1,1,1,0,0,0];
suction_new_term = (Sr_new .* Suction_new) .* [1,1,1,0,0,0];

% 有效应力起点：总应力 - 吸力项 - 孔压项（pw_term为负，减去负值=加绝对值）
sign_eff = sign - suction_old_term - pw_term;
dsuction_term = suction_new_term - suction_old_term;

sig_eff_tr_3 = pagemtimes(De, deps_3) + reshape(sign_eff', 6, 1, N) - reshape(dsuction_term', 6, 1, N);
sig_eff_tr = reshape(sig_eff_tr_3, 6, N)';

%% 3. 主应力分解（保留 for，3×3 eig 极快）
S = zeros(3,3,N);
S(1,1,:) = sig_eff_tr(:,1);  S(2,2,:) = sig_eff_tr(:,2);  S(3,3,:) = sig_eff_tr(:,3);
S(2,3,:) = sig_eff_tr(:,4);  S(3,2,:) = sig_eff_tr(:,4);
S(1,3,:) = sig_eff_tr(:,5);  S(3,1,:) = sig_eff_tr(:,5);
S(1,2,:) = sig_eff_tr(:,6);  S(2,1,:) = sig_eff_tr(:,6);

Q  = zeros(3,3,N);
ps = zeros(3,N);
for k = 1:N
    [V, D] = eig(S(:,:,k), 'vector');
    [ps(:,k), idx] = sort(D, 'descend');
    Q(:,:,k) = V(:,idx);
end

%% 4. 屈服判断（向量化）
f = (ps(1,:)-ps(3,:))' + (ps(1,:)+ps(3,:))'.*sf - 2*C.*cf;
plastic_mask = f > tol;
Np = sum(plastic_mask);

sig_eff = sig_eff_tr;
Dalg = De;
epsPn = epsn - epsEn;

%% 5–6. Return Mapping & sig 重构 & Dalg（掩码向量化，无 for k=1:Np）
if Np > 0
    idx_p = find(plastic_mask);
    Q_p   = Q(:,:,plastic_mask);
    ps_p  = ps(:,plastic_mask);
    C_p   = C(plastic_mask);   sf_p = sf(plastic_mask);  cf_p = cf(plastic_mask);
    sp_p  = sp_(plastic_mask); lam_p = lam(plastic_mask); mu_p = mu(plastic_mask);
    Np_loc = Np;

    % 主应力空间弹性矩阵 De3
    De3 = zeros(3,3,Np_loc);
    for k = 1:Np_loc
        De3(:,:,k) = [lam_p(k)+2*mu_p(k), lam_p(k), lam_p(k);
                      lam_p(k), lam_p(k)+2*mu_p(k), lam_p(k);
                      lam_p(k), lam_p(k), lam_p(k)+2*mu_p(k)];
    end

    % 角点容差
    maxS = max(abs(ps_p),[],1);
    ct   = 1e-6*maxS + 1e-9;

    % 统一为列向量，避免维度混乱
    s1 = ps_p(1,:)';  s2 = ps_p(2,:)';  s3 = ps_p(3,:)';
    sf_p = sf_p(:);  sp_p = sp_p(:);  C_p = C_p(:);  cf_p = cf_p(:);
    ct = ct(:);

    rc = abs(s1-s2) < ct;
    lc = abs(s2-s3) < ct;

    ps_new = ps_p;

    % ========== 光滑面（向量化）==========
    smooth_idx = find(~rc & ~lc);
    Ns = length(smooth_idx);
    if Ns > 0
        i = smooth_idx;
        ftr = (s1(i)-s3(i)) + (s1(i)+s3(i)).*sf_p(i) - 2*C_p(i).*cf_p(i);
        
        % 3×Ns：每行一个分量，每列一个点
        dg = [1+sp_p(i), zeros(Ns,1), -1+sp_p(i)]';
        df = [1+sf_p(i), zeros(Ns,1), -1+sf_p(i)]';
        
        Dedg = pagemtimes(De3(:,:,i), reshape(dg, 3, 1, Ns));
        Dedg = reshape(Dedg, 3, Ns);
        H = sum(df .* Dedg, 1);   % 1×Ns
        
        valid = abs(H) > 1e-12;
        dgam = zeros(1, Ns);
        dgam(valid) = ftr(valid)' ./ H(valid);
        
        ps_new(:,i) = ps_p(:,i) - dgam .* Dedg;
    end

    % ========== 右角点（向量化）==========
    right_idx = find(rc & ~lc);
    Nr = length(right_idx);
    if Nr > 0
        i = right_idx;
        f1 = (s1(i)-s3(i)) + (s1(i)+s3(i)).*sf_p(i) - 2*C_p(i).*cf_p(i);
        f2 = (s2(i)-s3(i)) + (s2(i)+s3(i)).*sf_p(i) - 2*C_p(i).*cf_p(i);
        F = [f1'; f2'];   % 2×Nr
        
        dg1 = [1+sp_p(i), zeros(Nr,1), -1+sp_p(i)]';
        dg2 = [zeros(Nr,1), 1+sp_p(i), -1+sp_p(i)]';
        df1 = [1+sf_p(i), zeros(Nr,1), -1+sf_p(i)]';
        df2 = [zeros(Nr,1), 1+sf_p(i), -1+sf_p(i)]';
        
        D1 = pagemtimes(De3(:,:,i), reshape(dg1, 3, 1, Nr)); D1 = reshape(D1, 3, Nr);
        D2 = pagemtimes(De3(:,:,i), reshape(dg2, 3, 1, Nr)); D2 = reshape(D2, 3, Nr);
        
        H = zeros(2,2,Nr);
        H(1,1,:) = sum(df1 .* D1, 1);
        H(1,2,:) = sum(df1 .* D2, 1);
        H(2,1,:) = sum(df2 .* D1, 1);
        H(2,2,:) = sum(df2 .* D2, 1);
        
        dg_sol = pagemldivide(H, reshape(F, 2, 1, Nr));
        dg_sol = reshape(dg_sol, 2, Nr);
        
        ps_new(:,i) = ps_p(:,i) - dg_sol(1,:).*D1 - dg_sol(2,:).*D2;
    end

    % ========== 左角点（向量化）==========
    left_idx = find(~rc & lc);
    Nl = length(left_idx);
    if Nl > 0
        i = left_idx;
        f1 = (s1(i)-s3(i)) + (s1(i)+s3(i)).*sf_p(i) - 2*C_p(i).*cf_p(i);
        f3 = (s1(i)-s2(i)) + (s1(i)+s2(i)).*sf_p(i) - 2*C_p(i).*cf_p(i);
        F = [f1'; f3'];   % 2×Nl
        
        dg1 = [1+sp_p(i), zeros(Nl,1), -1+sp_p(i)]';
        dg3 = [1+sp_p(i), -1+sp_p(i), zeros(Nl,1)]';
        df1 = [1+sf_p(i), zeros(Nl,1), -1+sf_p(i)]';
        df3 = [1+sf_p(i), -1+sf_p(i), zeros(Nl,1)]';
        
        D1 = pagemtimes(De3(:,:,i), reshape(dg1, 3, 1, Nl)); D1 = reshape(D1, 3, Nl);
        D3 = pagemtimes(De3(:,:,i), reshape(dg3, 3, 1, Nl)); D3 = reshape(D3, 3, Nl);
        
        H = zeros(2,2,Nl);
        H(1,1,:) = sum(df1 .* D1, 1);
        H(1,2,:) = sum(df1 .* D3, 1);
        H(2,1,:) = sum(df3 .* D1, 1);
        H(2,2,:) = sum(df3 .* D3, 1);
        
        dg_sol = pagemldivide(H, reshape(F, 2, 1, Nl));
        dg_sol = reshape(dg_sol, 2, Nl);
        
        ps_new(:,i) = ps_p(:,i) - dg_sol(1,:).*D1 - dg_sol(2,:).*D3;
    end

    % ========== 顶点（向量化，光滑面兜底）==========
    apex_idx = find(rc & lc);
    Na = length(apex_idx);
    if Na > 0
        i = apex_idx;
        ftr = (s1(i)-s3(i)) + (s1(i)+s3(i)).*sf_p(i) - 2*C_p(i).*cf_p(i);
        
        dg = [1+sp_p(i), zeros(Na,1), -1+sp_p(i)]';
        df = [1+sf_p(i), zeros(Na,1), -1+sf_p(i)]';
        
        Dedg = pagemtimes(De3(:,:,i), reshape(dg, 3, 1, Na));
        Dedg = reshape(Dedg, 3, Na);
        H = sum(df .* Dedg, 1);
        
        valid = abs(H) > 1e-12;
        dgam = zeros(1, Na);
        dgam(valid) = ftr(valid)' ./ H(valid);
        
        ps_new(:,i) = ps_p(:,i) - dgam .* Dedg;
    end

    % ========== 重构 sig（pagemtimes 向量化）==========
    Lambda = zeros(3,3,Np_loc);
    for k = 1:3
        Lambda(k,k,:) = ps_new(k,:);
    end
    temp = pagemtimes(Q_p, Lambda);
    temp = pagemtimes(temp, 'none', Q_p, 'ctranspose');
    
    sig_eff_p = [reshape(temp(1,1,:),[],1), reshape(temp(2,2,:),[],1), reshape(temp(3,3,:),[],1), ...
                 reshape(temp(2,3,:),[],1), reshape(temp(1,3,:),[],1), reshape(temp(1,2,:),[],1)];
    sig_eff(plastic_mask,:) = sig_eff_p;

    % ========== Dalg（向量化）==========
    n1 = reshape(Q_p(:,1,:), 3, Np_loc);
    n2 = reshape(Q_p(:,2,:), 3, Np_loc);
    n3 = reshape(Q_p(:,3,:), 3, Np_loc);
    
    N1 = [n1(1,:).^2; n1(2,:).^2; n1(3,:).^2; n1(2,:).*n1(3,:); n1(1,:).*n1(3,:); n1(1,:).*n1(2,:)];
    N2 = [n2(1,:).^2; n2(2,:).^2; n2(3,:).^2; n2(2,:).*n2(3,:); n2(1,:).*n2(3,:); n2(1,:).*n2(2,:)];
    N3 = [n3(1,:).^2; n3(2,:).^2; n3(3,:).^2; n3(2,:).*n3(3,:); n3(1,:).*n3(3,:); n3(1,:).*n2(2,:)];

    rc = abs(ps_new(1,:)-ps_new(2,:)) < ct';
    lc = abs(ps_new(2,:)-ps_new(3,:)) < ct';
    smooth_mask = ~rc & ~lc;
    right_mask  = rc & ~lc;
    left_mask   = ~rc & lc;
    apex_mask   = rc & lc;

    % --- 光滑面 ---
    if any(smooth_mask)
        ks = find(smooth_mask);
        idx_s = idx_p(ks);
        Ns = length(ks);
        
        sf_p_ks = sf_p(ks)';  sp_p_ks = sp_p(ks)';
        df_s = (1+sf_p_ks).*N1(:,ks) + (-1+sf_p_ks).*N3(:,ks);
        dg_s = (1+sp_p_ks).*N1(:,ks) + (-1+sp_p_ks).*N3(:,ks);
        
        Dedg_3 = pagemtimes(De(:,:,idx_s), reshape(dg_s, 6, 1, Ns));
        Dedg_s = reshape(Dedg_3, 6, Ns);
        H_s = sum(df_s .* Dedg_s, 1);
        valid = abs(H_s) > 1e-12;
        
        if any(valid)
            v = find(valid);
            idx_v = idx_s(v);  Nv = length(v);
            
            Dedg_v = Dedg_s(:,v);
            df_v = df_s(:,v);
            H_v = reshape(H_s(v), 1, []);
            
            Dedg_3v = reshape(Dedg_v, 6, 1, Nv);
            df_3v   = reshape(df_v, 1, 6, Nv);
            T = pagemtimes(Dedg_3v, df_3v);
            T = T ./ reshape(H_v, 1, 1, Nv);
            
            De_v = De(:,:,idx_v);
            Dalg(:,:,idx_v) = De_v - pagemtimes(T, De_v);
        end
    end

    % --- 右角点（向量化，显式 2×2 逆）---
    if any(right_mask)
        kr = find(right_mask);
        idx_r = idx_p(kr);
        Nr = length(kr);
        
        sf_p_kr = sf_p(kr)';  sp_p_kr = sp_p(kr)';
        df1 = (1+sf_p_kr).*N1(:,kr) + (-1+sf_p_kr).*N3(:,kr);
        df2 = (1+sf_p_kr).*N2(:,kr) + (-1+sf_p_kr).*N3(:,kr);
        dg1 = (1+sp_p_kr).*N1(:,kr) + (-1+sp_p_kr).*N3(:,kr);
        dg2 = (1+sp_p_kr).*N2(:,kr) + (-1+sp_p_kr).*N3(:,kr);
        
        D1 = pagemtimes(De(:,:,idx_r), reshape(dg1, 6, 1, Nr)); D1 = reshape(D1, 6, Nr);
        D2 = pagemtimes(De(:,:,idx_r), reshape(dg2, 6, 1, Nr)); D2 = reshape(D2, 6, Nr);
        
        H = zeros(2,2,Nr);
        H(1,1,:) = sum(df1 .* D1, 1);
        H(1,2,:) = sum(df1 .* D2, 1);
        H(2,1,:) = sum(df2 .* D1, 1);
        H(2,2,:) = sum(df2 .* D2, 1);
        
        detH = H(1,1,:).*H(2,2,:) - H(1,2,:).*H(2,1,:);
        iH = zeros(2,2,Nr);
        iH(1,1,:) = H(2,2,:)./detH;
        iH(1,2,:) = -H(1,2,:)./detH;
        iH(2,1,:) = -H(2,1,:)./detH;
        iH(2,2,:) = H(1,1,:)./detH;
        
        D_mat = cat(2, reshape(D1,6,1,Nr), reshape(D2,6,1,Nr));
        df_mat = cat(1, reshape(df1,1,6,Nr), reshape(df2,1,6,Nr));
        
        T = pagemtimes(pagemtimes(D_mat, iH), df_mat);
        De_r = De(:,:,idx_r);
        Dalg(:,:,idx_r) = De_r - pagemtimes(T, De_r);
    end

    % --- 左角点（向量化）---
    if any(left_mask)
        kl = find(left_mask);
        idx_l = idx_p(kl);
        Nl = length(kl);
        
        sf_p_kl = sf_p(kl)';  sp_p_kl = sp_p(kl)';
        df1 = (1+sf_p_kl).*N1(:,kl) + (-1+sf_p_kl).*N3(:,kl);
        df3 = (1+sf_p_kl).*N1(:,kl) + (-1+sf_p_kl).*N2(:,kl);
        dg1 = (1+sp_p_kl).*N1(:,kl) + (-1+sp_p_kl).*N3(:,kl);
        dg3 = (1+sp_p_kl).*N1(:,kl) + (-1+sp_p_kl).*N2(:,kl);
        
        D1 = pagemtimes(De(:,:,idx_l), reshape(dg1, 6, 1, Nl)); D1 = reshape(D1, 6, Nl);
        D3 = pagemtimes(De(:,:,idx_l), reshape(dg3, 6, 1, Nl)); D3 = reshape(D3, 6, Nl);
        
        H = zeros(2,2,Nl);
        H(1,1,:) = sum(df1 .* D1, 1);
        H(1,2,:) = sum(df1 .* D3, 1);
        H(2,1,:) = sum(df3 .* D1, 1);
        H(2,2,:) = sum(df3 .* D3, 1);
        
        detH = H(1,1,:).*H(2,2,:) - H(1,2,:).*H(2,1,:);
        iH = zeros(2,2,Nl);
        iH(1,1,:) = H(2,2,:)./detH;
        iH(1,2,:) = -H(1,2,:)./detH;
        iH(2,1,:) = -H(2,1,:)./detH;
        iH(2,2,:) = H(1,1,:)./detH;
        
        D_mat = cat(2, reshape(D1,6,1,Nl), reshape(D3,6,1,Nl));
        df_mat = cat(1, reshape(df1,1,6,Nl), reshape(df3,1,6,Nl));
        
        T = pagemtimes(pagemtimes(D_mat, iH), df_mat);
        De_l = De(:,:,idx_l);
        Dalg(:,:,idx_l) = De_l - pagemtimes(T, De_l);
    end

    % --- 顶点（向量化）---
    if any(apex_mask)
        ka = find(apex_mask);
        idx_a = idx_p(ka);
        Na = length(ka);
        
        sf_p_ka = sf_p(ka)';  sp_p_ka = sp_p(ka)';
        df_a = (1+sf_p_ka).*N1(:,ka) + (-1+sf_p_ka).*N3(:,ka);
        dg_a = (1+sp_p_ka).*N1(:,ka) + (-1+sp_p_ka).*N3(:,ka);
        
        Dedg_3 = pagemtimes(De(:,:,idx_a), reshape(dg_a, 6, 1, Na));
        Dedg_a = reshape(Dedg_3, 6, Na);
        H_a = sum(df_a .* Dedg_a, 1);
        valid = abs(H_a) > 1e-12;
        
        if any(valid)
            v = find(valid);
            idx_v = idx_a(v);  Nv = length(v);
            
            Dedg_v = Dedg_a(:,v);
            df_v = df_a(:,v);
            H_v = reshape(H_a(v), 1, []);
            
            Dedg_3v = reshape(Dedg_v, 6, 1, Nv);
            df_3v   = reshape(df_v, 1, 6, Nv);
            T = pagemtimes(Dedg_3v, df_3v);
            T = T ./ reshape(H_v, 1, 1, Nv);
            
            De_v = De(:,:,idx_v);
            Dalg(:,:,idx_v) = De_v - pagemtimes(T, De_v);
        end
    end
end

%% 7. 弹性应变与塑性应变（核心：无虚假体应变）
% 输出 sig 为总应力 = 有效应力 + 吸力项 + 孔压项
sig = sig_eff + suction_new_term + pw_term;

if Np == 0
    % 弹性步：塑性应变严格保持历史，弹性应变由总应变分解
    epsP = epsPn;
    epsE = eps - epsP;
else
    % 塑性步：epsE 必须由骨架净应力反推（总应力 - 孔压）
    % 净应力 = sig_eff + suction_new_term = sig - pw_term
    % Ce 是骨架柔度，只能响应骨架真实应力，不能含孔压
    sig_net = sig - pw_term;           % 扣除孔压后的骨架应力
    sig_net_3 = reshape(sig_net', 6, 1, N);
    epsE_3 = pagemtimes(Ce, sig_net_3);
    epsE = reshape(epsE_3, 6, N)';
    epsP = eps - epsE;
end

%% 8. STATEV
if isfield(Cal_par, 'Calculate_time') && Cal_par.Calculate_time > 1.0
    epsP_ref = STATEV(:,11:16);
    epsP_inc = epsP - epsP_ref;
    STATEV(:,10) = -sum(epsP_inc(:,1:3), 2);
    meanP = sum(epsP_inc(:,1:3), 2) / 3;
    dev = epsP_inc;
    dev(:,1:3) = epsP_inc(:,1:3) - meanP;
    STATEV(:,9) = sqrt( (2/3) * ( sum(dev(:,1:3).^2, 2) + 0.5*sum(epsP_inc(:,4:6).^2, 2) ) );
else
    STATEV(:,11:16) = epsP;
    STATEV(:,9)  = zeros(N,1);
    STATEV(:,10) = zeros(N,1);
end

if (abs(TIME-round(TIME)) < 1e-6)
    STATEV(:,18)= sig(:,2);
end

end