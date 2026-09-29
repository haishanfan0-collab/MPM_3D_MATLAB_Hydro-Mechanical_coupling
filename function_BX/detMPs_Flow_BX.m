function [fint,mpData,mesh] = detMPs_Flow_BX(H,mpData,mesh)
% 【核心-渗流场物质点本构与矩阵组装代码】，内置水头梯度计算、饱和度-吸力更新、
% 渗透系数张量构造、降雨入渗流量计算、容水度矩阵组装、整体渗流刚度/容量矩阵稀疏拼接
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行物质点法渗流场核心计算，涵盖节点水头插值获取物质点水头梯度（类比应变）、
% 基于 VG 模型计算饱和度与吸力、构造各向异性渗透系数张量、计算降雨入渗边界流量、
% 组装容水度矩阵（类比阻尼矩阵）与渗透刚度矩阵（类比切线刚度矩阵）、节点内流量
% 散度计算、整体稀疏矩阵拼接全过程，为湿-力耦合分析中渗流场求解的核心子函数。
%
%--------------------------------------------------------------------------
% [fint, mpData, mesh] = detMPs_Flow_BX(H, mpData, mesh)
%--------------------------------------------------------------------------
% 输入:
% H        - 背景网格节点总水头向量 [nDoF×1]
% mpData   - 物质点信息结构体数组，含 dSvp（形函数梯度）、Svp（形函数值）、
%            Flow_CP（VG参数 [Fa,Fm,Fn]）、Ksat（饱和渗透系数）、vp（体积）、
%            e（孔隙比）、H0（初始水头）、dH_sum（水头累积）、dH_sum_strain
%            （应变致水头变化）、Flow_P（降雨强度）、Flow_SIZE（特征尺寸）等
% mesh     - 背景网格结构体，输出时嵌入 Flow.Kt（渗流刚度）与 Flow.Ct（容水度）
%--------------------------------------------------------------------------
% 输出:
% fint     - 节点内流量向量 [nDoF×1]，类比力学中的节点内力
% mpData   - 更新后的物质点信息结构体数组，含 eps_Flow（水头梯度）、sig_Flow
%            （渗流应力 K∇H）、Sr（饱和度）、Saturated_YES（饱和标记）、
%            Flow_P_current（实际入渗流量）等
% mesh     - 更新后的网格结构体，mesh.Flow.Kt 为整体渗透刚度矩阵，
%            mesh.Flow.Ct 为整体容水度矩阵（均稀疏）
%--------------------------------------------------------------------------
% 关键算法:
%   1. 水头梯度: eps_Flow = G * H(ed)，通过 cumsum 前缀和向量化提取
%   2. 吸力计算: Suction = -(H0 + dH_sum + mpH + dH_sum_strain) * 9.81
%   3. 饱和度 VG 模型: Sr = [1 + (Suction/coeff_swrc)^n]^{-m}
%   4. 渗透系数: K_mp = diag(Ksat) * kr(Sr)，各向异性张量
%   5. 降雨入渗: Flow_P_current = min(Flow_P, k_SLOW * L * eps_FlownEtr0)
%   6. 容水度矩阵: cp = diag(sum(vp*c_scalar*(N*N'), 2))，行求和 lumping
%   7. 渗透刚度: kp = vp * (G' * K_mp * G)
%   8. 节点流量: fp = vp * (G' * sig_Flow)
%   9. 稀疏组装: mesh.Flow.Kt/Ct = sparse(krow, kcol, kval/cval)
%--------------------------------------------------------------------------

nn_all = [mpData.nn]';
nmp   = length(mpData);                                                     % number of material points
fint  = zeros(size(H));                                                     % zero internal force vector
nD = length(mpData(1).mpC);

npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]/nD^2);                                            % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow; cval=krow;                    % zero the stiffness information
eps_FlownEtr = [mpData.epsn_Flow]';

%% 尝试去除 for 循环
ed = [mpData.nIN];        % 物质点关联的节点（类比位移求解的节点集）
G_all = [mpData.dSvp];        % 形函数导数 ∇N (nD × nn)（类比位移求解的B矩阵）
vp = [mpData.vp]';         % 物质点体积
nn_sizeall = [mpData.nn_size]';

% ------------------- 1. 水头增量 + 历史累加（核心：位移类比水头） -------------------
% 结果：values(j) = uvw( ed_all(j) )
values = H(ed);  % 189796 x 1
% 逐元素加权（广播机制）
% G_scaled(:,j) = G_all(:,j) * values(j)
G_scaled = G_all .* values';  % 2 x 189796
% 列方向累积求和（前缀和）
Prefix = cumsum(G_scaled, 2);  % 2 x 189796
% 通过索引相减提取各段结果（向量化，无循环）
ss = nn_sizeall(:,1);
ee = nn_sizeall(:,2);
% 初始化结果矩阵 4 x 23530
Result = Prefix(:, ee);
% 处理起始位置大于1的段：减去前一个前缀和
% 利用逻辑掩码处理 s=1 的边界情况
mask = (ss > 1);
if any(mask)
    s_prev = ss - 1;
    Result(:, mask) = Result(:, mask) - Prefix(:, s_prev(mask));
end

eps_FlownEtr = eps_FlownEtr + Result';  % 累加至历史值

% 当前饱和度
LS = reshape([mpData.Flow_CP]',3,nmp)';
Flowpar_a = LS(:,1)/1000;
Flowpar_m = LS(:,2);
Flowpar_n = LS(:,3);

% 计算 mpH = N*H(ed);
N_all   = [mpData.Svp]';
mpH = N_all.*H(ed);
% 列方向累积求和（前缀和）
Prefix = cumsum(mpH, 1);  % 189796*1
% 初始化结果矩阵 4 x 23530
Result = Prefix(ee);
% 处理起始位置大于1的段：减去前一个前缀和
% 利用逻辑掩码处理 s=1 的边界情况
mask = (ss > 1);
if any(mask)
    s_prev = ss - 1;
    Result(mask) = Result(mask) - Prefix(s_prev(mask));
end
mpH = Result;

Suction_old =    ([mpData.H0]' + [mpData.dH_sum]' + [mpData.dH_sum_strain]')*-9.81;
Suction_new = max(([mpData.H0]' + [mpData.dH_sum]' + mpH + [mpData.dH_sum_strain]')*-9.81,0.00);

% 计算SWRC系数: m3 / e^{1/(m1*m2)}
e = [mpData.e]';
coeff_swrc = Flowpar_a ./ (e.^(1./(Flowpar_m.*Flowpar_n)));

% Sr = [1 + (suction / A)^{n}]^{-m}  ← 注意是除法
Sr = (1.0 + (Suction_new ./ coeff_swrc).^Flowpar_n).^(-Flowpar_m);
Sr =  min(1.0,Sr);

%% 完全饱和判定，饱和区域容水度为零
Saturated_YES = [mpData.Saturated_YES]';
Saturated_YES(Suction_new <= 0.000) = 1.0;
Saturated_YES(Suction_new > 0.000) = 0.0;

%% 相对渗透系数（VG模型）
kr = max(0.01,Sr.^0.5 .*(1-(1-Sr.^(1./Flowpar_m)).^Flowpar_m).^2);

% 渗透系数张量（类比弹性模量，负号适配达西定律，保留你的对角矩阵）
Ksat0 = [mpData.Ksat]';

Ksat = zeros(nD, nD, nmp);
Ksat(1,1,:) = Ksat0(:,1);  % 第一列放到 (1,1) 位置
Ksat(2,2,:) = Ksat0(:,2);  % 第二列放到 (2,2) 位置
if nD == 3
    Ksat(3,3,:) = Ksat0(:,3);  % 第二列放到 (2,2) 位置
end
kr = reshape(kr, 1, 1, nmp);
K_mp = pagemtimes(Ksat,kr);

%% 实际降雨渗透流量
Flow_P = [mpData.Flow_P]';
Flow_P_current = zeros(nmp,1);
Flow_P_YES = Flow_P~=0;
L=[mpData(Flow_P_YES).Flow_P]'./[mpData(Flow_P_YES).Flow_SIZE]';
if nD == 3
    Len = L.^0.5;
elseif nD == 2
    Len = L;
end

k_SLOW = max(0.01, [mpData(Flow_P_YES).Sr]'.^0.5 .* ...
    (1-(1-[mpData(Flow_P_YES).Sr]'.^(1./Flowpar_m(Flow_P_YES))).^Flowpar_m(Flow_P_YES)).^2) .* ...
    Ksat0(Flow_P_YES,2);

FHS_LSLS = reshape([mpData(Flow_P_YES).lp],nD,sum(Flow_P_YES))';

eps_FlownEtr0 = max(-1,min(100,2.*Suction_old(Flow_P_YES)./9.81./min(Len,FHS_LSLS(:,2))+1));
% 当前水头
Flow_P_current(Flow_P_YES)=min(Flow_P(Flow_P_YES),k_SLOW.*L.*eps_FlownEtr0);

%% 持水度计算 = e/(1+e) * (m*n/A) * Sr^{...} * (...)  ← 注意是除以A，饱和区域容水度为零
c_scalar = zeros(nmp,1);
c_scalar(~Saturated_YES) = (e(~Saturated_YES)./(e(~Saturated_YES)+1)) ...
    .* (Flowpar_m(~Saturated_YES) .*...
    Flowpar_n(~Saturated_YES) ./ coeff_swrc(~Saturated_YES)) .* ...
    Sr(~Saturated_YES).^((Flowpar_m(~Saturated_YES)+1)...
    ./Flowpar_m(~Saturated_YES)) .* (Sr(~Saturated_YES).^(-1./Flowpar_m(~Saturated_YES)) - 1)...
    .^((Flowpar_n(~Saturated_YES)-1)./Flowpar_n(~Saturated_YES));

%% 类比：固体力学内力=B^TσV → 渗流流量=G^T*(K∇H)*V
l0 = reshape(eps_FlownEtr', nD, 1, nmp);
mpV = pagemtimes(K_mp, l0);  % 类比应力σ（渗流中为K∇H，即"渗流应力"）
mpV =squeeze(mpV);

vp_3d = reshape(vp,1,1,nmp);
mpV_3d = reshape(mpV,nD,1,nmp);
vpc = vp .* c_scalar;
vpc_3d = reshape(vpc,1,1,nmp);

%% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 4
if nD ==2
    max_4_YES = (nn_all == 4);
    if sum(max_4_YES)>=4000
        cols_3d_4 = [ss(max_4_YES),ss(max_4_YES)+1,ss(max_4_YES)+2,ss(max_4_YES)+3];
        G_3d_4 = reshape(G_all(:, cols_3d_4'), nD, 4, sum(max_4_YES));
        % fp计算
        temp = pagemtimes(vp_3d(:,:,max_4_YES),'none',G_3d_4,'ctranspose');
        fp_4 = pagemtimes(temp,mpV_3d(:,:,max_4_YES));

        % cp计算
        N_3d_4 = reshape(N_all(cols_3d_4',:)', 1, 4, sum(max_4_YES));
        temp = pagemtimes(N_3d_4,'ctranspose',N_3d_4,'none');
        cp_temp = pagemtimes(vpc_3d(:,:,max_4_YES),temp);
        cp_temp = sum(cp_temp, 2); % 行求和 lumping

        cp_4 = zeros(4,4,sum(max_4_YES));
        cp_4(1,1,:)=cp_temp(1,1,:);
        cp_4(2,2,:)=cp_temp(2,1,:);
        cp_4(3,3,:)=cp_temp(3,1,:);
        cp_4(4,4,:)=cp_temp(4,1,:);

        % kp 计算
        temp = pagemtimes(vp_3d(:,:,max_4_YES),'none',G_3d_4,'ctranspose');
        temp = pagemtimes(temp,K_mp(:,:,max_4_YES));
        kp_4 = pagemtimes(temp,G_3d_4);
    else
        max_4_YES = (nn_all == -1);
    end

    %% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 6
    max_6_YES = (nn_all == 6);
    if sum(max_6_YES)>=2000
        cols_3d_6 = [ss(max_6_YES),ss(max_6_YES)+1,ss(max_6_YES)+2,ss(max_6_YES)+3,...
            ss(max_6_YES)+4,ss(max_6_YES)+5];
        G_3d_6 = reshape(G_all(:, cols_3d_6'), nD, 6, sum(max_6_YES));
        % fp计算
        temp = pagemtimes(vp_3d(:,:,max_6_YES),'none',G_3d_6,'ctranspose');
        fp_6 = pagemtimes(temp,mpV_3d(:,:,max_6_YES));

        % cp计算
        N_3d_6 = reshape(N_all(cols_3d_6',:)', 1, 6, sum(max_6_YES));
        temp = pagemtimes(N_3d_6,'ctranspose',N_3d_6,'none');
        cp_temp = pagemtimes(vpc_3d(:,:,max_6_YES),temp);
        cp_temp = sum(cp_temp, 2); % 行求和 lumping

        cp_6 = zeros(6,6,sum(max_6_YES));
        cp_6(1,1,:)=cp_temp(1,1,:);
        cp_6(2,2,:)=cp_temp(2,1,:);
        cp_6(3,3,:)=cp_temp(3,1,:);
        cp_6(4,4,:)=cp_temp(4,1,:);
        cp_6(5,5,:)=cp_temp(5,1,:);
        cp_6(6,6,:)=cp_temp(6,1,:);

        % kp 计算
        temp = pagemtimes(vp_3d(:,:,max_6_YES),'none',G_3d_6,'ctranspose');
        temp = pagemtimes(temp,K_mp(:,:,max_6_YES));
        kp_6 = pagemtimes(temp,G_3d_6);
    else
        max_6_YES = (nn_all == -1);
    end

    %% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 9
    max_9_YES = (nn_all == 9);
    if sum(max_9_YES)>=1000
        cols_3d_9 = [ss(max_9_YES),ss(max_9_YES)+1,ss(max_9_YES)+2,ss(max_9_YES)+3,...
            ss(max_9_YES)+4,ss(max_9_YES)+5,ss(max_9_YES)+6,ss(max_9_YES)+7,...
            ss(max_9_YES)+8];
        G_3d_9 = reshape(G_all(:, cols_3d_9'), nD, 9, sum(max_9_YES));
        % fp计算
        temp = pagemtimes(vp_3d(:,:,max_9_YES),'none',G_3d_9,'ctranspose');
        fp_9 = pagemtimes(temp,mpV_3d(:,:,max_9_YES));

        % cp计算
        N_3d_9 = reshape(N_all(cols_3d_9',:)', 1, 9, sum(max_9_YES));
        temp = pagemtimes(N_3d_9,'ctranspose',N_3d_9,'none');
        cp_temp = pagemtimes(vpc_3d(:,:,max_9_YES),temp);
        cp_temp = sum(cp_temp, 2); % 行求和 lumping

        cp_9 = zeros(9,9,sum(max_9_YES));
        cp_9(1,1,:)=cp_temp(1,1,:);
        cp_9(2,2,:)=cp_temp(2,1,:);
        cp_9(3,3,:)=cp_temp(3,1,:);
        cp_9(4,4,:)=cp_temp(4,1,:);
        cp_9(5,5,:)=cp_temp(5,1,:);
        cp_9(6,6,:)=cp_temp(6,1,:);
        cp_9(7,7,:)=cp_temp(7,1,:);
        cp_9(8,8,:)=cp_temp(8,1,:);
        cp_9(9,9,:)=cp_temp(9,1,:);

        % kp 计算
        temp = pagemtimes(vp_3d(:,:,max_9_YES),'none',G_3d_9,'ctranspose');
        temp = pagemtimes(temp,K_mp(:,:,max_9_YES));
        kp_9 = pagemtimes(temp,G_3d_9);
    else
        max_9_YES = (nn_all == -1);
    end
else
    max_4_YES = (nn_all == 8);
    if sum(max_4_YES)>=4000
        cols_3d_4 = [ss(max_4_YES),ss(max_4_YES)+1,ss(max_4_YES)+2,ss(max_4_YES)+3,...
                     ss(max_4_YES)+4,ss(max_4_YES)+5,ss(max_4_YES)+6,ss(max_4_YES)+7];
        G_3d_4 = reshape(G_all(:, cols_3d_4'), nD, 8, sum(max_4_YES));
        % fp计算
        temp = pagemtimes(vp_3d(:,:,max_4_YES),'none',G_3d_4,'ctranspose');
        fp_4 = pagemtimes(temp,mpV_3d(:,:,max_4_YES));

        % cp计算
        N_3d_4 = reshape(N_all(cols_3d_4',:)', 1, 8, sum(max_4_YES));
        temp = pagemtimes(N_3d_4,'ctranspose',N_3d_4,'none');
        cp_temp = pagemtimes(vpc_3d(:,:,max_4_YES),temp);
        cp_temp = sum(cp_temp, 2); % 行求和 lumping

        cp_4 = zeros(8,8,sum(max_4_YES));
        for FFi = 1:8
            cp_4(FFi,FFi,:)=cp_temp(FFi,1,:);
        end

        % kp 计算
        temp = pagemtimes(vp_3d(:,:,max_4_YES),'none',G_3d_4,'ctranspose');
        temp = pagemtimes(temp,K_mp(:,:,max_4_YES));
        kp_4 = pagemtimes(temp,G_3d_4);
    else
        max_4_YES = (nn_all == -1);
    end

    %% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 6
    max_6_YES = (nn_all == 12);
    if sum(max_6_YES)>=2000
        cols_3d_6 = [ss(max_6_YES),ss(max_6_YES)+1,ss(max_6_YES)+2,ss(max_6_YES)+3,...
                     ss(max_6_YES)+4,ss(max_6_YES)+5,ss(max_6_YES)+6,ss(max_6_YES)+7,...
                     ss(max_6_YES)+8,ss(max_6_YES)+9,ss(max_6_YES)+10,ss(max_6_YES)+11];
        G_3d_6 = reshape(G_all(:, cols_3d_6'), nD, 12, sum(max_6_YES));
        % fp计算
        temp = pagemtimes(vp_3d(:,:,max_6_YES),'none',G_3d_6,'ctranspose');
        fp_6 = pagemtimes(temp,mpV_3d(:,:,max_6_YES));

        % cp计算
        N_3d_6 = reshape(N_all(cols_3d_6',:)', 1, 12, sum(max_6_YES));
        temp = pagemtimes(N_3d_6,'ctranspose',N_3d_6,'none');
        cp_temp = pagemtimes(vpc_3d(:,:,max_6_YES),temp);
        cp_temp = sum(cp_temp, 2); % 行求和 lumping

        cp_6 = zeros(12,12,sum(max_6_YES));
        for FFi = 1:12
            cp_6(FFi,FFi,:)=cp_temp(FFi,1,:);
        end
        
        % kp 计算
        temp = pagemtimes(vp_3d(:,:,max_6_YES),'none',G_3d_6,'ctranspose');
        temp = pagemtimes(temp,K_mp(:,:,max_6_YES));
        kp_6 = pagemtimes(temp,G_3d_6);
    else
        max_6_YES = (nn_all == -1);
    end

    %% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 9
    max_9_YES = (nn_all == 27);
    if sum(max_9_YES)>=1000
        cols_3d_9 = [ss(max_9_YES),ss(max_9_YES)+1,ss(max_9_YES)+2,ss(max_9_YES)+3,...
                     ss(max_9_YES)+4,ss(max_9_YES)+5,ss(max_9_YES)+6,ss(max_9_YES)+7,...
                     ss(max_9_YES)+8,ss(max_9_YES)+9,ss(max_9_YES)+10,ss(max_9_YES)+11,...
                     ss(max_9_YES)+12,ss(max_9_YES)+13,ss(max_9_YES)+14,ss(max_9_YES)+15,...
                     ss(max_9_YES)+16,ss(max_9_YES)+17,ss(max_9_YES)+18,ss(max_9_YES)+19,...
                     ss(max_9_YES)+20,ss(max_9_YES)+21,ss(max_9_YES)+22,ss(max_9_YES)+23,...
                     ss(max_9_YES)+24,ss(max_9_YES)+25,ss(max_9_YES)+26];
        G_3d_9 = reshape(G_all(:, cols_3d_9'), nD, 27, sum(max_9_YES));
        % fp计算
        temp = pagemtimes(vp_3d(:,:,max_9_YES),'none',G_3d_9,'ctranspose');
        fp_9 = pagemtimes(temp,mpV_3d(:,:,max_9_YES));

        % cp计算
        N_3d_9 = reshape(N_all(cols_3d_9',:)', 1, 27, sum(max_9_YES));
        temp = pagemtimes(N_3d_9,'ctranspose',N_3d_9,'none');
        cp_temp = pagemtimes(vpc_3d(:,:,max_9_YES),temp);
        cp_temp = sum(cp_temp, 2); % 行求和 lumping

        cp_9 = zeros(27,27,sum(max_9_YES));
        for FFi = 1:27
            cp_9(FFi,FFi,:)=cp_temp(FFi,1,:);
        end
        % kp 计算
        temp = pagemtimes(vp_3d(:,:,max_9_YES),'none',G_3d_9,'ctranspose');
        temp = pagemtimes(temp,K_mp(:,:,max_9_YES));
        kp_9 = pagemtimes(temp,G_3d_9);
    else
        max_9_YES = (nn_all == -1);
    end

end
% disp(['4-',num2str(sum((nn_all == 4))), ' 5-', num2str(sum((nn_all == 5))), ...
%      ' 6-',num2str(sum((nn_all == 6))), ' 7-', num2str(sum((nn_all == 7))), ...
%      ' 8-',num2str(sum((nn_all == 8))), ' 9-', num2str(sum((nn_all == 9))),])

%% 循环赋值，仅内存拷贝
k_4_JS=1;  k_6_JS=1; k_9_JS=1;

for mp=1:nmp                                                                % material point loop
    if max_4_YES(mp)
        cp = cp_4(:,:,k_4_JS); kp = kp_4(:,:,k_4_JS); fp = fp_4(:,k_4_JS); k_4_JS = k_4_JS+1;
    elseif max_6_YES(mp)
        cp = cp_6(:,:,k_6_JS); kp = kp_6(:,:,k_6_JS); fp = fp_6(:,k_6_JS); k_6_JS = k_6_JS+1;
    elseif max_9_YES(mp)
        cp = cp_9(:,:,k_9_JS); kp = kp_9(:,:,k_9_JS); fp = fp_9(:,k_9_JS); k_9_JS = k_9_JS+1;
    else
        %----------------------------------------------------------------------
        cols = ss(mp):ee(mp);                    % 当前块列范围
        % 左乘：dXdx(:,:,i) * G_all(:,cols)
        G = G_all(:, cols);
        N = N_all(cols);
        fp = vp(mp) * (G' * mpV(:,mp));  % 节点流量（类比节点内力，保留你的detF）
        % ------------------- 4. 容水度矩阵（类比阻尼矩阵，保留你的公式） -------------------
        % 第二步：行求和（集中质量法核心，只加这1行！）
        % 或者一致形式（方式2，但通常需要 lumping 简化）
        cp_consistent = vp(mp) * c_scalar(mp) * (N * N');   % nn×nn矩阵
        cp = diag(sum(cp_consistent, 2)); % 行求和 lumping
        % ------------------- 5. 刚度矩阵（类比切线刚度，保留你的组装） -------------------
        kp = vp(mp) * (G' * K_mp(:,:,mp) * G);  % 渗透刚度（类比固体力学切线刚度）
    end
    % ------------------- 6. 矩阵组装（完全保留你的逻辑） -------------------
    ed0=ed(ss(mp):ee(mp));
    npDoF=(nn_all(mp))^2;                                                % no. entries in kp
    nnDoF=nn_all(mp);
              
    krow(npCnt+1:npCnt+npDoF)=repmat(ed0.',nnDoF,1);         
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed0  ,nnDoF,1);         
    kval(npCnt+1:npCnt+npDoF)=kp;           % 渗透刚度（类比切线刚度）
    cval(npCnt+1:npCnt+npDoF)=cp;     % 容水度（类比阻尼矩阵）
    fint(ed0)=fint(ed0)+fp;                   % 内流量（类比节点内力）
    npCnt=npCnt+npDoF;                      

    % ------------------- 7. 保存当前状态（类比位移求解的应力/应变保存） -------------------
    mpData(mp).eps_Flow  = eps_FlownEtr(mp,:)';    % 本次水头增量（类比应变增量）
    mpData(mp).sig_Flow  = mpV(:,mp);              % 本次"渗流应力"（K∇H，类比Cauchy应力）
    mpData(mp).Flow_P_current = Flow_P_current(mp);
    mpData(mp).Saturated_YES = Saturated_YES(mp);
    mpData(mp).Sr = Sr(mp);  % 添加这一行
end

% 构建全局矩阵（类比固体力学的整体刚度/阻尼矩阵）
nDoF=length(H);                                                             
mesh.Flow.Kt=sparse(krow,kcol,kval,nDoF,nDoF);  % 全局渗透刚度（类比整体刚度）
mesh.Flow.Ct=sparse(krow,kcol,cval,nDoF,nDoF);  % 全局容水度（类比整体阻尼）
end