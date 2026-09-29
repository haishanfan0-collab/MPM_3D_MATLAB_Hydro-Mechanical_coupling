function [fint,mpData,Stifiness] = detMPs_BX(uvw,mpData,c_par,Cal_par)

% 【核心-物质点本构与矩阵组装代码】，内置多类型本构模型计算、物质点应变更新、节点内力
% 散度计算、整体刚度/质量/阻尼矩阵拼接功能，该计算代码在源代码基础上重新梳理整理获得
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行物质点力学状态更新，涵盖物质点变形梯度与Hencky对数应变计算、多类型本构模型调用
% （线弹性、Von-Mises理想塑性、Mohr-Coulomb、吸力相关理想弹塑性、大周期弹塑性边界面模型）、
% Cauchy应力更新、节点内力组装、整体刚度/质量/阻尼矩阵稀疏拼接全过程，为物质点法核心子函数
%
%--------------------------------------------------------------------------
% [fint, mpData, Stifiness] = detMPs_BX(uvw, mpData, c_par, Cal_par)
%--------------------------------------------------------------------------
% 输入:
% uvw      - 背景网格节点位移向量 [nDoF×1]
% mpData   - 物质点信息结构体数组
% c_par    - 阻尼系数参数 [c1, c2]，用于组装阻尼矩阵 Ct = c1*Kt + c2*Mt
% Cal_par  - 求解控制参数结构体，含 method 字段（'static'/'dynamic'）
%--------------------------------------------------------------------------
% 输出:
% fint     - 节点内力向量 [nDoF×1]
% mpData   - 更新后的物质点信息结构体数组（含应力、应变、变形梯度等）
% Stifiness- 结构体，包含：
%            Kt - 整体切线刚度矩阵 [nDoF×nDoF]（稀疏）
%            Mt - 整体质量矩阵 [nDoF×nDoF]（稀疏，动力分析时）
%            Ct - 整体阻尼矩阵 [nDoF×nDoF]（稀疏，动力分析时）
%--------------------------------------------------------------------------
% 此子函数包含/调用子函数
% Hooke3d_BX          - 线弹性本构模型
% VMconst_BX          - Von-Mises理想弹塑性本构模型
% MohrCoulomb3d_BX   - Mohr-Coulomb弹塑性本构模型
% VMconst_FLOW_BX    - 吸力相关理想弹塑性本构模型
% Boundary_module_BX - 大周期弹塑性边界面本构模型
% formULstiff_BX     - 空间切线刚度矩阵 A 构造
% batch_det3x3       - 批量3×3矩阵行列式计算
%--------------------------------------------------------------------------
if sum(~isreal(uvw),'all')> 10
    error('dhajsifghlkjashgkja')
else
    uvw = real(uvw);
end

nmp   = length(mpData);                                                     % number of material points
fint  = zeros(size(uvw));                                                   % zero internal force vector

npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]);                                                 % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow; mval=krow; cval=krow;         % zero the stiffness information
ddF   = zeros(nmp,3,3);                                                     % derivative of duvw wrt. spatial position

nD = length(mpData(1).mpC);

if nD==1                                                                    % 1D case
    fPos=1;                                                                 % deformation gradient positions
    aPos=1;                                                                 % material stiffness matrix positions for global stiffness
    sPos=1;                                                                 % Cauchy stress components for internal force
elseif nD==2                                                                % 2D case (plane strain & stress)
    fPos=[1 5 4 2];
    aPos=[1 2 4 5];
    sPos=[1 2 4 4];
else                                                                        % 3D case
    fPos=[1 5 9 4 2 8 6 3 7];
    aPos=[1 2 3 4 5 6 7 8 9];
    sPos=[1 2 3 4 4 5 5 6 6];
end

%% 尝试去除 for 循环 / 参数配置
nIN_all = [mpData.nIN];
nn_all = [mpData.nn]';
ed  = repmat((nIN_all-1)*nD,nD,1)+repmat((1:nD).',1,sum(nn_all));           % degrees of freedom of nodes (matrix form)
ed  = reshape(ed,1,sum(nn_all)*nD);                                         % degrees of freedom of nodes (vector form)
dNx_all = [mpData.dSvp];
nn_sizeall = [mpData.nn_size]';

if nD==1                                                                    % 1D case
    G_all = dNx_all;                                                        % strain-displacement matrix
elseif nD==2                                                                % 2D case (plane strain & stress)
    G_all = zeros(4,nD*sum(nn_all));                                        % zero the strain-disp matrix (2D)
    G_all([1 3],1:nD:end)=dNx_all;                                          % strain-displacement matrix
    G_all([4 2],2:nD:end)=dNx_all;
else                                                                        % 3D case
    G_all=zeros(9,nD*sum(nn_all));                                          % zero the strain-disp matrix (3D)
    G_all([1 4 9],1:nD:end)=dNx_all;                                        % strain-displacement matrix
    G_all([5 2 6],2:nD:end)=dNx_all;
    G_all([8 7 3],3:nD:end)=dNx_all;
end

%% 根据节点位移计算节点应变-一次性
% 1. ddF计算
% 通过索引一次性提取所有需要的 uvw 值
% 结果：values(j) = uvw( ed_all(j) )
values = uvw(ed);  % 189796 x 1
% 逐元素加权（广播机制）
% G_scaled(:,j) = G_all(:,j) * values(j)
G_scaled = G_all .* values';  % 4 x 189796
% 列方向累积求和（前缀和）
Prefix = cumsum(G_scaled, 2);  % 4 x 189796
% 通过索引相减提取各段结果（向量化，无循环）
s = (nn_sizeall(:,1)-1)*nD + 1;
e = nn_sizeall(:,2)*nD;
% 初始化结果矩阵 4 x 23530
Result = Prefix(:, e);
% 处理起始位置大于1的段：减去前一个前缀和
% 利用逻辑掩码处理 s=1 的边界情况
mask = (s > 1);
if any(mask)
    s_prev = s - 1;
    Result(:, mask) = Result(:, mask) - Prefix(:, s_prev(mask));
end

ddF(:,fPos) = Result';

% 2. dF 计算（保持不变）
I = zeros(nmp, 3, 3);
for k = 1:3, I(:,k,k) = 1; end
dF = I + ddF;

% 3. F 计算（若需要，补充缺失的乘法）
Fn_p = reshape([mpData.Fn], 3, 3, nmp);
dF_p = permute(dF, [2, 3, 1]);
F = pagemtimes(dF_p, Fn_p);

% 4. epsEn 处理（保持不变）
epsEn = [mpData.epsEn]';
epsEn(:,4:6) = 0.5*epsEn(:,4:6);
idx = [1 4 6; 4 2 5; 6 5 3];
epsEn = reshape(epsEn(:, idx(:)), nmp, 3, 3);

% 5. 第一次特征值分解（epsEn）
epsEn_paged = permute(epsEn, [2, 3, 1]);
[V, D] = pageeig(epsEn_paged);  % V,D 是 3×3×N

% 6. BeT 计算（使用第一次的 V）
dF_p = permute(dF, [2,3,1]);
Lambda_p = zeros(3,3,nmp);
for k = 1:3, Lambda_p(k,k,:) = exp(2 * D(k,k,:)); end

temp = pagemtimes(V, Lambda_p);
temp = pagemtimes(temp, 'none', V, 'ctranspose');  % 注意：这里 V 是 3×3×N
temp = pagemtimes(dF_p, temp);
BeT = pagemtimes(temp, 'none', dF_p, 'ctranspose');  % BeT 是 3×3×N

% 7. 第二次特征值分解（BeT）
[V, D] = pageeig(BeT);  % V,D 被覆盖为 BeT 的特征系统

% 8. epsEtr 计算（使用第二次的 V，关键修正！）
Lambda_p = zeros(3,3,nmp);
for k = 1:3, Lambda_p(k,k,:) = 0.5 * log(D(k,k,:)); end

epsEtr = pagemtimes(V, Lambda_p);
epsEtr = pagemtimes(epsEtr, 'none', V, 'ctranspose');  % 修正：用 V 而非 V_p
epsEtr = permute(epsEtr, [3, 1, 2]);
epsEtr = [1 1 1 2 2 2] .* epsEtr(:,[1 5 9 2 6 3]);  % N×6

%==========================================================================
% 计算当前总应变 (Hencky 对数应变)
%==========================================================================
% 1. C = F' * F
C = pagemtimes(F, 'ctranspose', F, 'none');  % 3×3×N[V, D] = eig(C);
% 2. 特征值分解 (批量)
[V, D] = pageeig(C);  % V: 3×3×N (特征向量), D: 3×3×N (对角特征值)
% 3. 提取特征值并处理数值稳定性
% 提取对角元 (主拉伸的平方 λ²)，确保为正避免 log(0)
lambda_sq = zeros(nmp, 3);
for k = 1:3
    lambda_sq(:,k) = squeeze(D(k,k,:));  % 明确 squeeze 确保维度为 nmp×1; 
end
lambda_sq = max(real(lambda_sq), 1e-20);  % N×3
% 4. 主对数应变 (批量)
eps_principal = 0.5 * log(lambda_sq);  % N×3
% 5. 谱重构回张量形式 (批量)
% 构造批量对角矩阵 Lambda (3×3×N)
Lambda_p = zeros(3,3,nmp);
for k = 1:3
    Lambda_p(k,k,:) = eps_principal(:,k);
end
% eps_tensor = V * Lambda * V'
eps_tensor = pagemtimes(V, Lambda_p);
eps_tensor = pagemtimes(eps_tensor, 'none', V, 'ctranspose');  % 3×3×N
% 6. 转换为 Voigt 格式 (N×6)
% 将 3×3×N reshape 为 N×9 (列优先: 11,21,31,12,22,32,13,23,33)
eps_tensor = permute(eps_tensor, [3, 1, 2]);
eps_flat = reshape(eps_tensor, nmp, 9);
% 提取 [11,22,33,12,13,23] 对应列索引 1,5,9,4,7,8
% 注意：剪应变位置 12,13,23 对应张量分量，后续乘2转为工程剪应变
eps_total = eps_flat(:, [1 5 9 2 6 3]);  % N×6
% 7. 张量剪应变 → 工程剪应变 (批量)
eps_total(:, 4:6) = 2 * eps_total(:, 4:6);

% 组装材料矩阵
FHS_cmType = [mpData.cmType];
mCst_1 = [mpData(FHS_cmType==1).mCst]';  mCst_1 = reshape(mCst_1,2,length(mCst_1)/2); mCst_1=mCst_1';
mCst_2 = [mpData(FHS_cmType==2).mCst]';  mCst_2 = reshape(mCst_2,3,length(mCst_2)/3); mCst_2=mCst_2';
mCst_3 = [mpData(FHS_cmType==3).mCst]';  mCst_3 = reshape(mCst_3,6,length(mCst_3)/6); mCst_3=mCst_3';
mCst_4 = [mpData(FHS_cmType==4).mCst]';  mCst_4 = reshape(mCst_4,4,length(mCst_4)/4); mCst_4=mCst_4';
mCst_5 = [mpData(FHS_cmType==5).mCst]';  mCst_5 = reshape(mCst_5,32,length(mCst_5)/32); mCst_5=mCst_5';

%----------------------------------------------------------------------     % Constitutive model
D = zeros(6,6,nmp);
Ksig = zeros(nmp,6);
epsE = zeros(nmp,6);
epsP = zeros(nmp,6);
sign = [mpData.sign]';
epsn = [mpData.epsn]';
epsEn = [mpData.epsEn]';
STATEV = [mpData.STATEV]';
STATEVn = [mpData.STATEVn]';
e0 = [mpData.e0]';
Flow_CP = reshape([mpData.Flow_CP]',3,nmp)';
Mc = [mpData.Gravity_Water_Content_old]';
dMc = [mpData.Gravity_Water_Content]' - [mpData.Gravity_Water_Content_old]';
Pore_Pressure = [mpData.Pore_Pressure]'*1000;
Out = []; 

%% 本构1 线弹性本构
if ~isempty(mCst_1)
    [D0,Ksig0,epsE0,epsP0,STATEV0]=Hooke3d_BX(epsEtr(FHS_cmType==1,:),mCst_1,STATEV(FHS_cmType==1,:),Cal_par);
    D(:,:,FHS_cmType==1) = D0;
    Ksig(FHS_cmType==1,:) = Ksig0;
    epsE(FHS_cmType==1,:) = epsE0;
    epsP(FHS_cmType==1,:) = epsP0;
    STATEV(FHS_cmType==1,:) = STATEV0;
end

%% 本构2 理想塑性本构
if ~isempty(mCst_2)
    [D0,Ksig0,epsE0,epsP0]=VMconst_BX(epsEtr(FHS_cmType==2,:),mCst_2);
    D(:,:,FHS_cmType==2) = D0;
    Ksig(FHS_cmType==2,:) = Ksig0;
    epsE(FHS_cmType==2,:) = epsE0;
    epsP(FHS_cmType==2,:) = epsP0;
end

%% 本构3 摩尔库伦本构
if ~isempty(mCst_3)
    modulus = [mpData.modulus]';
    [D0,Ksig0,epsE0,epsP0,Out]=MohrCoulomb3d_BX(eps_total(FHS_cmType==3,:), ...
        mCst_3, STATEV(FHS_cmType==3,:),STATEVn(FHS_cmType==3,:), ...
        Cal_par, epsn(FHS_cmType==3,:), ...
        sign(FHS_cmType==3,:),epsEn(FHS_cmType==3,:),e0(FHS_cmType==3,:),Flow_CP(FHS_cmType==3,:),...
        Mc(FHS_cmType==3,:),dMc(FHS_cmType==3,:),Pore_Pressure(FHS_cmType==3,:));

    STATEV(FHS_cmType==3,:) = Out;
    D(:,:,FHS_cmType==3) = D0;
    Ksig(FHS_cmType==3,:) = Ksig0;
    epsE(FHS_cmType==3,:) = epsE0;
    epsP(FHS_cmType==3,:) = epsP0;
end

%% 本构4 待定本构-现初步设置：屈服应力与吸力相关的理想弹塑性本构
if ~isempty(mCst_4)
    Suction = [mpData.Scution_Flow]';Suction = Suction(FHS_cmType==4,:);
    [D0,Ksig0,epsE0,epsP0]=VMconst_FLOW_BX(epsEtr(FHS_cmType==4,:),mCst_4,Suction); 
    D(:,:,FHS_cmType==4) = D0;
    Ksig(FHS_cmType==4,:) = Ksig0;
    epsE(FHS_cmType==4,:) = epsE0;
    epsP(FHS_cmType==4,:) = epsP0;
end

%% 本构5 大周期弹塑性边界面模型
if ~isempty(mCst_5)
    modulus = [mpData.modulus]';

    [D0,Ksig0,epsE0,epsP0,Out]=Boundary_module_BX(eps_total(FHS_cmType==5,:),...
                               mCst_5,Cal_par,sign(FHS_cmType==5,:),epsn(FHS_cmType==5,:),...
                               epsEn(FHS_cmType==5,:),STATEVn(FHS_cmType==5,:),...
                               Mc(FHS_cmType==5,:),dMc(FHS_cmType==5,:),STATEV(FHS_cmType==5,:),Pore_Pressure(FHS_cmType==5,:));          % elastic behaviour

    D(:,:,FHS_cmType==5) = D0;
    Ksig(FHS_cmType==5,:) = Ksig0;
    try
        epsE(FHS_cmType==5,:) = epsE0;
        epsP(FHS_cmType==5,:) = epsP0;
    catch
    end

    if ~isempty(Out)
        modulus(FHS_cmType==5,:) = Out{1};
        STATEV(FHS_cmType==5,:) = Out{2};
    end

end

%% A矩阵构建，
det_F = batch_det3x3(F);
sig = Ksig./det_F;  
A   = formULstiff_BX(F,D,sig,BeT);                                          % spatial tangent stiffness matrix

iF = pagemldivide(dF_p, eye(3));  % 或者 iF = pageinv(dF);
dXdx = zeros(9, 9, nmp);
% 提取各分量变为 1×1×N，便于广播到 9×9×N 的对应切片
F11 = iF(1,1,:); F21 = iF(2,1,:); F31 = iF(3,1,:);
F12 = iF(1,2,:); F22 = iF(2,2,:); F32 = iF(3,2,:);
F13 = iF(1,3,:); F23 = iF(2,3,:); F33 = iF(3,3,:);

% 第1行: [F11, 0, 0, F21, 0, 0, 0, 0, F31]
dXdx(1,1,:) = F11; dXdx(1,4,:) = F21; dXdx(1,9,:) = F31;
% 第2行: [0, F22, 0, 0, F12, F32, 0, 0, 0]
dXdx(2,2,:) = F22; dXdx(2,5,:) = F12; dXdx(2,6,:) = F32;
% 第3行: [0, 0, F33, 0, 0, 0, F23, F13, 0]
dXdx(3,3,:) = F33; dXdx(3,7,:) = F23; dXdx(3,8,:) = F13;
% 第4行: [F12, 0, 0, F22, 0, 0, 0, 0, F32]
dXdx(4,1,:) = F12; dXdx(4,4,:) = F22; dXdx(4,9,:) = F32;
% 第5行: [0, F21, 0, 0, F11, F31, 0, 0, 0]
dXdx(5,2,:) = F21; dXdx(5,5,:) = F11; dXdx(5,6,:) = F31;
% 第6行: [0, F23, 0, 0, F13, F33, 0, 0, 0]
dXdx(6,2,:) = F23; dXdx(6,5,:) = F13; dXdx(6,6,:) = F33;
% 第7行: [0, 0, F32, 0, 0, 0, F22, F12, 0]
dXdx(7,3,:) = F32; dXdx(7,7,:) = F22; dXdx(7,8,:) = F12;
% 第8行: [0, 0, F31, 0, 0, 0, F21, F11, 0]
dXdx(8,3,:) = F31; dXdx(8,7,:) = F21; dXdx(8,8,:) = F11;
% 第9行: [F13, 0, 0, F23, 0, 0, 0, 0, F33]
dXdx(9,1,:) = F13; dXdx(9,4,:) = F23; dXdx(9,9,:) = F33;

dXdx = dXdx(aPos,aPos,:);

dets = batch_det3x3(dF_p);
A = A(aPos,aPos,:);
vp = [mpData.vp]';

vp_dets = vp.*dets;
vp_dets = reshape(vp_dets,1,1,nmp);

%% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 4
if nD == 2
    max_4_YES = (nn_all == 4);
    sig_3d = reshape(sig(:,sPos)',4,1,nmp);
    if sum(max_4_YES)>=4000
        cols_3d_4 = [s(max_4_YES),s(max_4_YES)+1,s(max_4_YES)+2,s(max_4_YES)+3,...
            s(max_4_YES)+4,s(max_4_YES)+5,s(max_4_YES)+6,s(max_4_YES)+7];
        G_3d_4 = reshape(G_all(:, cols_3d_4'), 4, 8, sum(max_4_YES));
        G_3d_4 = pagemtimes(dXdx(:,:,max_4_YES), G_3d_4);

        % 节点刚度矩阵
        kp_3d_ls = pagemtimes(G_3d_4,'ctranspose',A(:,:,max_4_YES),'none');
        kp_3d_4 = pagemtimes(kp_3d_ls,G_3d_4);
        kp_3d_4 = pagemtimes(vp_dets(:,:,max_4_YES),kp_3d_4);

        fp_3d_4 = squeeze(pagemtimes(G_3d_4,'ctranspose',sig_3d(:,:,max_4_YES),'none'));
        fp_3d_4 = reshape(vp_dets(max_4_YES),1,sum(max_4_YES)) .* fp_3d_4;
    else
        max_4_YES = (nn_all == -1);
    end

    %% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 6
    max_6_YES = (nn_all == 6);
    if sum(max_6_YES)>=2000
        cols_3d_6 = [s(max_6_YES),s(max_6_YES)+1,s(max_6_YES)+2,s(max_6_YES)+3,...
            s(max_6_YES)+4,s(max_6_YES)+5,s(max_6_YES)+6,s(max_6_YES)+7,...
            s(max_6_YES)+8,s(max_6_YES)+9,s(max_6_YES)+10,s(max_6_YES)+11];
        G_3d_6 = reshape(G_all(:, cols_3d_6'), 4, 12, sum(max_6_YES));
        G_3d_6 = pagemtimes(dXdx(:,:,max_6_YES), G_3d_6);

        % 节点刚度矩阵
        kp_3d_ls = pagemtimes(G_3d_6,'ctranspose',A(:,:,max_6_YES),'none');
        kp_3d_6 = pagemtimes(kp_3d_ls,G_3d_6);
        kp_3d_6 = pagemtimes(vp_dets(:,:,max_6_YES),kp_3d_6);

        fp_3d_6 = squeeze(pagemtimes(G_3d_6,'ctranspose',sig_3d(:,:,max_6_YES),'none'));
        fp_3d_6 = reshape(vp_dets(max_6_YES),1,sum(max_6_YES)) .* fp_3d_6;
    else
        max_6_YES = (nn_all == -1);
    end

    %% 尝试去除for循环，仅赋值(对于尺寸是4的直接矩阵运算，其他部分不矩阵运算）- 9
    max_9_YES = (nn_all == 9);
    if sum(max_9_YES)>=1000
        cols_3d_9 = [s(max_9_YES),s(max_9_YES)+1,s(max_9_YES)+2,s(max_9_YES)+3,...
            s(max_9_YES)+4,s(max_9_YES)+5,s(max_9_YES)+6,s(max_9_YES)+7,...
            s(max_9_YES)+8,s(max_9_YES)+9,s(max_9_YES)+10,s(max_9_YES)+11,...
            s(max_9_YES)+12,s(max_9_YES)+13,s(max_9_YES)+14,s(max_9_YES)+15,...
            s(max_9_YES)+16,s(max_9_YES)+17];
        G_3d_9 = reshape(G_all(:, cols_3d_9'), 4, 18, sum(max_9_YES));
        G_3d_9 = pagemtimes(dXdx(:,:,max_9_YES), G_3d_9);

        % 节点刚度矩阵
        kp_3d_ls = pagemtimes(G_3d_9,'ctranspose',A(:,:,max_9_YES),'none');
        kp_3d_9 = pagemtimes(kp_3d_ls,G_3d_9);
        kp_3d_9 = pagemtimes(vp_dets(:,:,max_9_YES),kp_3d_9);

        fp_3d_9 = squeeze(pagemtimes(G_3d_9,'ctranspose',sig_3d(:,:,max_9_YES),'none'));
        fp_3d_9 = reshape(vp_dets(max_9_YES),1,sum(max_9_YES)) .* fp_3d_9;
    else
        max_9_YES = (nn_all == -1);
    end
    % disp(['4-',num2str(sum((nn_all == 4))), ' 5-', num2str(sum((nn_all == 5))), ...
    %      ' 6-',num2str(sum((nn_all == 6))), ' 7-', num2str(sum((nn_all == 7))), ...
    %      ' 8-',num2str(sum((nn_all == 8))), ' 9-', num2str(sum((nn_all == 9))),])
else
    sig_3d = reshape(sig(:,sPos)',9,1,nmp);
    max_4_YES = (nn_all == 8);
    if sum(max_4_YES)>=4000
        cols_3d_4 = [s(max_4_YES),s(max_4_YES)+1,s(max_4_YES)+2,s(max_4_YES)+3,...
            s(max_4_YES)+4,s(max_4_YES)+5,s(max_4_YES)+6,s(max_4_YES)+7,...
            s(max_4_YES)+8,s(max_4_YES)+9,s(max_4_YES)+10,s(max_4_YES)+11,...
            s(max_4_YES)+12,s(max_4_YES)+13,s(max_4_YES)+14,s(max_4_YES)+15,...
            s(max_4_YES)+16,s(max_4_YES)+17,s(max_4_YES)+18,s(max_4_YES)+19,...
            s(max_4_YES)+20,s(max_4_YES)+21,s(max_4_YES)+22,s(max_4_YES)+23];
        G_3d_4 = reshape(G_all(:, cols_3d_4'), 9, 24, sum(max_4_YES));
        G_3d_4 = pagemtimes(dXdx(:,:,max_4_YES), G_3d_4);

        % 节点刚度矩阵
        kp_3d_ls = pagemtimes(G_3d_4,'ctranspose',A(:,:,max_4_YES),'none');
        kp_3d_4 = pagemtimes(kp_3d_ls,G_3d_4);
        kp_3d_4 = pagemtimes(vp_dets(:,:,max_4_YES),kp_3d_4);

        fp_3d_4 = squeeze(pagemtimes(G_3d_4,'ctranspose',sig_3d(:,:,max_4_YES),'none'));
        fp_3d_4 = reshape(vp_dets(max_4_YES),1,sum(max_4_YES)) .* fp_3d_4;
    else
        max_4_YES = (nn_all == -1);
    end
    %% 12
    max_6_YES = (nn_all == 12);
    if sum(max_6_YES)>=2000
        cols_3d_6 = [s(max_6_YES),s(max_6_YES)+1,s(max_6_YES)+2,s(max_6_YES)+3,...
            s(max_6_YES)+4,s(max_6_YES)+5,s(max_6_YES)+6,s(max_6_YES)+7,...
            s(max_6_YES)+8,s(max_6_YES)+9,s(max_6_YES)+10,s(max_6_YES)+11,...
            s(max_6_YES)+12,s(max_6_YES)+13,s(max_6_YES)+14,s(max_6_YES)+15,...
            s(max_6_YES)+16,s(max_6_YES)+17,s(max_6_YES)+18,s(max_6_YES)+19,...
            s(max_6_YES)+20,s(max_6_YES)+21,s(max_6_YES)+22,s(max_6_YES)+23,...
            s(max_6_YES)+24,s(max_6_YES)+25,s(max_6_YES)+26,s(max_6_YES)+27,...
            s(max_6_YES)+28,s(max_6_YES)+29,s(max_6_YES)+30,s(max_6_YES)+31,...
            s(max_6_YES)+32,s(max_6_YES)+33,s(max_6_YES)+34,s(max_6_YES)+35];
        G_3d_6 = reshape(G_all(:, cols_3d_6'), 9, 36, sum(max_6_YES));
        G_3d_6 = pagemtimes(dXdx(:,:,max_6_YES), G_3d_6);

        % 节点刚度矩阵
        kp_3d_ls = pagemtimes(G_3d_6,'ctranspose',A(:,:,max_6_YES),'none');
        kp_3d_6 = pagemtimes(kp_3d_ls,G_3d_6);
        kp_3d_6 = pagemtimes(vp_dets(:,:,max_6_YES),kp_3d_6);

        fp_3d_6 = squeeze(pagemtimes(G_3d_6,'ctranspose',sig_3d(:,:,max_6_YES),'none'));
        fp_3d_6 = reshape(vp_dets(max_6_YES),1,sum(max_6_YES)) .* fp_3d_6;
    else
        max_6_YES = (nn_all == -1);
    end
    %% 27
    max_9_YES = (nn_all == 27);
    if sum(max_9_YES)>=1000
        cols_3d_9 = [s(max_9_YES),s(max_9_YES)+1,s(max_9_YES)+2,s(max_9_YES)+3,...
            s(max_9_YES)+4,s(max_9_YES)+5,s(max_9_YES)+6,s(max_9_YES)+7,...
            s(max_9_YES)+8,s(max_9_YES)+9,s(max_9_YES)+10,s(max_9_YES)+11,...
            s(max_9_YES)+12,s(max_9_YES)+13,s(max_9_YES)+14,s(max_9_YES)+15,...
            s(max_9_YES)+16,s(max_9_YES)+17,s(max_9_YES)+18,s(max_9_YES)+19,...
            s(max_9_YES)+20,s(max_9_YES)+21,s(max_9_YES)+22,s(max_9_YES)+23,...
            s(max_9_YES)+24,s(max_9_YES)+25,s(max_9_YES)+26,s(max_9_YES)+27,...
            s(max_9_YES)+28,s(max_9_YES)+29,s(max_9_YES)+30,s(max_9_YES)+31,...
            s(max_9_YES)+32,s(max_9_YES)+33,s(max_9_YES)+34,s(max_9_YES)+35,...
            s(max_9_YES)+36,s(max_9_YES)+37,s(max_9_YES)+38,s(max_9_YES)+39,...
            s(max_9_YES)+40,s(max_9_YES)+41,s(max_9_YES)+42,s(max_9_YES)+43,...
            s(max_9_YES)+44,s(max_9_YES)+45,s(max_9_YES)+46,s(max_9_YES)+47,...
            s(max_9_YES)+48,s(max_9_YES)+49,s(max_9_YES)+50,s(max_9_YES)+51,...
            s(max_9_YES)+52,s(max_9_YES)+53,s(max_9_YES)+54,s(max_9_YES)+55,...
            s(max_9_YES)+56,s(max_9_YES)+57,s(max_9_YES)+58,s(max_9_YES)+59,...
            s(max_9_YES)+60,s(max_9_YES)+61,s(max_9_YES)+62,s(max_9_YES)+63,...
            s(max_9_YES)+64,s(max_9_YES)+65,s(max_9_YES)+66,s(max_9_YES)+67,...
            s(max_9_YES)+68,s(max_9_YES)+69,s(max_9_YES)+70,s(max_9_YES)+71,...
            s(max_9_YES)+72,s(max_9_YES)+73,s(max_9_YES)+74,s(max_9_YES)+75,...
            s(max_9_YES)+76,s(max_9_YES)+77,s(max_9_YES)+78,s(max_9_YES)+79,...
            s(max_9_YES)+80];
        G_3d_9 = reshape(G_all(:, cols_3d_9'), 9, 81, sum(max_9_YES));
        G_3d_9 = pagemtimes(dXdx(:,:,max_9_YES), G_3d_9);

        % 节点刚度矩阵
        kp_3d_ls = pagemtimes(G_3d_9,'ctranspose',A(:,:,max_9_YES),'none');
        kp_3d_9 = pagemtimes(kp_3d_ls,G_3d_9);
        kp_3d_9 = pagemtimes(vp_dets(:,:,max_9_YES),kp_3d_9);

        fp_3d_9 = squeeze(pagemtimes(G_3d_9,'ctranspose',sig_3d(:,:,max_9_YES),'none'));
        fp_3d_9 = reshape(vp_dets(max_9_YES),1,sum(max_9_YES)) .* fp_3d_9;
    else
        max_9_YES = (nn_all == -1);
    end
end
%% 循环赋值，仅内存拷贝
k_4_JS=1;  k_6_JS=1; k_9_JS=1;
for mp=1:nmp                                                                % material point loop
    if max_4_YES(mp)
        kp = kp_3d_4(:,:,k_4_JS); fp = fp_3d_4(:,k_4_JS); k_4_JS = k_4_JS+1;
    elseif max_6_YES(mp)
        kp = kp_3d_6(:,:,k_6_JS); fp = fp_3d_6(:,k_6_JS); k_6_JS = k_6_JS+1;
    elseif max_9_YES(mp)
        kp = kp_3d_9(:,:,k_9_JS); fp = fp_3d_9(:,k_9_JS); k_9_JS = k_9_JS+1;
    else
        %----------------------------------------------------------------------
        cols = s(mp):e(mp);                                                 % 当前块列范围
        % 左乘：dXdx(:,:,i) * G_all(:,cols)
        G = dXdx(:,:,mp) * G_all(:, cols);
        % 节点刚度矩阵
        kp = vp_dets(mp)*(G.'*A(:,:,mp)*G);                                 % material point stiffness contribution
        % 节点应力矩阵
        fp = vp_dets(mp)*(G.'*sig(mp,sPos)');                               % internal force contribution
    end

    mpData(mp).F    = F(:,:,mp);                                            % store deformation gradient
    mpData(mp).sig  = Ksig(mp,:)';                                          % store Cauchy stress
    mpData(mp).epsE = epsE(mp,:)';                                          % store elastic logarithmic strain
    mpData(mp).epsP = epsP(mp,:)';
    mpData(mp).eps = eps_total(mp,:)';

    if ~isempty(Out)
        mpData(mp).modulus = modulus(mp);
        mpData(mp).STATEV = STATEV(mp,:)';
    end

    ed0=ed(s(mp):e(mp));
    npDoF=(nn_all(mp)*nD)^2;                                                % no. entries in kp
    nnDoF=nn_all(mp)*nD;                                                    % no. DoF in kp
    krow(npCnt+1:npCnt+npDoF)=repmat(ed0.',nnDoF,1);                        % row position storage
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed0  ,nnDoF,1);                        % column position storage

    if strcmp(Cal_par.method, 'static')
        % 矩阵组装
        kval(npCnt+1:npCnt+npDoF)=kp;                                       % stiffness storage
        fint(ed0)=fint(ed0)+fp;                                             % internal force contribution

    elseif strcmp(Cal_par.method, 'dynamic')
        % 节点质量矩阵
        m_Matrix = mpData(mp).mpM * mpData(mp).Svp / sum(mpData(mp).Svp) / nD;
        m_Matrix = diag(repelem(m_Matrix, nD));
        % 矩阵组装
        kval(npCnt+1:npCnt+npDoF)=kp;                                       % stiffness storage
        mval(npCnt+1:npCnt+npDoF)=m_Matrix;                                 % mass storage
        cval(npCnt+1:npCnt+npDoF)=c_par(1)*kp + c_par(2)*m_Matrix;          % damp storage
        fint(ed0)=fint(ed0)+fp;                                             % internal force contribution
    else
        error('method only can be ''dynamic'' and ''static'' !!!')
    end
    npCnt=npCnt+npDoF;                                                      % number of entries in Kt
end

nDoF=length(uvw);                                                           % number of degrees of freedom
Stifiness.Kt=sparse(krow,kcol,kval,nDoF,nDoF);                              % form the global stiffness matrix
Stifiness.Mt=sparse(krow,kcol,mval,nDoF,nDoF);                              % form the global maxx matrix
Stifiness.Ct=sparse(krow,kcol,cval,nDoF,nDoF);                              % form the global damp matrix
end

function d = batch_det3x3(F)
% BATCH_DET3X3 批量计算 3×3 矩阵行列式（萨吕法则向量化）
%   F: 3×3×N 数组
%   d: N×1 列向量（与 arrayfun(@det) 输出维度一致）
N = size(F, 3);

% 一次性 reshape 为 9×N，列优先顺序：
% [F11; F21; F31; F12; F22; F32; F13; F23; F33]
f = reshape(F, 9, N);

% 萨吕法则向量化计算（6 项展开）
d = (f(1,:).*f(5,:).*f(9,:) ...  % F11*F22*F33
   + f(4,:).*f(8,:).*f(3,:) ...  % F12*F23*F31
   + f(7,:).*f(2,:).*f(6,:) ...  % F13*F21*F32
   - f(7,:).*f(5,:).*f(3,:) ...  % F13*F22*F31
   - f(1,:).*f(8,:).*f(6,:) ...  % F11*F23*F32
   - f(4,:).*f(2,:).*f(9,:));    % F12*F21*F33

% 确保输出为列向量（与 arrayfun 行为一致）
d = d(:);
end