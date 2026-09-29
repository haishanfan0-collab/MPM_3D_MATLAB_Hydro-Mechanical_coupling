function [fint, mpData, Stifiness] = detMPs_BX_explicit(V_old, dt, mpData, c_par, Cal_par)

% 【核心-物质点显式动力本构与矩阵组装代码】
% 基于上一步速度 V_old 直接计算速度梯度，积分得到应变增量，调用本构模型更新应力，
% 组装节点内力、整体质量矩阵及纯质量比例阻尼矩阵。不含刚度矩阵。
% 应变更新与位移更新解耦：应变基于 V_old，位移由主程序通过 ExplicitSolve 更新。
%--------------------------------------------------------------------------
% 作者: FAN Haishan (Modified for Explicit USF Format)
% 日期: 07/05/2026
% 描述:
% 显式USF格式物质点更新：V_old → 速度梯度 → 增量位移梯度 ddF = (∂v/∂X)*dt → 
% 变形梯度更新 F = (I+ddF)*Fn → Hencky应变 → 本构 → 内力 fint。
% 节点位移增量 dU 由主程序通过 ExplicitSolve 独立计算，不传入本子函数。
%
%--------------------------------------------------------------------------
% [fint, mpData, Stifiness] = detMPs_BX_Explicit(V_old, dt, mpData, c_par, Cal_par)
%--------------------------------------------------------------------------
% 输入:
% V_old    - 背景网格节点速度向量 [nDoF×1]（上一步速度 v_{n-1/2}，已知量）
% dt       - 时间步长
% mpData   - 物质点信息结构体数组
% c_par    - 阻尼系数参数 [c1, c2]，显式中仅 c2 有效（Ct = c2*Mt）
% Cal_par  - 求解控制参数结构体
%--------------------------------------------------------------------------
% 输出:
% fint     - 节点内力向量 [nDoF×1]
% mpData   - 更新后的物质点信息结构体数组
% Stifiness- 结构体：
%            Kt - 空稀疏矩阵 [nDoF×nDoF]（显式不使用，保持兼容性）
%            Mt - 整体质量矩阵 [nDoF×nDoF]（稀疏）
%            Ct - 整体阻尼矩阵 [nDoF×nDoF]（稀疏，Ct = c_par(2)*Mt）

nmp   = length(mpData);
fint  = zeros(size(V_old));

npCnt = 0;
tnSMe = sum([mpData.nSMe]);
krow  = zeros(tnSMe,1); kcol=krow; mval=krow; cval=krow;

nD = length(mpData(1).mpC);

if nD==1
    fPos=1; 
    aPos=1; 
    sPos=1;
elseif nD==2
    fPos=[1 5 4 2]; 
    aPos=[1 2 4 5]; 
    sPos=[1 2 4 4];
else
    fPos=[1 5 9 4 2 8 6 3 7]; 
    aPos=[1 2 3 4 5 6 7 8 9]; 
    sPos=[1 2 3 4 4 5 5 6 6];
end

%% 参数配置（与隐式版本一致）
nIN_all = [mpData.nIN];
nn_all = [mpData.nn]';
ed  = repmat((nIN_all-1)*nD,nD,1)+repmat((1:nD).',1,sum(nn_all));
ed  = reshape(ed,1,sum(nn_all)*nD);
dNx_all = [mpData.dSvp];
nn_sizeall = [mpData.nn_size]';

if nD==1
    G_all = dNx_all;
elseif nD==2
    G_all = zeros(4,nD*sum(nn_all));
    G_all([1 3],1:nD:end)=dNx_all;
    G_all([4 2],2:nD:end)=dNx_all;
else
    G_all=zeros(9,nD*sum(nn_all));
    G_all([1 4 9],1:nD:end)=dNx_all;
    G_all([5 2 6],2:nD:end)=dNx_all;
    G_all([8 7 3],3:nD:end)=dNx_all;
end

%% ================== 核心修改：输入 V_old 而非位移增量 ==================
% 显式USF格式：V_old 是已知量，直接计算速度梯度，乘以 dt 得增量位移梯度
% 这样应变更新基于 V_old，与主程序中 ExplicitSolve 求解的 dU 完全解耦

values = V_old(ed);                         % 提取节点速度
G_scaled = G_all .* values';                % G * V_old = ∂v/∂X（参考构型速度梯度）
Prefix = cumsum(G_scaled, 2);
s = (nn_sizeall(:,1)-1)*nD + 1;
e = nn_sizeall(:,2)*nD;
Result = Prefix(:, e);
mask = (s > 1);
if any(mask)
    s_prev = s - 1;
    Result(:, mask) = Result(:, mask) - Prefix(:, s_prev(mask));
end

% 速度梯度（参考构型）∂v/∂X
ddF_v = zeros(nmp, 3, 3);
ddF_v(:,fPos) = Result';

% 增量位移梯度 = 速度梯度 * dt
ddF = ddF_v * dt;

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
C = pagemtimes(F, 'ctranspose', F, 'none');
[V_eig, D] = pageeig(C);
lambda_sq = zeros(nmp, 3);
for k = 1:3
    lambda_sq(:,k) = squeeze(D(k,k,:));
end
lambda_sq = max(real(lambda_sq), 1e-20);
eps_principal = 0.5 * log(lambda_sq);
Lambda_p = zeros(3,3,nmp);
for k = 1:3
    Lambda_p(k,k,:) = eps_principal(:,k);
end
eps_tensor = pagemtimes(V_eig, Lambda_p);
eps_tensor = pagemtimes(eps_tensor, 'none', V_eig, 'ctranspose');
eps_tensor = permute(eps_tensor, [3, 1, 2]);
eps_flat = reshape(eps_tensor, nmp, 9);
eps_total = eps_flat(:, [1 5 9 2 6 3]);
eps_total(:, 4:6) = 2 * eps_total(:, 4:6);

% 组装材料矩阵
FHS_cmType = [mpData.cmType];
mCst_1 = [mpData(FHS_cmType==1).mCst]';  mCst_1 = reshape(mCst_1,2,length(mCst_1)/2); mCst_1=mCst_1';
mCst_2 = [mpData(FHS_cmType==2).mCst]';  mCst_2 = reshape(mCst_2,3,length(mCst_2)/3); mCst_2=mCst_2';
mCst_3 = [mpData(FHS_cmType==3).mCst]';  mCst_3 = reshape(mCst_3,6,length(mCst_3)/6); mCst_3=mCst_3';
mCst_4 = [mpData(FHS_cmType==4).mCst]';  mCst_4 = reshape(mCst_4,4,length(mCst_4)/4); mCst_4=mCst_4';
mCst_5 = [mpData(FHS_cmType==5).mCst]';  mCst_5 = reshape(mCst_5,32,length(mCst_5)/32); mCst_5=mCst_5';

% 本构调用（不接收 D0）
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
    [~,Ksig0,epsE0,epsP0,STATEV0]=Hooke3d_BX(epsEtr(FHS_cmType==1,:),mCst_1,STATEV(FHS_cmType==1,:),Cal_par);
    Ksig(FHS_cmType==1,:) = Ksig0;
    epsE(FHS_cmType==1,:) = epsE0;
    epsP(FHS_cmType==1,:) = epsP0;
    STATEV(FHS_cmType==1,:) = STATEV0;
end

%% 本构2 理想塑性本构
if ~isempty(mCst_2)
    [~,Ksig0,epsE0,epsP0]=VMconst_BX(epsEtr(FHS_cmType==2,:),mCst_2);
    Ksig(FHS_cmType==2,:) = Ksig0;
    epsE(FHS_cmType==2,:) = epsE0;
    epsP(FHS_cmType==2,:) = epsP0;
end

%% 本构3 摩尔库伦本构
if ~isempty(mCst_3)
    [~,Ksig0,epsE0,epsP0,STATEV0]=MohrCoulomb3d_BX(eps_total(FHS_cmType==3,:), ...
        mCst_3, STATEV(FHS_cmType==3,:),STATEVn(FHS_cmType==3,:), ...
        Cal_par, epsn(FHS_cmType==3,:), ...
        sign(FHS_cmType==3,:),epsEn(FHS_cmType==3,:),e0(FHS_cmType==3,:),Flow_CP(FHS_cmType==3,:),...
        Mc(FHS_cmType==3,:),dMc(FHS_cmType==3,:),Pore_Pressure(FHS_cmType==3,:));
    STATEV(FHS_cmType==3,:) = STATEV0;
    Ksig(FHS_cmType==3,:) = Ksig0;
    epsE(FHS_cmType==3,:) = epsE0;
    epsP(FHS_cmType==3,:) = epsP0;
end

%% 本构4 吸力相关理想弹塑性本构
if ~isempty(mCst_4)
    Suction = [mpData.Scution_Flow]'; Suction = Suction(FHS_cmType==4,:);
    [~,Ksig0,epsE0,epsP0]=VMconst_FLOW_BX(epsEtr(FHS_cmType==4,:),mCst_4,Suction); 
    Ksig(FHS_cmType==4,:) = Ksig0;
    epsE(FHS_cmType==4,:) = epsE0;
    epsP(FHS_cmType==4,:) = epsP0;
end

%% 本构5 大周期弹塑性边界面模型
if ~isempty(mCst_5)
    modulus = [mpData.modulus]';
    [~,Ksig0,epsE0,epsP0,Out]=Boundary_module_BX(eps_total(FHS_cmType==5,:),...
                               mCst_5,Cal_par,sign(FHS_cmType==5,:),epsn(FHS_cmType==5,:),...
                               epsEn(FHS_cmType==5,:),STATEVn(FHS_cmType==5,:),...
                               Mc(FHS_cmType==5,:),dMc(FHS_cmType==5,:),STATEV(FHS_cmType==5,:),Pore_Pressure(FHS_cmType==5,:));
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

%% 内力计算（显式中不计算刚度矩阵，仅保留内力组装）
det_F = batch_det3x3(F);
sig = Ksig ./ det_F;  % Cauchy应力

iF = pagemldivide(dF_p, eye(3));
dXdx = zeros(9, 9, nmp);
F11 = iF(1,1,:); F21 = iF(2,1,:); F31 = iF(3,1,:);
F12 = iF(1,2,:); F22 = iF(2,2,:); F32 = iF(3,2,:);
F13 = iF(1,3,:); F23 = iF(2,3,:); F33 = iF(3,3,:);

dXdx(1,1,:) = F11; dXdx(1,4,:) = F21; dXdx(1,9,:) = F31;
dXdx(2,2,:) = F22; dXdx(2,5,:) = F12; dXdx(2,6,:) = F32;
dXdx(3,3,:) = F33; dXdx(3,7,:) = F23; dXdx(3,8,:) = F13;
dXdx(4,1,:) = F12; dXdx(4,4,:) = F22; dXdx(4,9,:) = F32;
dXdx(5,2,:) = F21; dXdx(5,5,:) = F11; dXdx(5,6,:) = F31;
dXdx(6,2,:) = F23; dXdx(6,5,:) = F13; dXdx(6,6,:) = F33;
dXdx(7,3,:) = F32; dXdx(7,7,:) = F22; dXdx(7,8,:) = F12;
dXdx(8,3,:) = F31; dXdx(8,7,:) = F21; dXdx(8,8,:) = F11;
dXdx(9,1,:) = F13; dXdx(9,4,:) = F23; dXdx(9,9,:) = F33;

dXdx = dXdx(aPos,aPos,:);
dets = batch_det3x3(dF_p);
vp = [mpData.vp]';
vp_dets = vp .* dets;
vp_dets = reshape(vp_dets,1,1,nmp);

%% 批量内力计算（删除刚度矩阵 kp 计算）
max_4_YES = (nn_all == 4);
if nD == 2
sig_3d = reshape(sig(:,sPos)',4,1,nmp);
else
    sig_3d = reshape(sig(:,sPos)',9,1,nmp);
end
if sum(max_4_YES)>=500
    cols_3d_4 = [s(max_4_YES),s(max_4_YES)+1,s(max_4_YES)+2,s(max_4_YES)+3,...
        s(max_4_YES)+4,s(max_4_YES)+5,s(max_4_YES)+6,s(max_4_YES)+7];
    G_3d_4 = reshape(G_all(:, cols_3d_4'), 4, 8, sum(max_4_YES));
    G_3d_4 = pagemtimes(dXdx(:,:,max_4_YES), G_3d_4);
    fp_3d_4 = squeeze(pagemtimes(G_3d_4,'ctranspose',sig_3d(:,:,max_4_YES),'none'));
    fp_3d_4 = reshape(vp_dets(max_4_YES),1,sum(max_4_YES)) .* fp_3d_4;
else
    max_4_YES = (nn_all == -1);
end

max_6_YES = (nn_all == 6);
if sum(max_6_YES)>=500
    cols_3d_6 = [s(max_6_YES),s(max_6_YES)+1,s(max_6_YES)+2,s(max_6_YES)+3,...
        s(max_6_YES)+4,s(max_6_YES)+5,s(max_6_YES)+6,s(max_6_YES)+7,...
        s(max_6_YES)+8,s(max_6_YES)+9,s(max_6_YES)+10,s(max_6_YES)+11];
    G_3d_6 = reshape(G_all(:, cols_3d_6'), 4, 12, sum(max_6_YES));
    G_3d_6 = pagemtimes(dXdx(:,:,max_6_YES), G_3d_6);
    fp_3d_6 = squeeze(pagemtimes(G_3d_6,'ctranspose',sig_3d(:,:,max_6_YES),'none'));
    fp_3d_6 = reshape(vp_dets(max_6_YES),1,sum(max_6_YES)) .* fp_3d_6;
else
    max_6_YES = (nn_all == -1);
end

max_9_YES = (nn_all == 9);
if sum(max_9_YES)>=500
    cols_3d_9 = [s(max_9_YES),s(max_9_YES)+1,s(max_9_YES)+2,s(max_9_YES)+3,...
        s(max_9_YES)+4,s(max_9_YES)+5,s(max_9_YES)+6,s(max_9_YES)+7,...
        s(max_9_YES)+8,s(max_9_YES)+9,s(max_9_YES)+10,s(max_9_YES)+11,...
        s(max_9_YES)+12,s(max_9_YES)+13,s(max_9_YES)+14,s(max_9_YES)+15,...
        s(max_9_YES)+16,s(max_9_YES)+17];
    G_3d_9 = reshape(G_all(:, cols_3d_9'), 4, 18, sum(max_9_YES));
    G_3d_9 = pagemtimes(dXdx(:,:,max_9_YES), G_3d_9);
    fp_3d_9 = squeeze(pagemtimes(G_3d_9,'ctranspose',sig_3d(:,:,max_9_YES),'none'));
    fp_3d_9 = reshape(vp_dets(max_9_YES),1,sum(max_9_YES)) .* fp_3d_9;
else
    max_9_YES = (nn_all == -1);
end

%% 循环赋值：仅组装内力、质量矩阵、质量比例阻尼矩阵
k_4_JS=1;  k_6_JS=1; k_9_JS=1;
for mp=1:nmp
    if max_4_YES(mp)
        fp = fp_3d_4(:,k_4_JS); k_4_JS = k_4_JS+1;
    elseif max_6_YES(mp)
        fp = fp_3d_6(:,k_6_JS); k_6_JS = k_6_JS+1;
    elseif max_9_YES(mp)
        fp = fp_3d_9(:,k_9_JS); k_9_JS = k_9_JS+1;
    else
        cols = s(mp):e(mp);
        G = dXdx(:,:,mp) * G_all(:, cols);
        fp = vp_dets(mp)*(G.'*sig(mp,sPos)');
    end

    % 保存物质点状态
    %mpData(mp).F    = F(:,:,mp);
    mpData(mp).sig  = Ksig(mp,:)';
    mpData(mp).epsE = epsE(mp,:)';
    mpData(mp).epsP = epsP(mp,:)';
    mpData(mp).eps  = eps_total(mp,:)';

    if ~isempty(Out)
        mpData(mp).modulus = modulus(mp);
        mpData(mp).STATEV = STATEV(mp,:)';
    end

    ed0=ed(s(mp):e(mp));
    npDoF=(nn_all(mp)*nD)^2;
    nnDoF=nn_all(mp)*nD;
    krow(npCnt+1:npCnt+npDoF)=repmat(ed0.',nnDoF,1);
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed0  ,nnDoF,1);

    % 显式动力：仅质量矩阵与纯质量比例阻尼
    m_Matrix = mpData(mp).mpM * mpData(mp).Svp / sum(mpData(mp).Svp) / nD;
    m_Matrix = diag(repelem(m_Matrix, nD));
    mval(npCnt+1:npCnt+npDoF)=m_Matrix;
    cval(npCnt+1:npCnt+npDoF)=c_par(2)*m_Matrix;
    fint(ed0)=fint(ed0)+fp;

    npCnt=npCnt+npDoF;
end

nDoF=length(V_old);
Stifiness.Kt=sparse(nDoF,nDoF);
Stifiness.Mt=sparse(krow,kcol,mval,nDoF,nDoF);
Stifiness.Ct=sparse(krow,kcol,cval,nDoF,nDoF);
end

function d = batch_det3x3(F)
N = size(F, 3);
f = reshape(F, 9, N);
d = (f(1,:).*f(5,:).*f(9,:) ...
   + f(4,:).*f(8,:).*f(3,:) ...
   + f(7,:).*f(2,:).*f(6,:) ...
   - f(7,:).*f(5,:).*f(3,:) ...
   - f(1,:).*f(8,:).*f(6,:) ...
   - f(4,:).*f(2,:).*f(9,:));
d = d(:);
end