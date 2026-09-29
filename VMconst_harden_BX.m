function [Dalg,sig,epsE,epsP] = VMconst_harden_BX(epsEtr,mCst,eps_total)

%von Mises linear elastic perfectly plastic constitutive model
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   16/05/2016
% Description:
% von Mises perfect plasticity constitutive model with an implicit backward
% Euler stress integration algorithm based on the following thesis:
%
% Coombs, W.M. (2011). Finite deformation of particulate geomaterials: 
% frictional and anisotropic Critical State elasto-plasticity. School of 
% Engineering and Computing Sciences. Durham University. PhD.
%
%--------------------------------------------------------------------------
% [Dalg,sigma,epsE] = VMCONST(epsEtr,mCst)
%--------------------------------------------------------------------------
% Input(s):
% epsEtr - trial elastic strain (6,1)
% mCst   - material constants 
%--------------------------------------------------------------------------
% Ouput(s);
% sig    - Cauchy stress (6,1)
% epsE   - elastic strain (6,1)
% Dalg   - algorithmic consistent tangent (6,6)
%--------------------------------------------------------------------------
% See also:
% YILEDFUNCDERIVATIVES - yield function 1st and 2nd derivatives
%--------------------------------------------------------------------------
E = mCst(:, 1);  
v = mCst(:, 2);  
rhoY = mCst(:, 3);
fc_k = mCst(:, 4);

rhoY = min(rhoY*1.20,max(rhoY*0.80,rhoY-sum(eps_total,2)*fc_k));
N = size(epsEtr, 1);
tol = 1e-9;
maxit = 5;

bm1 = [1; 1; 1; 0; 0; 0];

%% 1. 构造批量弹性刚度 De (6×6×N) 和柔度 Ce (6×6×N)
factor = E ./ ((1+v) .* (1-2*v));

% De = factor * (v*bm1*bm1' + (1-2*v)*diag([1,1,1,0.5,0.5,0.5]))
De = zeros(6, 6, N);
for k = 1:N
    De(:,:,k) = factor(k) * (v(k)*(bm1*bm1') + (1-2*v(k))*diag([1 1 1 0.5 0.5 0.5]));
end

% Ce = [-(v)*ones(3)+(1+v)*eye(3), 0; 0, 2*(1+v)*eye(3)] / E
Ce = zeros(6, 6, N);
for k = 1:N
    Ce_top = -v(k)*ones(3) + (1+v(k))*eye(3);
    Ce_bot = 2*(1+v(k))*eye(3);
    Ce(:,:,k) = [Ce_top, zeros(3); zeros(3), Ce_bot] / E(k);
end

%% 2. 弹性试算应力 sig = De * epsEtr
epsEtr_3 = reshape(epsEtr', 6, 1, N);  % 6×1×N
sig_3 = pagemtimes(De, epsEtr_3);       % 6×1×N
sig = reshape(sig_3, 6, N)';          % N×6

%% 3. 屈服判断
mean_sig = sum(sig(:, 1:3), 2) / 3;
s = sig - mean_sig * bm1';
J2 = (sum(s(:, 1:3).^2, 2) + sum(s(:, 4:6).^2, 2)) / 2;
f = sqrt(2*J2) ./ rhoY - 1;

%% 4. 初始化
epsE = epsEtr;          % N×6，默认保持试算值 (弹性情况)
Dalg = De;              % 6×6×N，默认弹性刚度
dgam = zeros(N, 1);     % N×1，塑性乘子初始化为0

% 逻辑掩码分离弹塑性点
plastic_mask = f > tol;     % N×1
elastic_mask = ~plastic_mask; % N×1

% 若全为弹性，直接计算 epsP 并返回
if ~any(plastic_mask)
    epsP_3 = pagemldivide(De, sig_3) - epsEtr_3;
    epsP = reshape(epsP_3, 6, N)';
    return;
end

%% 5. 塑性点 Newton-Raphson 同步迭代
% 提取塑性点子集 (仅对这些点进行迭代)
Np = sum(plastic_mask);
idx_p = find(plastic_mask);

% 初始化塑性点变量 (6×Np)
epsE_p = epsE(plastic_mask, :)';  % 初始为试算应变
dgam_p = zeros(1, Np);            % 1×Np，初始为0

% 初始残差 b (7×Np): 前6行为0，第7行为f
b = zeros(7, Np);
b(7, :) = f(plastic_mask)';

% 活动点标记 (未收敛的点)
active = true(1, Np);

for it = 1:maxit
    if ~any(active)
        break;
    end
    
    n_active = sum(active);
    idx_active_local = find(active);      % 在活动子集中的索引 (1~Np)
    idx_active_global = idx_p(active);    % 在全局 N 中的索引
    
    % 提取活动点数据
    epsE_active = epsE_p(:, active);     % 6×n_active
    dgam_active = dgam_p(active);         % 1×n_active
    
    % 计算当前应力 sig = De * epsE (仅活动塑性点)
    De_active = De(:,:,idx_active_global);
    sig_active = pagemtimes(De_active, reshape(epsE_active, 6, 1, n_active));
    sig_active = reshape(sig_active, 6, n_active);
    
    % 计算屈服函数导数 (仅活动点)
    [df, ddf] = yieldFuncDerivatives_batch(sig_active, rhoY(idx_active_global)');
    
    % 构造 Hessian A (7×7×n_active)
    A = zeros(7, 7, n_active);
    for k = 1:n_active
        i = idx_active_global(k);
        dg = dgam_active(k);
        % A = [I + dg*ddf*De, df; df'*De, 0]
        A(1:6, 1:6, k) = eye(6) + dg * ddf(:,:,k) * De(:,:,i);
        A(1:6, 7, k) = df(:, k);
        A(7, 1:6, k) = (df(:,k)' * De(:,:,i))';  % (df' * De)'
        A(7, 7, k) = 0;
    end
    
    % 求解线性系统 A * dx = -b
    b_active = b(:, active);
    dx = -pagemldivide(A, reshape(b_active, 7, 1, n_active));
    dx = reshape(dx, 7, n_active);
    
    % 更新未知数
    epsE_p(:, active) = epsE_p(:, active) + dx(1:6, :);
    dgam_p(active) = dgam_p(active) + dx(7, :);
    
    % 重新计算残差 (仅活动点)
    sig_active = pagemtimes(De(:,:,idx_active_global), reshape(epsE_p(:, active), 6, 1, n_active));
    sig_active = reshape(sig_active, 6, n_active);
    
    s_active = sig_active - sum(sig_active(1:3, :), 1)/3 .* bm1;
    J2_active = (sum(s_active(1:3,:).^2, 1) + sum(s_active(4:6,:).^2, 1))/2;
    f_new = sqrt(2*J2_active) ./ rhoY(idx_active_global)' - 1;
    
    [df_new, ~] = yieldFuncDerivatives_batch(sig_active, rhoY(idx_active_global)');
    
    % 更新残差 b
    b(:, active) = [epsE_p(:, active) - epsEtr(idx_active_global, :)' + ...
                    dgam_p(active) .* df_new; 
                    f_new];
    
    % 检查收敛
    res_norm = sqrt(sum(b(1:6, active).^2, 1));
    f_abs = abs(b(7, active));
    converged_now = (res_norm < tol) & (f_abs < tol);
    
    % 更新活动标记
    active(active) = ~converged_now;
end

% 回写塑性点结果到全局数组
epsE(plastic_mask, :) = epsE_p';
dgam(plastic_mask) = dgam_p';  % 保存最终 dgam 供后续使用

%% 6. 计算塑性点算法切线模量 Dalg (仅塑性点)
% Dalg = B(1:6,1:6), 其中 B = inv([Ce + dgam*ddf, df; df', 0])
for k = 1:Np
    i = idx_p(k);  % 全局索引
    % 重新计算最终状态的导数 (或从最后一次迭代保存)
    sig_final = epsE(i, :)';
    [df_final, ddf_final] = yieldFuncDerivatives_single(sig_final, rhoY(i));
    
    % 构造 KKT 矩阵并求逆
    KKT = [Ce(:,:,i) + dgam(i)*ddf_final, df_final; 
           df_final', 0];
    B = inv(KKT);
    Dalg(:,:,i) = B(1:6, 1:6);
end

%% 7. 计算 epsP = Dalg\sig - epsE (所有点)
epsP_3 = pagemldivide(Dalg, sig_3) - reshape(epsE', 6, 1, N);
epsP = reshape(epsP_3, 6, N)';

end

%% 辅助函数：批量 yieldFuncDerivatives
function [df, ddf] = yieldFuncDerivatives_batch(sig, rhoY)
% sig: 6×Np, rhoY: 1×Np
Np = size(sig, 2);
bm1 = [1; 1; 1; 0; 0; 0];

mean_sig = sum(sig(1:3, :), 1) / 3;
s = sig - bm1 * mean_sig;
J2 = (sum(s(1:3,:).^2, 1) + sum(s(4:6,:).^2, 1)) / 2;

dj2 = s;
dj2(4:6, :) = 2 * dj2(4:6, :);

df = dj2 ./ (rhoY .* sqrt(2*J2));

ddf = zeros(6, 6, Np);
for k = 1:Np
    term1 = ([eye(3)-ones(3)/3, zeros(3); zeros(3), 2*eye(3)]) / sqrt(2*J2(k));
    term2 = (dj2(:,k) * dj2(:,k)') / (2*J2(k))^(3/2);
    ddf(:,:,k) = (term1 - term2) / rhoY(k);
end
end

%% 辅助函数：单点 yieldFuncDerivatives (用于最后计算 Dalg)
function [df, ddf] = yieldFuncDerivatives_single(sig, rhoY)
bm1 = [1; 1; 1; 0; 0; 0];
s = sig - sum(sig(1:3))/3 * bm1;
J2 = (s'*s + s(4:6)'*s(4:6))/2;
dj2 = s; dj2(4:6) = 2*dj2(4:6);
ddj2 = [eye(3)-ones(3)/3, zeros(3); zeros(3), 2*eye(3)];
df = dj2 / (rhoY * sqrt(2*J2));
ddf = 1/rhoY * (ddj2/sqrt(2*J2) - (dj2*dj2')/(2*J2)^(3/2));
end