function [fint,mpData,mesh] = detMPs_Flow(H,mpData,mesh)

%Stiffness and internal force calculation for all material points
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to determine the stiffness contribution of a particle to the
% nodes that it influences based on a Updated Lagrangian finite deformation
% formulation.  The function also returns the stresses at the particles and
% the internal force contribution.  This function allows for elasto-
% plasticity at the material points.  The functionis applicable to 1, 2 and
% 3 dimensional problems without modification as well as different material
% point methods and background meshes.
%
%--------------------------------------------------------------------------
% [fint,Kt,mpData] = DETMPS(uvw,mpData)
%--------------------------------------------------------------------------
% Input(s):
% uvw    - nodal displacements that influence the MP (nn*nD,1)
% mpData - material point structured array. The following fields are
%          required by the function:
%           - dSvp  : basis function derivatives (nD,nn)
%           - nIN   : background mesh nodes associated with the MP (1,nn)
%           - Fn    : previous deformation gradient (3,3)
%           - epsEn : previous elastic logarithmic strain (6,1)
%           - mCst  : material constants
%           - vp    : material point volume (1)
%           - nSMe  : number stiffness matrix entries
% nD     - number of dimensions
%--------------------------------------------------------------------------
% Ouput(s);
% fint   - global internal force vector
% Kt     - global stiffness matrix
% mpData - material point structured array (see above).  The following
%          fields are updated by the function:
%           - F     : current deformation gradient (3,3)
%           - sig   : current Cauchy stress (6,1)
%           - epsE  : current elastic logarithmic strain (6,1)
%--------------------------------------------------------------------------
% See also:
% FORMULSTIFF      - updated Lagrangian material stiffness calculation
% HOOKE3D          - linear elastic constitutive model
% VMCONST          - von Mises elasto-plastic constitutive model
%--------------------------------------------------------------------------

nmp   = length(mpData);                                                     % number of material points
fint  = zeros(size(H));                                                     % zero internal force vector
nD = length(mpData(1).mpC);

npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]/nD^2);                                                 % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow; cval=krow;             % zero the stiffness information

for mp=1:nmp                                                                % material point loop
    %------------------- 通过H求解deltaH，类似于通过位移增量更新 ---------------------------------------------------
    ed = mpData(mp).nIN;        % 物质点关联的节点（类比位移求解的节点集）
    G = mpData(mp).dSvp;        % 形函数导数 ∇N (nD × nn)（类比位移求解的B矩阵）
    vp = mpData(mp).vp;         % 物质点体积
    
    % ------------------- 1. 水头增量 + 历史累加（核心：位移类比水头） -------------------
    dH_mp = G*H(ed);  % 本次分析步物质点水头增量（类比位移增量du，你之前用ddF，可替换）
    eps_FlownEtr = mpData(mp).epsn_Flow + dH_mp;  % 累加至历史值
    % 当前饱和度
    Flowpar_a = mpData(mp).Flow_CP(1)/1000;
    Flowpar_m = mpData(mp).Flow_CP(2);
    Flowpar_n = mpData(mp).Flow_CP(3);
    N   = mpData(mp).Svp;
    mpH = N*H(ed);
    Suction =    (mpData(mp).H0 + mpData(mp).dH_sum + mpData(mp).dH_sum_strain)*-9.81;
    H_suction0 = (mpData(mp).H0 + mpData(mp).dH_sum + mpH + mpData(mp).dH_sum_strain)*-9.81;
    H_suction = max(H_suction0,0.00);
    
    % % 历史水头梯度累加（epsn_Flow存储累积值，类比累积位移）
    % if mpData(mp).Flow_P > 0
    %     % 梯度修正
    %     L=mpData(mp).Flow_P/mpData(mp).Flow_SIZE;
    %     if strcmp(mpData(mp).Flow_Position,'Slope')
    %         eps_FlownEtr(1) = max(0,min(100,2*Suction/9.81/min(L,mpData(1).lp(2))));
    %         eps_FlownEtr(2) = max(0,min(100,2*Suction/9.81/min(L,mpData(1).lp(2))+1));
    %     else
    %         eps_FlownEtr(2) = max(0,min(100,2*Suction/9.81/min(L,mpData(1).lp(2))+1));
    %     end
    % end

    % 计算SWRC系数: m3 / e^{1/(m1*m2)}
    e = mpData(mp).e;
    coeff_swrc = Flowpar_a / (e^(1/(Flowpar_m*Flowpar_n)));

    % Sr = [1 + (suction / A)^{n}]^{-m}  ← 注意是除法
    Sr = (1.0 + (H_suction ./ coeff_swrc).^Flowpar_n).^(-Flowpar_m);
    Sr =  min(1.0,Sr);

    if H_suction0 <= 0.000
        mpData(mp).Saturated_YES = 1.0;
    else
        mpData(mp).Saturated_YES = 0.0;
    end

    % ------------------- 2. 非饱和渗流本构（保留你的逻辑） -------------------
    % 相对渗透系数（VG模型，保留）
    kr = max(0.01,Sr^0.5 *(1-(1-Sr^(1/Flowpar_m))^Flowpar_m)^2);
    
    % if Sr==1.0 && mpData(mp).mpC(2)<0.025
    %     % 坐标
    %     mpC = mpData(mp).mpC;
    %     LS = (mpC-mesh.coord);
    %     LS = LS(:,1).^2+LS(:,2).^2;
    %     [~,N] = min(LS);
    %     bc_Flow = unique([mesh.bc_Flow(:,1);N]);
    %     mesh.bc_Flow = [bc_Flow,zeros(size(bc_Flow))];
    % end
    
    % 渗透系数张量（类比弹性模量，负号适配达西定律，保留你的对角矩阵）
    K_mp = diag(mpData(mp).Ksat*kr);

    if mpData(mp).Flow_P~=0
        % 梯度修正
        L=mpData(mp).Flow_P/mpData(mp).Flow_SIZE;
        k_SLOW = max(0.01,mpData(mp).Sr^0.5 *(1-(1-mpData(mp).Sr^(1/Flowpar_m))^Flowpar_m)^2)*mpData(mp).Ksat(2);
        eps_FlownEtr0 = max(-1,min(100,2*Suction/9.81/min(L,mpData(1).lp(2))+1));
        % 当前水头
        mpData(mp).Flow_P_current=min(mpData(mp).Flow_P,k_SLOW*L*eps_FlownEtr0);
    else
        mpData(mp).Flow_P_current=0.0;
    end

    if mpData(mp).Saturated_YES == 1.0
        
        c_scalar = 0.0;
    else
        % 持水度 = e/(1+e) * (m*n/A) * Sr^{...} * (...)  ← 注意是除以A
        c_scalar = (e/(e+1)) * (Flowpar_m * Flowpar_n / coeff_swrc) * ...
            Sr.^((Flowpar_m+1)/Flowpar_m) .* (Sr.^(-1/Flowpar_m) - 1).^((Flowpar_n-1)/Flowpar_n);
    end

    % ------------------- 3. 渗流"节点力"（流量）计算（类比固体力学节点力） -------------------
    % 类比：固体力学内力=B^TσV → 渗流流量=G^T*(K∇H)*V
    mpV = K_mp * eps_FlownEtr;  % 类比应力σ（渗流中为K∇H，即"渗流应力"）
    fp = vp * (G' * mpV);  % 节点流量（类比节点内力，保留你的detF）

    % ------------------- 4. 容水度矩阵（类比阻尼矩阵，保留你的公式） -------------------
    % 第二步：行求和（集中质量法核心，只加这1行！）
    % 或者一致形式（方式2，但通常需要 lumping 简化）
    cp_consistent = vp * c_scalar * (mpData(mp).Svp' * mpData(mp).Svp);   % nn×nn矩阵
    cp = diag(sum(cp_consistent, 2)); % 行求和 lumping

    % ------------------- 5. 刚度矩阵（类比切线刚度，保留你的组装） -------------------
    kp = vp * (G' * K_mp * G);  % 渗透刚度（类比固体力学切线刚度）

    % ------------------- 6. 矩阵组装（完全保留你的逻辑） -------------------
    npDoF=(size(ed,1)*size(ed,2))^2;        
    nnDoF=size(ed,1)*size(ed,2);            
    krow(npCnt+1:npCnt+npDoF)=repmat(ed.',nnDoF,1);         
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed  ,nnDoF,1);         
    kval(npCnt+1:npCnt+npDoF)=kp;           % 渗透刚度（类比切线刚度）
    cval(npCnt+1:npCnt+npDoF)=cp;     % 容水度（类比阻尼矩阵）
    npCnt=npCnt+npDoF;                      
    fint(ed)=fint(ed)+fp;                   % 内流量（类比节点内力）

    % ------------------- 7. 保存当前状态（类比位移求解的应力/应变保存） -------------------
    mpData(mp).eps_Flow  = eps_FlownEtr;    % 本次水头增量（类比应变增量）
    mpData(mp).sig_Flow  = mpV;             % 本次"渗流应力"（K∇H，类比Cauchy应力）
end

% 构建全局矩阵（类比固体力学的整体刚度/阻尼矩阵）
nDoF=length(H);                                                             
mesh.Flow.Kt=sparse(krow,kcol,kval,nDoF,nDoF);  % 全局渗透刚度（类比整体刚度）
mesh.Flow.Ct=sparse(krow,kcol,cval,nDoF,nDoF);  % 全局容水度（类比整体阻尼）
end