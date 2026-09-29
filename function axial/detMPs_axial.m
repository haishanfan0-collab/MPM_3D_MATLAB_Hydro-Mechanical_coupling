function [fint,mpData,Stifiness] = detMPs_axial(uvw,mpData,c_par,Cal_par)

%Stiffness and internal force calculation for axisymmetric material points
%--------------------------------------------------------------------------
% FINAL VERSION: 
% 1. Reconstructs A matrix from 3D to axisymmetric format [1 5 9 4]
% 2. Merged static/dynamic with only final assembly differing
%--------------------------------------------------------------------------

nmp   = length(mpData);                                                     
fint  = zeros(size(uvw));                                                   

npCnt = 0;                                                                  
tnSMe = sum([mpData.nSMe]);                                                 
krow  = zeros(tnSMe,1); kcol=krow; kval=krow; mval=krow; cval=krow;         

nD = length(mpData(1).mpC);

% 轴对称索引（F11=rr, F22=zz, F33=theta, F12=rz）
aPos=[1 2 3 4];        
sPos=[1 2 3 4];        

%% 
for mp=1:nmp                                                                

    nIN = mpData(mp).nIN;                                                   
    dNx = mpData(mp).Gvp;       % [dN/dr; dN/dz] (2×nn)
    Tvp = mpData(mp).Tvp;        % N/r (1×nn)
    nn  = size(dNx,2);                                                      
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                    
    ed  = reshape(ed,1,nn*nD);                                              

    % 当前粒子半径
    rp = mpData(mp).mpC(1);
    lp = mpData(mp).lp(1);
    
    %==========================================================================
    % G矩阵组装（4行：rr, zz, theta, rz）
    %==========================================================================
    G=zeros(4,nD*nn);
    G(1,1:nD:end)=dNx(1,:);        % rr
    G(2,2:nD:end)=dNx(2,:);        % zz  
    G(3,1:nD:end)=Tvp;             % theta = u/r
    G(4,1:nD:end)=dNx(2,:);        % rz
    G(4,2:nD:end)=dNx(1,:);

    %==========================================================================
    % 变形梯度计算
    %==========================================================================
    u_vec = uvw(ed);
    u_radial = u_vec(1:2:end);
    u_axial  = u_vec(2:2:end);
    
    eps_rr = dNx(1,:) * u_radial;
    eps_zz = dNx(2,:) * u_axial;
    eps_theta = Tvp * u_radial;
    
    ddF = zeros(3);
    ddF(1,1) = eps_rr;         % F11
    ddF(2,2) = eps_zz;         % F22
    ddF(1,2) = dNx(2,:) * u_radial;  % F12
    ddF(2,1) = dNx(1,:) * u_axial;   % F21
    ddF(3,3) = eps_theta;      % F33

    dF = eye(3) + ddF;
    F  = dF*mpData(mp).Fn;
    
    %% 应变计算（UL格式） 仅弹性应变
    epsEn  = mpData(mp).epsEn; 
    epsEn(4:6) = 0.5*epsEn(4:6);
    epsEn  = epsEn([1 4 6; 4 2 5; 6 5 3]);
    [V,D]  = eig(epsEn);
    BeT    = dF*(V*diag(exp(2*diag(D)))*V.')*dF.';
    [V,D]  = eig(BeT);
    epsEtr = 0.5*V*diag(log(diag(D)))*V.';
    epsEtr = diag([1 1 1 2 2 2])*epsEtr([1 5 9 2 6 3]).';

    %==========================================================================
    %% 总应变计算 (Hencky 对数应变，直接由 F 计算)
    %==========================================================================
    % 右 Cauchy-Green 张量
    C = F'*F;
    
    % 谱分解
    [V, D] = eig(C);
    
    % 特征值处理 (避免数值误差)
    lambda_sq = max(diag(D), 1e-20);
    
    % 主对数应变 (Hencky 应变主值)
    eps_principal = 0.5*log(lambda_sq);
    
    % 组装完整对数应变张量
    eps_tensor = V*diag(eps_principal)*V';
    
    % 转换为 Voigt 格式 (ABAQUS顺序: 11,22,33,12,13,23)
    % 先提取张量分量，后2倍得到工程剪应变
    eps_total = [eps_tensor(1,1);      % rr
                 eps_tensor(2,2);      % zz  
                 eps_tensor(3,3);      % theta-theta
                 eps_tensor(1,2);      % rz (张量分量)
                 eps_tensor(1,3);      % r-theta (通常为0，轴对称)
                 eps_tensor(2,3)];     % z-theta (通常为0，轴对称)
             
    % 转换为工程剪应变 (2*张量剪应变)
    eps_total(4:6) = 2*eps_total(4:6);

    %==========================================================================
    %% 本构计算（静力/动力共用）
    %==========================================================================
    if mpData(mp).cmType == 1
        [D,Ksig,epsE,epsP]=Hooke3d(epsEtr,mpData(mp));
    elseif mpData(mp).cmType == 2
        [D,Ksig,epsE,epsP]=VMconst(epsEtr,mpData(mp));
    elseif mpData(mp).cmType == 3
        [D,Ksig,epsE,epsP]=VMconst_harden(epsEtr,mpData(mp));
    elseif mpData(mp).cmType == 4
        [D,Ksig,epsE,epsP]=VMconst_FLOW(epsEtr,mpData(mp));
    elseif mpData(mp).cmType == 5
        [D,Ksig,epsE,epsP,mpData(mp)] = Boundary_module(eps_total,mpData(mp),Cal_par);
    end

    sig = Ksig/det(F);
    A = formULstiff(F,D,sig,BeT);  % 这是为3D设计的A矩阵

    %==========================================================================
    % Updated Lagrangian 映射
    %==========================================================================
    iF   = dF\eye(3);
    dXdx = [iF(1) 0     0     iF(2) 0     0     0     0     iF(3) ;
            0     iF(5) 0     0     iF(4) iF(6) 0     0     0     ;
            0     0     iF(9) 0     0     0     iF(8) iF(7) 0     ;
            iF(4) 0     0     iF(5) 0     0     0     0     iF(6) ;
            0     iF(2) 0     0     iF(1) iF(3) 0     0     0     ;
            0     iF(8) 0     0     iF(7) iF(9) 0     0     0     ;
            0     0     iF(6) 0     0     0     iF(5) iF(4) 0     ;
            0     0     iF(3) 0     0     0     iF(2) iF(1) 0     ;
            iF(7) 0     0     iF(8) 0     0     0     0     iF(9)];
    
    G_transformed = dXdx(aPos,aPos)*G;
    
    % 轴对称体积（使用当前构型半径）
    axisym_factor = 2*pi*max(rp, lp*0.01);  
    current_volume = axisym_factor * mpData(mp).vp * det(dF);
    
    % 刚度和内力（静力/动力共用）
    kp = current_volume * (G_transformed.'*A(aPos,aPos)*G_transformed);
    fp = current_volume * (G_transformed.'*sig(sPos));

    % 存储更新
    mpData(mp).F    = F;
    mpData(mp).sig  = Ksig;
    mpData(mp).epsE = epsE;
    mpData(mp).epsP = epsP;

    %==========================================================================
    % 矩阵组装（静力/动力分支）
    %==========================================================================
    npDoF=(size(ed,1)*size(ed,2))^2;
    nnDoF=size(ed,1)*size(ed,2);
    
    krow(npCnt+1:npCnt+npDoF)=repmat(ed.',nnDoF,1);
    kcol(npCnt+1:npCnt+npDoF)=repmat(ed,nnDoF,1);
    kval(npCnt+1:npCnt+npDoF)=kp;
    
    % 动力模式下额外计算质量矩阵和阻尼矩阵
    if strcmp(Cal_par.method, 'dynamic')
        mpMass = mpData(mp).mpM * axisym_factor;
        m_vec = mpMass * mpData(mp).Svp / sum(mpData(mp).Svp);
        m_Matrix = diag(repelem(m_vec, nD));
        c_Matrix = c_par(1)*kp + c_par(2)*m_Matrix;
        
        mval(npCnt+1:npCnt+npDoF)=m_Matrix;
        cval(npCnt+1:npCnt+npDoF)=c_Matrix;
    end
    
    npCnt=npCnt+npDoF;
    fint(ed)=fint(ed)+fp;
end

nDoF=length(uvw);
Stifiness.Kt=sparse(krow,kcol,kval,nDoF,nDoF);
Stifiness.Mt=sparse(krow,kcol,mval,nDoF,nDoF);
Stifiness.Ct=sparse(krow,kcol,cval,nDoF,nDoF);
end