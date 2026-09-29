function [fint,mpData,Stifiness] = detMPs(uvw,mpData,c_par,Cal_par)

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
fint  = zeros(size(uvw));                                                   % zero internal force vector

npCnt = 0;                                                                  % counter for the number of entries in Kt
tnSMe = sum([mpData.nSMe]);                                                 % total number of stiffness matrix entries
krow  = zeros(tnSMe,1); kcol=krow; kval=krow; mval=krow; cval=krow;         % zero the stiffness information
ddF   = zeros(3,3);                                                           % derivative of duvw wrt. spatial position

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

for mp=1:nmp                                                                % material point loop
    nIN = mpData(mp).nIN;                                                   % nodes associated with the material point
    dNx = mpData(mp).dSvp;                                                  % basis function derivatives (start of lstp)
    nn  = size(dNx,2);                                                      % no. dimensions & no. nodes
    ed  = repmat((nIN-1)*nD,nD,1)+repmat((1:nD).',1,nn);                    % degrees of freedom of nodes (matrix form)
    ed  = reshape(ed,1,nn*nD);                                              % degrees of freedom of nodes (vector form)

    if nD==1                                                                % 1D case
        G=dNx;                                                              % strain-displacement matrix
    elseif nD==2                                                            % 2D case (plane strain & stress)
        G=zeros(4,nD*nn);                                                   % zero the strain-disp matrix (2D)
        G([1 3],1:nD:end)=dNx;                                              % strain-displacement matrix
        G([4 2],2:nD:end)=dNx;
    else                                                                    % 3D case
        G=zeros(9,nD*nn);                                                   % zero the strain-disp matrix (3D)
        G([1 4 9],1:nD:end)=dNx;                                            % strain-displacement matrix
        G([5 2 6],2:nD:end)=dNx;
        G([8 7 3],3:nD:end)=dNx;
    end

    ddF(fPos) = G*uvw(ed);                                                  % spatial gradient (start of lstp) of displacements
    dF     = (eye(3)+ddF);                                                  % deformation gradient increment
    F      = dF*mpData(mp).Fn;                                              % deformation gradient
    epsEn  = mpData(mp).epsEn; epsEn(4:6) = 0.5*epsEn(4:6);                 % previous elastic strain
    epsEn  = epsEn([1 4 6; 4 2 5; 6 5 3]);                                  % matrix form of previous elastic strain
    [V,D]  = eig(epsEn);                                                    % eigen values and vectors of the elastic strain
    BeT    = dF*(V*diag(exp(2*diag(D)))*V.')*dF.';                          % trial left Cauchy-Green strain
    [V,D]  = eig(BeT);                                                      % eigen values and vectors of the trial left Cauchy-Green strain
    epsEtr = 0.5*V*diag(log(diag(D)))*V.';                                  % trial elastic strain (tensor form)
    epsEtr = diag([1 1 1 2 2 2])*epsEtr([1 5 9 2 6 3]).';                   % trial elastic strain (vector form)

    %==========================================================================
    %% 1. 计算当前总应变 (Hencky 对数应变)
    %==========================================================================
    C = F' * F;                                 % 右 Cauchy-Green
    [V, D] = eig(C);
    lambda_sq = max(real(diag(D)), 1e-20);

    % 主对数应变
    eps_principal = 0.5*log(lambda_sq);         % ln(λ)

    % 组装张量 → Voigt 格式 [11,22,33,12,13,23]
    eps_tensor = V*diag(eps_principal)*V';
    eps_total = [eps_tensor(1,1); eps_tensor(2,2); eps_tensor(3,3);
        eps_tensor(1,2); eps_tensor(1,3); eps_tensor(2,3)];
    eps_total(4:6) = 2*eps_total(4:6);          % 工程剪应变

    %---------------------------------------------------------------------- % Constitutive model
    if mpData(mp).cmType == 1
        [D,Ksig,epsE,epsP]=Hooke3d(epsEtr,mpData(mp));                      % elastic behaviour
    elseif mpData(mp).cmType == 2
        [D,Ksig,epsE,epsP]=VMconst(epsEtr,mpData(mp));                      % elasto-plastic behaviour (von Mises)
    elseif mpData(mp).cmType == 3
        [D,Ksig,epsE,epsP]=VMconst_harden(epsEtr,mpData(mp));               % elasto-plastic behaviour (von Mises with harden)
    elseif mpData(mp).cmType == 4
        [D,Ksig,epsE,epsP]=VMconst_FLOW(epsEtr,mpData(mp));               % elasto-plastic behaviour (von Mises with harden)
    elseif mpData(mp).cmType == 5
        [D,Ksig,epsE,epsP,mpData(mp)] = Boundary_module(eps_total,mpData(mp),Cal_par);
    end
    %----------------------------------------------------------------------

    sig = Ksig/det(F);                                                      % Cauchy stress
    A   = formULstiff(F,D,sig,BeT);                                         % spatial tangent stiffness matrix

    iF   = dF\eye(3);                                                       % inverse deformation gradient increment
    dXdx = [iF(1) 0     0     iF(2) 0     0     0     0     iF(3) ;         % start of loadstep to current configuration
        0     iF(5) 0     0     iF(4) iF(6) 0     0     0     ;         % derivative mapping matrix
        0     0     iF(9) 0     0     0     iF(8) iF(7) 0     ;
        iF(4) 0     0     iF(5) 0     0     0     0     iF(6) ;
        0     iF(2) 0     0     iF(1) iF(3) 0     0     0     ;
        0     iF(8) 0     0     iF(7) iF(9) 0     0     0     ;
        0     0     iF(6) 0     0     0     iF(5) iF(4) 0     ;
        0     0     iF(3) 0     0     0     iF(2) iF(1) 0     ;
        iF(7) 0     0     iF(8) 0     0     0     0     iF(9)];
    G  = dXdx(aPos,aPos)*G;                                                 % derivatives of basis functions (current)

    % 节点刚度矩阵
    kp = mpData(mp).vp*det(dF)*(G.'*A(aPos,aPos)*G);                        % material point stiffness contribution
    % 节点应力矩阵
    fp = mpData(mp).vp*det(dF)*(G.'*sig(sPos));                             % internal force contribution

    mpData(mp).F    = F;                                                    % store deformation gradient
    mpData(mp).sig  = Ksig;                                                  % store Cauchy stress
    mpData(mp).epsE = epsE;                                                 % store elastic logarithmic strain
    mpData(mp).epsP = epsP;

    if strcmp(Cal_par.method, 'static')
        % 矩阵组装
        npDoF=(size(ed,1)*size(ed,2))^2;                                        % no. entries in kp
        nnDoF=size(ed,1)*size(ed,2);                                            % no. DoF in kp
        krow(npCnt+1:npCnt+npDoF)=repmat(ed.',nnDoF,1);                         % row position storage
        kcol(npCnt+1:npCnt+npDoF)=repmat(ed  ,nnDoF,1);                         % column position storage
        kval(npCnt+1:npCnt+npDoF)=kp;                                           % stiffness storage
        npCnt=npCnt+npDoF;                                                      % number of entries in Kt
        fint(ed)=fint(ed)+fp;                                                   % internal force contribution

    elseif strcmp(Cal_par.method, 'dynamic')
        % 节点质量矩阵
        m_Matrix = mpData(mp).mpM * mpData(mp).Svp / sum(mpData(mp).Svp) / nD;
        m_Matrix = diag(repelem(m_Matrix, nD));
        % 阻尼矩阵
        c_Matrix = c_par(1)*kp + c_par(2)*m_Matrix;

        % 矩阵组装
        npDoF=(size(ed,1)*size(ed,2))^2;                                        % no. entries in kp
        nnDoF=size(ed,1)*size(ed,2);                                            % no. DoF in kp
        krow(npCnt+1:npCnt+npDoF)=repmat(ed.',nnDoF,1);                         % row position storage
        kcol(npCnt+1:npCnt+npDoF)=repmat(ed  ,nnDoF,1);                         % column position storage
        kval(npCnt+1:npCnt+npDoF)=kp;                                           % stiffness storage
        mval(npCnt+1:npCnt+npDoF)=m_Matrix;                                     % mass storage
        cval(npCnt+1:npCnt+npDoF)=c_Matrix;                                     % damp storage
        npCnt=npCnt+npDoF;                                                      % number of entries in Kt
        fint(ed)=fint(ed)+fp;                                                   % internal force contribution
    else
        error('method only can be ''dynamic'' and ''static'' !!!')
    end
end

nDoF=length(uvw);                                                           % number of degrees of freedom
Stifiness.Kt=sparse(krow,kcol,kval,nDoF,nDoF);                                        % form the global stiffness matrix
Stifiness.Mt=sparse(krow,kcol,mval,nDoF,nDoF);                                        % form the global maxx matrix
Stifiness.Ct=sparse(krow,kcol,cval,nDoF,nDoF);                                        % form the global damp matrix
end