function [mesh,mpData,HVA,HVA_old,oobf_Flow,fd_Flow,H,frct_Flow,Cal_par] = ...
          AMPLE_Solve_Flow_Main(mesh,mpData,HVA_old,...  % 物质点/网格参数
          oobf_Flow,fd_Flow,...                          % 外力数据
          H,frct_Flow,...                                % 核心出装
          Cal_par)                                       % 计算参数
    [nodes,~] = size(mesh.coord);
    
    %% 关键：获取节点deltaH / 考虑由于变形引起的水头变化，每一次的 dH/dx 和 dH/dy 都是现算的
    if Cal_par.NRit_Flow==0
        nmp   = length(mpData);                                                     % number of material points
        for mp=1:nmp                                                                % material point loop
            %------------------- 通过H求解deltaH，类似于通过位移增量更新 ---------------------------------------------------
            ed = mpData(mp).nIN;        % 物质点关联的节点（类比位移求解的节点集）
            G = mpData(mp).dSvp;        % 形函数导数 ∇N (nD × nn)（类比位移求解的B矩阵）
            mpData(mp).epsn_Flow = G*HVA_old.H(ed);  % 根据总水头获取水头变化率
            % mpData(mp).epsn_Flow(2) = 1 + mpData(mp).epsn_Flow(2);
        end
    end

    %% 求解渗流
    [dH,dQ,HVA.dHdt] = model_Solve_Flow(mesh.bc_Flow, mesh.Flow.Kt, mesh.Flow.Ct, oobf_Flow,...
                                          Cal_par.NRit_Flow, fd_Flow, Cal_par.dt,...
                                          H, HVA_old.dHdt);            % linear solver

    Cal_par.NRit_Flow = Cal_par.NRit_Flow+1;                                           % increment the NR counter
    H  = H+dH;                                                                         % update displacements
    frct_Flow = frct_Flow+dQ;                                                          % update reaction forces
    [Qint,mpData,mesh] = detMPs_Flow(H,mpData,mesh);                                   % global stiffness & internal force
    
    Q_fext = detExtFlaw(nodes,mpData);
    Q_FORCE = Q_fext+frct_Flow;

    fint_all = (Qint+mesh.Flow.Ct*HVA.dHdt);
    % 
    % check = [mesh.coord,Q_FORCE,Qint,fint_all,HVA.dHdt,H];
    % check_Down=check(check(:,2)~=0,:);
    % check_Left=check(check(:,1)~=0,:);
    Cal_par.fErr_Flow = norm(Q_FORCE-fint_all)/norm(Q_FORCE+eps);                      % normalised oobf error
    % disp(num2str( Cal_par.fErr_Flow ))
    if Cal_par.NRit_Flow<5
        oobf_Flow = Q_FORCE-fint_all;
    else
        oobf_Flow = oobf_Flow*0.1+(Q_FORCE-fint_all)*0.9;
    end
fprintf(1,'%s %2i %s: %8.3e, %s: %8.3e, %s: %8.3e \n','  iteration',Cal_par.NRit_Flow,' NR error',Cal_par.fErr_Flow,'max(oobf)',norm(abs(Q_FORCE-fint_all)), 'max(H)', norm(abs(H)));   % text output to screen (NR error)
end


%% 求解主程序
function [dH,dQ,dHdt] = model_Solve_Flow(bc, Kt, Ct, oobf, NRit, fd, dt, H, dHdt_old)
%Linear solver
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   23/01/2019
% Description:
% Function to solve the linear system of equations for the increment in
% displacements and reaction forces.  The linear system is only solved for
% the first Newton-Raphson iteration (NRit>0) onwards as the zeroth
% iteration is required to construct the stiffness matrix based on the
% positions of the material points at the start of the loadstep.  This is
% different from the finite element method where the stiffness matrix from
% the last iteration from the previous loadstep can be used for the zeroth
% iteration. 
%
% In the case of non-zero displacement boundary conditions, these are
% applied when NRit = 1 and then the displacements for these degrees of
% freedom are fixed for the remaining iterations.
%
%--------------------------------------------------------------------------
% [duvw,drct] = LINSOLVE(bc,Kt,oobf,NRit,fd)
%--------------------------------------------------------------------------
% Input(s):
% bc    - boundary conditions (*,2)
% Kt    - global stiffness matrix (nDoF,nDoF)
% oobf  - out of balance force vector (nDoF,1)
% NRit  - Newton-Raphson iteration counter (1)
% fd    - free degrees of freedom (*,1)
%--------------------------------------------------------------------------
% Ouput(s);
% duvw  - displacement increment (nDoF,1)
% drct  - reaction force increment (nDoF,1)
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

nDoF = length(oobf);                                                        % number of degrees of freedom 
dH = zeros(nDoF,1);                                                         % zero reaction increment
dQ = zeros(nDoF,1);

if (NRit)>0
    ga = 1.00;
    dH(bc(:,1))=(1+sign(1-NRit))*bc(:,2);                                   % apply non-zero boundary conditions

    % --- 1. 等效刚度矩阵 ---
    K_hat = Ct/(ga*dt) + Kt;

    % --- 3. 边界条件处理 ---
    F_hat_mod = oobf - K_hat(:, bc(:,1)) * dH(bc(:,1));

    % 求解自由自由度上的位移
    dH(fd) = K_hat(fd, fd) \ F_hat_mod(fd);

    % --- 5. 更新加速度和速度 (必须使用求得的 u_new) ---
    dHdt = (H+dH)/(ga*dt);
    % --- 6. 计算约束自由度上的总反力 (可选) ---
    dQ(bc(:,1))=K_hat(bc(:,1),:)*dH-oobf(bc(:,1));                          % determine reaction forces 
else
    dHdt = dHdt_old;                                                        % v_{n+1}
end

end