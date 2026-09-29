function [mpData] = updateMPs(uvw,duvw,dduvw,mpData,mesh)

%Material point update: stress, position and volume
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   29/01/2019
% Description:
% Function to update the material point positions and volumes (and domain
% lengths for GIMPM).  The function also updates the previously converged
% value of the deformation gradient and the logarithmic elastic strain at
% each material point based on the converged value and calculates the total
% displacement of each material point.
%
% For the generalised interpolation material point method the domain
% lengths are updated according to the stretch tensor following the
% approach of:
% Charlton, T.J., Coombs, W.M. & Augarde, C.E. (2017). iGIMP: An implicit
% generalised interpolation material point method for large deformations.
% Computers and Structures 190: 108-125.
%
%--------------------------------------------------------------------------
% [mpData] = UPDATEMPS(uvw,mpData)
%--------------------------------------------------------------------------
% Input(s):
% uvw    - nodal displacements (nodes*nD,1)
% mpData - material point structured array.  The function requires:
%           - mpC : material point coordinates
%           - Svp : basis functions
%           - F   : deformation gradient
%           - lp0 : initial domain lenghts (GIMPM only)
%--------------------------------------------------------------------------
% Ouput(s);
% mpData - material point structured array.  The function modifies:
%           - mpC   : material point coordinates
%           - vp    : material point volume
%           - epsEn : converged elastic strain
%           - Fn    : converged deformation gradient
%           - u     : material point total displacement
%           - lp    : domain lengths (GIMPM only)
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------
t = [1 5 9];                                                                % stretch components for domain updating
nmp = length(mpData);                                                       % number of material points
nD  = length(mpData(1).mpC);                                                % number of dimensions
Name_mp = fieldnames(mpData);

for mp=1:nmp
    nIN = mpData(mp).nIN;                                                   % nodes associated with material point
    nn  = length(nIN);                                                      % number nodes
    N   = mpData(mp).Svp;                                                   % basis functions
    F   = mpData(mp).F;                                                     % deformation gradient
    ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom

    mpU = N*uvw(ed);                                                        % material point displacement
    dmpU = N*duvw(ed);                                                      % material point displacement
    ddmpU = N*dduvw(ed);                                                    % material point displacement

    mpData(mp).mpC   = mpData(mp).mpC + mpU;                                % update material point coordinates
    mpData(mp).vp    = det(F)*mpData(mp).vp0;                               % update material point volumes
    mpData(mp).epsEn = mpData(mp).epsE;                                     % update material point elastic strains
    mpData(mp).epsPn = mpData(mp).epsP;                                     % update material point plastic strains
    mpData(mp).Fn    = mpData(mp).F;                                        % update material point deformation gradients

    mpData(mp).u     = mpData(mp).u + mpU.';                                % update material point displacements
    mpData(mp).du     = dmpU';
    mpData(mp).ddu     = ddmpU';

    mpData(mp).sign = mpData(mp).sig;
    if mpData(mp).mpType == 2                                               % GIMPM only (update domain lengths)
        [V,D] = eig(F.'*F);                                                 % eigen values and vectors F'F
        U     = V*sqrt(D)*V.';                                              % material stretch matrix
        mpData(mp).lp = (mpData(mp).lp0).*U(t(1:nD));                       % update domain lengths
        mpData(mp).lp = min(mpData(mp).lp,mpData(mp).lp0*5.00);
    end

    try
        mpData(mp).STATEVn = mpData(mp).STATEV;
    catch
    end

    %% 渗流模块更新
    if sum(strcmp(Name_mp,'Flow_Type'))>0 && mpData(mp).Gravity_Water_Content>0.0
        mpH = N*mesh.H(nIN);                                                % material point displacement
        dmpH = N*mesh.HVA.dHdt(nIN);                                        % material point displacement
        mpData(mp).dH   = dmpH;
        mpData(mp).epsn_Flow = mpData(mp).eps_Flow;                         % update material point elastic strains
        mpData(mp).sign_Flow = mpData(mp).sig_Flow;                         % update material point plastic strains=

        % 更新水头
        mpData(mp).dH_sum    = mpData(mp).dH_sum + mpH;                     % update material point displacements
        % mpData(mp).dH_sum    = min(-mpData(mp).H0-mpData(mp).dH_sum_strain,mpData(mp).dH_sum);

        mpData(mp).Scution_Flow = (mpData(mp).H0 + mpData(mp).dH_sum + mpData(mp).dH_sum_strain)*-9.81;
        mpData(mp).Scution_Flow = max(mpData(mp).Scution_Flow,0.00);
        
        if mpData(mp).Flow_Type == 99
            mpData(mp).Scution_Flow = 0.0;
            mpData(mp).dH_sum    = 0.0;
            mpData(mp).Sr = 1.00;
            mpData(mp).dH = 0;
        end

        Fa = mpData(mp).Flow_CP(1)/1000;
        Fm = mpData(mp).Flow_CP(2);
        Fn = mpData(mp).Flow_CP(3);
        e = mpData(mp).e;  % 提取孔隙比

        % 计算SWRC系数: m3 / e^{1/(m1*m2)}
        coeff_swrc = Fa / (e^(1/(Fm*Fn)));

        % 正确的饱和度计算（除法，不是乘法！）
        suction = mpData(mp).Scution_Flow;  % 注意拼写错误
        mpData(mp).Sr = (1.0 + (suction ./ coeff_swrc).^Fn).^(-Fm);
        mpData(mp).Sr = min(1.0, max(0.01, mpData(mp).Sr));  % 添加下限保护

        if mpData(mp).Sr >= 0.999
            mpData(mp).Saturated_YES = 1.0;
        else
            mpData(mp).Saturated_YES = 0.0;
        end

        % 更新干密度
        mpData(mp).Roud = mpData(mp).Gs/(1+mpData(mp).e);

        % 计算质量含水率
        mpData(mp).Gravity_Water_Content = (mpData(mp).Sr .* mpData(mp).e) ./ mpData(mp).Gs*1000;

        % 更新物质点质量
        mpData(mp).mpM = mpData(mp).vp*(1+mpData(mp).Gravity_Water_Content)*mpData(mp).Roud;

        try
            mpData(mp).STATEVn(7) = mpData(mp).Scution_Flow;
        catch
        end

    end

end
end