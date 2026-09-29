function [mesh,mpData,VA_old,NRit,uvw,frct,fd,bounday_imfo] = updataLS_mp(duvw,dduvw,mpData,mesh,Cal_par,uvw)

FHS_LS = [];
method = Cal_par.method;
nmp = length(mpData);                                                       % number of material points
nD  = length(mpData(1).mpC);                                                % number of dimensions

if strcmp(method, 'dynamic')
    for mp=1:nmp
        nIN = mpData(mp).nIN;                                                   % nodes associated with material point
        nn  = length(nIN);                                                      % number nodes
        N   = mpData(mp).Svp;                                                   % basis functions
        ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom

        dmpU = N*duvw(ed);                                                      % material point displacement
        ddmpU = N*dduvw(ed);                                                    % material point displacement

        if max(abs(dmpU))>1e2 && max(abs(ddmpU))>1e3
            FHS_LS(end+1)=mp;
        end
    end
end

if strcmp(method, 'static')
    for mp=1:nmp
        F   = mpData(mp).F;
        if mpData(mp).mpType == 2                                               % GIMPM only (update domain lengths)
            [V,D] = eig(F.'*F);                                                 % eigen values and vectors F'F
            U     = V*sqrt(D)*V.';                                              % material stretch matrix
        end
        if norm(U)>5.00
            FHS_LS(end+1)=mp;
        end
    end
end

if ~isempty(FHS_LS)
    % 更新位置
    for i=1:length(FHS_LS)
        mp=FHS_LS(i);
        nIN = mpData(mp).nIN;                                                   % nodes associated with material point
        nn  = length(nIN);                                                      % number nodes
        N   = mpData(mp).Svp;                                                   % basis functions
        ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
        mpU = N*uvw(ed);                                                        % material point displacement
        mpData(mp).u     = mpData(mp).u + mpU.';                                % update material point displacements
    end

    disp('  ==>>> 部分物质点速度/加速度/变形过大，已识别异常物质点，已适当调整物质点位置')
else
    % 更新位置
    for mp=1:nmp
        nIN = mpData(mp).nIN;                                                   % nodes associated with material point
        nn  = length(nIN);                                                      % number nodes
        N   = mpData(mp).Svp;                                                   % basis functions
        ed  = repmat((nIN.'-1)*nD,1,nD)+repmat((1:nD),nn,1);                    % nodal degrees of freedom
        mpU = N*uvw(ed);                                                        % material point displacement
        mpData(mp).u     = mpData(mp).u + mpU.';                                % update material point displacements
    end

    disp('  ==>>> 部分物质点速度/加速度/变形过大，未识别异常物质点，已调整全部物质点位置')
end

[mesh,mpData,VA_old,bounday_imfo] = elemMPinfo(mesh,mpData,Cal_par);
NRit = 0;
[nodes,nD] = size(mesh.coord);
nDoF = nodes*nD;                                                            % total number of degrees of freedom
frct = zeros(nDoF,1);                                                       % zero the reaction forces
uvw  = zeros(nDoF,1);                                                       % zero the displacements
fd   = detFDoFs(mesh);                                                      % free degrees of freedom

end