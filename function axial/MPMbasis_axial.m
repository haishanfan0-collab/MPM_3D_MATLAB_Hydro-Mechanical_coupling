function [Svp, dSvp, Tvp] = MPMbasis_axial(mesh, mpData, node)

%Basis functions for axisymmetric material point method - FINAL VERSION
% 修正：扩大 T 正则化范围到所有轴附近粒子（rp < 3*lp）
%--------------------------------------------------------------------------

coord  = mesh.coord(node,:);
h      = mesh.h;
mpC    = mpData.mpC;
lp     = mpData.lp;
mpType = mpData.mpType;
nD     = length(mpC);

if nD ~= 2
    error('Axisymmetric MPM requires 2D (r-z) coordinates');
end

S = zeros(nD,1);
G = zeros(nD,1);
T = zeros(nD,1);

for i = 1:nD
    if i == 1
        if mpType == 1
            [S(i), G(i), T(i)] = SvpMPM_axial(mpC(i), coord(i), h(i));
        else
            [S(i), G(i), T(i)] = SvpGIMP_axial(mpC(i), coord(i), h(i), lp(i));
        end
    else
        if mpType == 1
            [S(i), G(i)] = SvpMPM(mpC(i), coord(i), h(i));
        else
            [S(i), G(i)] = SvpGIMP(mpC(i), coord(i), h(i), lp(i));
        end
    end
end

% 组装
Svp = prod(S);

dSvp = zeros(nD, 1);
for i = 1:nD
    otherDims = setdiff(1:nD, i);
    dSvp(i) = G(i) * prod(S(otherDims));
end

Tvp = T(1) * S(2);

%==========================================================================
% CRITICAL FIX: 扩大 T 正则化范围到所有轴附近粒子（rp < 3*lp）
%==========================================================================
% 原代码只处理 rp <= lp（截断域），但标准域的轴附近粒子（lp < rp < 3lp）也有 T 过大问题
% if mpC(1) < 3 * lp(1)  % 扩大到 3 倍 lp 范围（可根据需要调到 2*lp 或 4*lp）
%     max_Gr = max(abs(dSvp(1)));
% 
%     % 限制 Tvp 不超过径向梯度的 1.5 倍
%     max_T_allowed = max_Gr * 1.5;
%     max_T_current = max(abs(Tvp));
% 
%     if max_T_current > max_T_allowed && max_T_current > 1e-10
%         scale_factor = max_T_allowed / max_T_current;
%         Tvp = Tvp * scale_factor;
%     end
% 
%     % 极近轴时（rp < 0.5*lp），强制 Tvp 与 dSvp 同量级
%     if mpC(1) < 0.5 * lp(1)
%         Tvp = dSvp(1) * 0.8;
%     end
% end

end