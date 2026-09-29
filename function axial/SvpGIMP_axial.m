function [S, G, T] = SvpGIMP_axial(rp, rv, h, lp)

% 轴对称 uGIMP 形函数 - 平滑修正版
% 关键改进：截断域（rp < 0.5*lp）内 T/G 从 0 二次方增长到 0.5
% 效果：轴心处 T=0（消除奇异性），边界处 T/G=0.5（连续过渡）

tol = 1e-12;
if rp < 0, rp = 0; end
if lp <= 0, lp = eps; end
if h <= 0, error('Grid spacing h must be positive'); end
if rp < tol, rp = tol; end

delta = rp - rv;
n = rv / h;           % n_i = r_i / h
lp_t = 2 * lp / h;    % 无量纲粒子半宽
xi = 2 * delta / h;   % 无量纲坐标

S = 0; G = 0; T = 0;

% 截断域判定（Nairn Sec. 4.2）
isTruncated = (rp < 0.5 * lp);

% =========================================================================
% 标准域（rp >= 0.5*lp）：Nairn 解析解
% =========================================================================
if ~isTruncated
    denom = 2*n + xi;
    if abs(denom) < tol, denom = tol; end
    factor = 2/h;
    
    if xi > -(2 + lp_t) && xi < -(2 - lp_t)          % Region 1
        S = (2 + lp_t + xi)^2/(8*lp_t) * (1 - (2*(1-lp_t) + xi)/(3*denom));
        dSdxi = (2 + lp_t + xi)/(4*lp_t) * (1 - (2 - lp_t + xi)/(2*denom));
        G = dSdxi * factor;
    elseif xi >= -(2 - lp_t) && xi < -lp_t           % Region 2
        S = (2 + xi)/2 + lp_t^2/(6*denom);
        G = 0.5 * factor;
    elseif xi >= -lp_t && xi < lp_t                  % Region 3 (Case C)
        S = ((4 - lp_t)*lp_t - xi^2)/(4*lp_t) + xi*(xi^2 - 3*lp_t^2)/(12*lp_t*denom);
        
        % 子积分 G
        r_left = rp - lp; r_right = rp + lp;
        gauss_xi = [-sqrt(3/5), 0, sqrt(3/5)];
        gauss_w = [5/9, 8/9, 5/9];
        
        G_sum = 0; w_sum = 0;
        for i = 1:3
            r_i = (r_left + r_right)/2 + gauss_xi(i) * lp;
            if r_i < tol, r_i = tol; end
            delta_i = r_i - rv;
            if abs(delta_i) < h
                dN = (delta_i > 0)*(-1/h) + (delta_i < 0)*(1/h);
                weight = r_i * gauss_w(i) * lp;
                G_sum = G_sum + dN * weight;
                w_sum = w_sum + weight;
            end
        end
        
        if w_sum > tol
            G_gauss = G_sum / w_sum;
            dSdxi_th = -xi/(2*lp_t) - (lp_t^2 - xi^2)/(4*lp_t*denom);
            G_th = dSdxi_th * factor;
            if abs(G_gauss) > abs(G_th)*1.2
                G = G_th;
            else
                G = G_gauss;
            end
        end
        if abs(xi) < tol, G = 0; end
        
    elseif xi >= lp_t && xi < (2 - lp_t)             % Region 4
        S = (2 - xi)/2 - lp_t^2/(6*denom);
        G = -0.5 * factor;
    elseif xi >= (2 - lp_t) && xi < (2 + lp_t)       % Region 5
        S = (2 + lp_t - xi)^2/(8*lp_t) * (1 + (2*(1-lp_t) - xi)/(3*denom));
        dSdxi = -(2 + lp_t - xi)/(4*lp_t) * (1 + (2 - lp_t - xi)/(2*denom));
        G = dSdxi * factor;
    else
        return;
    end
    
    % 近轴区（rp < lp）硬截断
    if rp < lp
        T = sign(G) * abs(G) * 0.5;
    else
        T = 2*S/(h*denom);
    end
    
    if abs(G) > 1.2/h, G = sign(G)*1.2/h; end

% =========================================================================
% 截断域（rp < 0.5*lp）：平滑 T 增长（关键改进）
% =========================================================================
else
    lp_eff = rp + lp;
    rp_bar = lp_eff/2;
    r_max = rp + lp;
    
    if lp_eff < tol, return; end
    
    if (-lp) < delta && delta <= (lp) && rv >= 0 && rv <= r_max
        % 左右积分计算 S 和 G（保持原算法）
        if rv <= h
            S_left = (1-rv/h)*rv^2/2 + rv^3/(3*h);
            G_left = rv^2/(2*h);
        else
            a = max(rv-h, 0); b = rv;
            S_left = (b^3-a^3)/(3*h) - (rv-h)*(b^2-a^2)/(2*h);
            G_left = (b^2-a^2)/(2*h);
        end
        
        if (rv + h) >= r_max
            a = rv; b = r_max;
            S_right = (rv+h)*(b^2-a^2)/(2*h) - (b^3-a^3)/(3*h);
            G_right = -(b^2-a^2)/(2*h);
        else
            a = rv; b = rv+h;
            S_right = (rv+h)*(b^2-a^2)/(2*h) - (b^3-a^3)/(3*h);
            G_right = -(b^2-a^2)/(2*h);
        end
        
        denom_S = lp_eff * rp_bar;
        denom_G = 2 * lp_eff * rp_bar;
        
        if abs(denom_S) > tol && abs(denom_G) > tol
            S = (S_left + S_right) / denom_S;
            G = (G_left + G_right) / denom_G;
            
            % 关键改进：截断域内 T/G 二次方增长（0 → 0.5）
            % 当 rp=0 时，T=0（消除奇异性）
            % 当 rp=0.5*lp 时，T/G=0.5（与标准域连续）
            ratio = 0.5 * (2*rp/lp)^2;  % 二次方增长：0 到 0.5
            if abs(G) > tol
                T = sign(G) * abs(G) * ratio;
            else
                T = 0;
            end
        end
    end
end

if isnan(S) || isnan(G) || isnan(T), S = 0; G = 0; T = 0; end

end