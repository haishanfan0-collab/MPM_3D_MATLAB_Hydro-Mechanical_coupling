function [mpData] = updateMPs_BX(uvw,duvw,dduvw,mpData,mesh)
% 【核心-物质点状态更新代码】，内置变形梯度更新、体积/孔隙比计算、节点插值位移/速度/
% 加速度、GIMPM域长度演化、渗流吸力-饱和度-密度耦合更新全过程
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行物质点状态全面更新，涵盖基于变形梯度的体积与孔隙比计算、右拉伸张量谱分解更新
% GIMPM域长度、节点位移/速度/加速度/水头cumsum前缀和向量化插值、物质点坐标与历史变量
% 更新、渗流模块吸力-饱和度SWRC重新计算、干密度-质量含水率-物质点质量耦合更新全过程，
% 为物质点法每步迭代后的核心更新子函数。
%
%--------------------------------------------------------------------------
% [mpData] = updateMPs_BX(uvw, duvw, dduvw, mpData, mesh)
%--------------------------------------------------------------------------
% 输入:
% uvw    - 背景网格节点位移增量向量 [nDoF×1]
% duvw   - 背景网格节点速度向量 [nDoF×1]
% dduvw  - 背景网格节点加速度向量 [nDoF×1]
% mpData - 物质点信息结构体数组，含 mpC（坐标）、F（变形梯度）、Fn（历史变形梯度）、
%          vp0（初始体积）、e0（初始孔隙比）、lp0（初始域长）、epsE/epsP（弹塑性应变）、
%          sig（应力）、u/du/ddu（位移/速度/加速度）、Gs（土粒比重）、Flow_CP（VG参数）、
%          H0（初始水头）、dH_sum（水头累积）、dH_sum_strain（应变致水头变化）、
%          Gravity_Water_Content（含水率）、mpM（质量）等
% mesh   - 背景网格结构体，含 H（节点总水头）、HVA.dHdt（节点水头变化率）
%--------------------------------------------------------------------------
% 输出:
% mpData - 更新后的物质点结构体数组，含更新后的 mpC（坐标）、vp（体积）、e（孔隙比）、
%          lp（域长度）、u/du/ddu（位移/速度/加速度）、epsEn/epsPn/Fn（历史变量）、
%          sign（历史应力）、Sr（饱和度）、Scution_Flow（吸力）、Roud（干密度）、
%          Gravity_Water_Content（含水率）、mpM（质量）、Saturated_YES（饱和标记）等
%--------------------------------------------------------------------------
% 关键算法:
%   1. 变形梯度: F = dF * Fn，det_F = det(F)，vp = det_F * vp0
%   2. 孔隙比: e = (1+e0)*det_F - 1
%   3. GIMPM域长: 右拉伸张量谱分解 U = V*sqrt(D)*V'，lp = lp0 .* diag(U)
%   4. 节点插值: cumsum前缀和提取位移/速度/加速度/水头/水头率（向量化无循环）
%   5. 坐标更新: mpC = mpC + du，u = u + du
%   6. 历史变量: epsEn=epsE, epsPn=epsP, Fn=F, sign=sig, epsn=eps
%   7. 吸力计算: Suction = -(H0 + dH_sum + dH_sum_strain) * 9.81
%   8. SWRC更新: Sr = min([1+(Suction/coeff_swrc)^n]^{-m}, 1.0)
%   9. 密度更新: Roud = Gs/(1+e)，含水率: w = Sr*e/Gs*1000
%  10. 质量更新: mpM = vp*(1+w)*Roud
%--------------------------------------------------------------------------
% 此子函数包含子函数:
% batch_det3x3 - 批量3×3矩阵行列式计算（萨吕法则向量化）
%--------------------------------------------------------------------------
nD  = length(mpData(1).mpC);                                                % number of dimensions
if nD ==2
    t = [1 5];
    LSa = 2;
else
    t = [1 5 9];
    LSa = 3;
end                                                                         % stretch components for domain updating
nmp = length(mpData);                                                       % number of material points
Name_mp = fieldnames(mpData);

F   = reshape([mpData.F],3,3,nmp);                                          % deformation gradient
det_F = batch_det3x3(F);
vp = det_F.*[mpData.vp0]';
e_new = (1+[mpData.e0]').*det_F-1;

[V,D] = pageeig(pagemtimes(F, 'ctranspose', F, 'none'));  % 注意：这里 V 是 3×3×N
D = sqrt(D);

% lp 多维矩阵一次成型
temp = pagemtimes(V, D);
temp = pagemtimes(temp, 'none', V, 'ctranspose');
temp = permute(temp, [3, 1, 2]);
lp = reshape([mpData.lp0]',LSa,nmp)' .* temp(:,t);

% 预提取渗流参数（常量部分，用于循环内计算）
Gs_all = [mpData.Gs]';        % 土粒比重
Flow_CP_all = reshape([mpData.Flow_CP]', 3, nmp)';  % [Fa, Fm, Fn]

%% 牛-一次性插值计算（去除for）
nIN_all = [mpData.nIN];
ed_all  = repmat((nIN_all-1)*nD,nD,1)+(1:nD).';

uvw_all = uvw(ed_all);
duvw_all = duvw(ed_all);
dduvw_all = dduvw(ed_all);
H_all = mesh.H(nIN_all)';
dH_all = mesh.HVA.dHdt(nIN_all)';

N_all = [mpData.Svp];
% G_scaled(:,j) = G_all(:,j) * values(j)
uvw_scaled = uvw_all .* N_all;  % 2 x 189796
duvw_scaled = duvw_all .* N_all;  % 2 x 189796
dduvw_scaled = dduvw_all .* N_all;  % 2 x 189796
H_scaled = H_all .* N_all;  % 1 x 189796
dH_scaled = dH_all .* N_all;  % 1 x 189796

% 列方向累积求和（前缀和）
Prefix_uvw = cumsum(uvw_scaled, 2);  % 2 x 189796
Prefix_duvw = cumsum(duvw_scaled, 2);  % 2 x 189796
Prefix_dduvw = cumsum(dduvw_scaled, 2);  % 2 x 189796
Prefix_H = cumsum(H_scaled, 2);  % 2 x 189796
Prefix_dH = cumsum(dH_scaled, 2);  % 2 x 189796

nn_sizeall = [mpData.nn_size]';
% 通过索引相减提取各段结果（向量化，无循环）
s = nn_sizeall(:,1);
e = nn_sizeall(:,2);
% 初始化结果矩阵 4 x 23530
Result_uvw = Prefix_uvw(:,e);
Result_duvw = Prefix_duvw(:,e);
Result_dduvw = Prefix_dduvw(:,e);
Result_H = Prefix_H(:,e);
Result_dH = Prefix_dH(:,e);
% 处理起始位置大于1的段：减去前一个前缀和
% 利用逻辑掩码处理 s=1 的边界情况

mask = (s > 1);
if any(mask)
    s_prev = s - 1;
    Result_uvw(:, mask) = Result_uvw(:, mask) - Prefix_uvw(:, s_prev(mask));
    Result_duvw(:, mask) = Result_duvw(:, mask) - Prefix_duvw(:, s_prev(mask));
    Result_dduvw(:, mask) = Result_dduvw(:, mask) - Prefix_dduvw(:, s_prev(mask));
    Result_H(:, mask) = Result_H(:, mask) - Prefix_H(:, s_prev(mask));
    Result_dH(:, mask) = Result_dH(:, mask) - Prefix_dH(:, s_prev(mask));
end

%% 仅简单计算+赋值
for mp=1:nmp
    mpU = Result_uvw(:,mp)';                                                % material point displacement
    dmpU = Result_duvw(:,mp)';                                              % material point displacement
    ddmpU = Result_dduvw(:,mp)';                                            % material point displacement

    mpData(mp).mpC   = mpData(mp).mpC + mpU;                                % update material point coordinates
    mpData(mp).vp    = vp(mp);                                              % update material point volumes
    mpData(mp).epsEn = mpData(mp).epsE;                                     % update material point elastic strains
    mpData(mp).epsPn = mpData(mp).epsP;                                     % update material point plastic strains
    mpData(mp).Fn    = mpData(mp).F;                                        % update material point deformation gradients
    mpData(mp).epsn = mpData(mp).eps;

    mpData(mp).u     = mpData(mp).u + mpU.';                                % update material point displacements
    mpData(mp).du     = dmpU';
    mpData(mp).ddu     = ddmpU';

    mpData(mp).sign = mpData(mp).sig;
    mpData(mp).e = e_new(mp);
    mpData(mp).lp = lp(mp,:);                                           % update domain lengths
    mpData(mp).lp = min(mpData(mp).lp,mpData(mp).lp0*5.00);

    try
        mpData(mp).STATEVn = mpData(mp).STATEV;
    catch
    end

    %% 渗流模块更新（修复Sr滞后问题）
    if sum(strcmp(Name_mp,'Flow_Type'))>0 && mpData(mp).Gravity_Water_Content>0.0
        mpData(mp).dH   = Result_dH(mp);
        mpData(mp).epsn_Flow = mpData(mp).eps_Flow;                         % update material point elastic strains
        mpData(mp).sign_Flow = mpData(mp).sig_Flow;                         % update material point plastic strains=

        % 更新水头
        mpData(mp).dH_sum    = mpData(mp).dH_sum + Result_H(mp);            % update material point displacements
        
        % 计算当前吸力（必须在此步骤后，确保使用最新值）
        Suction = (mpData(mp).H0 + mpData(mp).dH_sum + mpData(mp).dH_sum_strain)*-9.81;
        if Suction>=-1e-18
            mpData(mp).Scution_Flow = max(0.0,Suction);
            mpData(mp).Pore_Pressure = 0.0;
        else
            mpData(mp).Scution_Flow = 0.0;
            mpData(mp).Pore_Pressure = max(0.0,-Suction);
        end

        % if mpData(mp).Flow_Type == 99
        %     mpData(mp).Scution_Flow = 0.0;
        %     mpData(mp).dH_sum    = 0.0;
        %     mpData(mp).Sr = 1.00;
        %     mpData(mp).dH = 0;
        % else
            % 提取当前MP的SWRC参数
            Fa = Flow_CP_all(mp,1)/1000;
            Fm = Flow_CP_all(mp,2);
            Fn = Flow_CP_all(mp,3);
            
            % 计算SWRC系数
            coeff_swrc = Fa / (e_new(mp)^(1/(Fm*Fn)));
            
            % 基于当前吸力重新计算饱和度（关键修复：消除滞后）
            suction = mpData(mp).Scution_Flow;
            mpData(mp).Sr = min((1.0 + (suction / coeff_swrc)^Fn)^(-Fm),1.0);
            
            % 更新饱和状态标记
            mpData(mp).Saturated_YES = (mpData(mp).Sr >= 0.999);
        % end

        % 更新干密度（基于当前孔隙比）
        mpData(mp).Roud = Gs_all(mp)/(1+e_new(mp));

        % 计算质量含水率（基于更新后的Sr）
        mpData(mp).Gravity_Water_Content = (mpData(mp).Sr * e_new(mp)) / Gs_all(mp) * 1000;

        % 更新物质点质量（基于当前体积和更新后的含水率）
        mpData(mp).mpM = vp(mp)*(1+mpData(mp).Gravity_Water_Content)*mpData(mp).Roud;

    end

end
end

function d = batch_det3x3(F)
% BATCH_DET3X3 批量计算 3×3 矩阵行列式（萨吕法则向量化）
%   F: 3×3×N 数组
%   d: N×1 列向量（与 arrayfun(@det) 输出维度一致）
N = size(F, 3);

% 一次性 reshape 为 9×N，列优先顺序：
% [F11; F21; F31; F12; F22; F32; F13; F23; F33]
f = reshape(F, 9, N);

% 萨吕法则向量化计算（6 项展开）
d = (f(1,:).*f(5,:).*f(9,:) ...  % F11*F22*F33
   + f(4,:).*f(8,:).*f(3,:) ...  % F12*F23*F31
   + f(7,:).*f(2,:).*f(6,:) ...  % F13*F21*F32
   - f(7,:).*f(5,:).*f(3,:) ...  % F13*F22*F31
   - f(1,:).*f(8,:).*f(6,:) ...  % F11*F23*F32
   - f(4,:).*f(2,:).*f(9,:));     % F12*F21*F33

% 确保输出为列向量（与 arrayfun 行为一致）
d = d(:);
end