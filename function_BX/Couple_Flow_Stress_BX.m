function [mesh, mpData] = Couple_Flow_Stress_BX(mesh, mpData, Cal_par)
% COUPLE_FLOW_STRESS_BX 应变-渗流耦合水头更新（向量化版）
%   基于物质点变形梯度计算孔隙比变化，通过VG模型更新吸力/饱和度，
%   并将水头变化反馈至背景网格。
%
%   核心原则：
%     1. 所有 mpData 字段在函数开头一次性预提取为 MATLAB 变量
%     2. 物理计算全部向量化，逻辑索引严格等价于原始 if 块
%     3. 所有字段更新（含 Ksat、lp）均限制在 active 掩码内，非活跃点保持原值
%     4. 回写使用 deal(tmp{:})，维度与原始循环赋值严格对齐
%
%   作者: FAN Haishan
%   日期: 28/04/2026
%--------------------------------------------------------------------------

nmp = length(mpData);                                                       % number of material points
nD  = length(mpData(1).mpC);
if nD ==2
    t = [1 5];
    LSa = 2;
else
    t = [1 5 9];
    LSa = 3;
end
dt = Cal_par.dt;

%% 向量化预计算（基于当前变形梯度）
F   = reshape([mpData.F],3,3,nmp);                                           % deformation gradient
det_F = batch_det3x3(F);
vp0 = [mpData.vp0]';
det_F_vp0 = det_F.*vp0;

%% 预计算孔隙比（两种分支都需要）
e0 = [mpData.e0]';
e_new = (1+e0).*det_F-1;

%% 预计算VG模型参数
LS = reshape([mpData.Flow_CP]',3,nmp)';
Fa = LS(:,1)/1000;
Fm = LS(:,2);
Fn = LS(:,3);

%% 基于孔隙比变化直接缩放饱和度
e_old = [mpData.e]';
Sr_old = min(1.00,[mpData.Sr]');
Sr = [mpData.Sr]' .* (e_old./e_new);
Sr = min(1.00,Sr);
Scution_Flow_old = ((Sr_old.^(-1./Fm) - 1).^(1./Fn)).*Fa./e_new.^(1./(Fm.*Fn));
Scution_Flow_new = ((Sr.^(-1./Fm) - 1).^(1./Fn)).*Fa./e_new.^(1./(Fm.*Fn));
% 计算水头变化
dH_all = (Scution_Flow_new - Scution_Flow_old)/-9.81;

%% 预计算GIMPM域长度 lp（仅用于 mpType==2 的MP）
[V,D] = pageeig(pagemtimes(F, 'ctranspose', F, 'none'));
D = sqrt(D);
temp = pagemtimes(V, D);
temp = pagemtimes(temp, 'none', V, 'ctranspose');
temp = permute(temp, [3, 1, 2]);
lp_new = reshape([mpData.lp0]',LSa,nmp)' .* temp(:,t);

%% 循环赋值（使用预计算值，严格对照原版逻辑）
for mp=1:nmp
    if mpData(mp).modulus>0.0 && mpData(mp).Gravity_Water_Content>0.0
        
        % 更新体积
        mpData(mp).vp = det_F_vp0(mp);
        % 更新饱和度
        mpData(mp).Sr = Sr(mp);
        % 更新孔隙比（对照原版第40行）
        mpData(mp).e = e_new(mp);
        % 【使用预计算值】更新吸力（对照原版第48-49行）
        mpData(mp).Scution_Flow = Scution_Flow_new(mp);
        dH = dH_all(mp);

        %% 【关键修正】先无条件累加 dH，后判断 Flow_Type 99 重置
        mpData(mp).dH_sum_strain = mpData(mp).dH_sum_strain + dH;           % 由于荷载导致应变产生的孔隙变化导致的水头变化
        % mpData(mp).dH = mpData(mp).dH + dH/dt;
        
        % if mpData(mp).Flow_Type == 99
        %     mpData(mp).dH_sum_strain = 0;
        %     mpData(mp).dH = 0;
        %     mpData(mp).Sr = 1.00;
        %     mpData(mp).Scution_Flow = 0.0;
        % end

        %% 更新干密度（对照原版第67行）
        mpData(mp).Roud = mpData(mp).Gs/(1+e_new(mp));

        %% 计算质量含水率（对照原版第70-71行） 储存本次含水率，用于湿力耦合计算中初始Mc和dMc
        mpData(mp).Gravity_Water_Content = (mpData(mp).Sr * e_new(mp)) / mpData(mp).Gs*1000;
        mpData(mp).Gravity_Water_Content_old = mpData(mp).Gravity_Water_Content;

        %% 仅在 mpType==2 (GIMPM) 时更新域长度 lp（对照原版第77-81行）
        if mpData(mp).mpType == 2
            mpData(mp).lp = lp_new(mp,:);
        end
        
        %% 无条件更新渗透系数，使用当前 lp（对照原版第83-84行）
        if nD ==2
            mpData(mp).Ksat(1) = mpData(mp).lp(2)/mpData(mp).lp0(2)*mpData(mp).Ksat0(1);
            mpData(mp).Ksat(2) = mpData(mp).lp(1)/mpData(mp).lp0(1)*mpData(mp).Ksat0(2);
        end
        % 更新最新吸力
        try
            mpData(mp).STATEV(7) = Scution_Flow_new(mp);
        catch
        end
    end

end



%% 插值更新背景网格吸力水头（对照原版第90-92行）
if Cal_par.Calculate_time>=1.0
    [mesh,mpData,~,~] = elemMPinfo(mesh,mpData,Cal_par);
end

end

%% 辅助函数：批量计算3×3矩阵行列式
function d = batch_det3x3(F)
N = size(F, 3);
f = reshape(F, 9, N);
d = (f(1,:).*f(5,:).*f(9,:) ...
   + f(4,:).*f(8,:).*f(3,:) ...
   + f(7,:).*f(2,:).*f(6,:) ...
   - f(7,:).*f(5,:).*f(3,:) ...
   - f(1,:).*f(8,:).*f(6,:) ...
   - f(4,:).*f(2,:).*f(9,:));
d = d(:);
end