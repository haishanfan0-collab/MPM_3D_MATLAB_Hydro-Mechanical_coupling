function postPro(mpData,mesh,Cal_par,uvw,VA,SAVE_NAME,Plot_OFF)
% 【后处理与可视化输出代码】，内置单相/多相物质点 VTK 文件生成、背景网格 VTK 文件生成、
% 场变量自动提取与合并输出全过程
%--------------------------------------------------------------------------
% 作者: FAN Haishan
% 日期: 28/04/2026
% 描述:
% 进行物质点法计算结果后处理与可视化输出，支持单相单组物质点与多相多组物质点两种模式。
% 自动提取物质点应力、位移、速度、加速度、吸力、含水率、水头变化、塑性应变等场变量，
% 以及背景网格节点位移与速度场，调用 VTK 文件生成器输出为 Paraview 可视化格式。
%
%--------------------------------------------------------------------------
% postPro(mpData, mesh, Cal_par, uvw, VA, SAVE_NAME)
%--------------------------------------------------------------------------
% 输入:
% mpData    - 物质点信息结构体数组（单组模式）或结构体（多组模式，字段为各相名称），
%             含 sig, mpC, u, du, ddu, modulus, STATEV, Scution_Flow,
%             Gravity_Water_Content, dH_sum, dH_sum_strain 等
% mesh      - 背景网格结构体（单组模式）或含多组子结构体的结构体（多组模式），
%             含 coord（节点坐标）、etpl（单元拓扑）等
% Cal_par   - 计算参数结构体，含 lstp_plot（当前输出步号）
% uvw       - 节点位移向量 [nDoF×1]（单组模式）或含多组位移的结构体（多组模式）
% VA        - 节点速度信息结构体，含 dudt（单组模式）或含多组速度的子结构体（多组模式）
% SAVE_NAME - 输出文件夹名字符串，文件保存至 output - SAVE_NAME/ 目录下
%--------------------------------------------------------------------------
% 输出:
% 无显式返回值，副作用为在 output - SAVE_NAME/ 目录下生成：
%   mpData_XXX.vtk - 物质点 VTK 文件（XXX 为 Cal_par.lstp_plot 步号）
%   mesh_XXX.vtk   - 背景网格 VTK 文件
%--------------------------------------------------------------------------
% 此子函数调用子函数:
% makeVtkMP - 物质点 VTK 文件生成器，接收坐标、应力、位移、速度、吸力、含水率等场变量
% makeVtk   - 背景网格 VTK 文件生成器，接收节点坐标、单元拓扑、位移、速度等场变量
%--------------------------------------------------------------------------
if nargin<7
    Plot_OFF = true;
end
Name_mp = fieldnames(mpData);

if sum(strcmp(Name_mp,'mpType'))>0
    if Plot_OFF
        mpData([mpData.mpType]==3)=[];
    end
    [~,nD] = size(mesh.coord);
    nmp  = length(mpData);
    sig = reshape([mpData.sig],6,nmp)';                                         % all material point stresses (nmp,6)
    mpC = reshape([mpData.mpC],nD,nmp)';                                        % all material point coordinates (nmp,nD)
    mpU = [mpData.u]';                                                          % all material point displacements
    modulus = [mpData.modulus]';
    mp_dudt = [mpData.du]';                                                     % all material point displacements
    mp_du2dt2 = [mpData.ddu]';                                                  % all material point displacements
    STATEV = [mpData.STATEV]';
    plastic_strain = STATEV(:,9:10);
    Add_s2 = -(sig(:,2)/1000-STATEV(:,18)/1000);
    try
        suction = [mpData.Scution_Flow]';
        WC = [mpData.Gravity_Water_Content]';
        H_Change_Seepage = [mpData.dH_sum]';
        H_Change_Stress = [mpData.dH_sum_strain]';
    catch
        suction = zeros(nmp,1);
        WC = zeros(nmp,1);
        H_Change_Seepage = zeros(nmp,1);
        H_Change_Stress = zeros(nmp,1);
    end
    mpDataName = sprintf(['output - ',SAVE_NAME,'/mpData_%i.vtk'],Cal_par.lstp_plot);        % MP output data file name
    makeVtkMP(mpC,sig,mpU,mp_dudt,mp_du2dt2,modulus,suction,WC,H_Change_Seepage,H_Change_Stress,mpDataName,plastic_strain,Add_s2);                   % generate material point VTK file

    meshName = sprintf(['output - ',SAVE_NAME,'/mesh_%i.vtk'],Cal_par.lstp_plot);            % MP output data file name
    makeVtk(mesh.coord,mesh.etpl,uvw,VA.dudt,meshName);                            % generate mesh VTK file
else
    [~,nD] = size(mesh.(Name_mp{1}).coord);
    for i = 1:length(Name_mp)
        %mpData.(Name_mp{i})([mpData.mpType]==3)=[];
        %mpData.(Name_mp{i})([mpData.mpType]==3)=[];
        nmp  = length(mpData.(Name_mp{i}));
        if i ==1
            STATEV = [mpData.(Name_mp{i}).STATEV]';
            sig = reshape([mpData.(Name_mp{i}).sig],6,nmp)';                                         % all material point stresses (nmp,6)
            
            mpC = reshape([mpData.(Name_mp{i}).mpC],nD,nmp)';                                        % all material point coordinates (nmp,nD)
            mpU = [mpData.(Name_mp{i}).u]';                                                          % all material point displacements
            modulus = [mpData.(Name_mp{i}).modulus]';
            mp_dudt = [mpData.(Name_mp{i}).du]';                                                     % all material point displacements
            mp_du2dt2 = [mpData.(Name_mp{i}).ddu]';                                                  % all material point displacements
            UVW = uvw.(Name_mp{i});
            DUDT = VA.(Name_mp{i}).dudt;
            
            plastic_strain = STATEV(:,9:10);
            try
                suction = [mpData.(Name_mp{i}).Scution_Flow]';
                WC = [mpData.(Name_mp{i}).Gravity_Water_Content]';
                H_Change = [mpData.(Name_mp{i}).dH_sum]'+[mpData.(Name_mp{i}).dH_sum_strain]';
            catch
                suction = zeros(nmp,1);
                WC = zeros(nmp,1);
                H_Change = zeros(nmp,1);
            end
        else
            sig = [sig;reshape([mpData.(Name_mp{i}).sig],6,nmp)'];                                   % all material point stresses (nmp,6)
            mpC = [mpC;reshape([mpData.(Name_mp{i}).mpC],nD,nmp)'];                                  % all material point coordinates (nmp,nD)
            mpU = [mpU;[mpData.(Name_mp{i}).u]'];                                                    % all material point displacements
            modulus = [modulus;[mpData.(Name_mp{i}).modulus]'];
            mp_dudt = [mp_dudt;[mpData.(Name_mp{i}).du]'];                                           % all material point displacements
            mp_du2dt2 = [mp_du2dt2;[mpData.(Name_mp{i}).ddu]'];                                      % all material point displacements
            UVW = UVW + uvw.(Name_mp{i});
            DUDT = DUDT + VA.(Name_mp{i}).dudt;
            STATEV = [mpData.(Name_mp{i}).STATEV]';
            plastic_strain = [plastic_strain;STATEV(:,9:10)];

            try
                suction = [suction;[mpData.(Name_mp{i}).Scution_Flow]'];
                WC = [WC;mpData.(Name_mp{i}).Gravity_Water_Content]';
                H_Change = [H_Change;mpData.(Name_mp{i}).dH_sum + [mpData.(Name_mp{i}).dH_sum_strain]]';
            catch
                suction = [suction;zeros(nmp,1)];
                WC = [WC;zeros(nmp,1)];
                H_Change = [H_Change;zeros(nmp,1)]';
            end
        end
    end
    mpDataName = sprintf(['output - ',SAVE_NAME,'/mpData_%i.vtk'],Cal_par.lstp_plot);        % MP output data file name
    makeVtkMP(mpC,sig,mpU,mp_dudt,mp_du2dt2,modulus,suction,WC,H_Change,mpDataName,plastic_strain,-(sig(:,2)/1000-STATEV(:,18)/1000));                   % generate material point VTK file

    meshName = sprintf(['output - ',SAVE_NAME,'/mesh_%i.vtk'],Cal_par.lstp_plot);            % MP output data file name
    makeVtk(mesh.(Name_mp{1}).coord,mesh.(Name_mp{1}).etpl,UVW,DUDT,meshName);                            % generate mesh VTK file
end
end