function makeVtkMP(mpC,sig,uvw,dudt,du2dt2,mCst,suction,WC,H_Change_Seepage,H_Change_Stress,mpFileName,plastic_strain,Add_s2)


[nmp,nD]=size(mpC);                                                         % number of material points and dimensions

fid=fopen(mpFileName,'wt');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'MATLAB generated vtk file, WMC\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %i double\n',nmp);

%% position output 
if nD<3
    mpC = [mpC zeros(nmp,3-nD)];  
end
fprintf(fid,'%f %f %f \n',mpC');
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %i\n',nmp);

%% stress output
fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,1)/1000);
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,2)/1000);
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,3)/1000);
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,4)/1000);
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_yz FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,5)/1000);
fprintf(fid,'\n');

fprintf(fid,'SCALARS sigma_zx FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',sig(:,6)/1000);
fprintf(fid,'\n');


%% displacement output
if nD==3
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,2));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_z FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,3));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS u_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',vecnorm(uvw, 2, 2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',dudt(:,1));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',dudt(:,2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_z FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',dudt(:,3));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',vecnorm(dudt, 2, 2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS du2dt2_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',vecnorm(du2dt2, 2, 2));
    fprintf(fid,'\n');
elseif nD==2
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,1));
    fprintf(fid,'\n');
    
    fprintf(fid,'SCALARS u_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw(:,2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS u_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',vecnorm(uvw, 2, 2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',dudt(:,1));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_y FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',dudt(:,2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',vecnorm(dudt, 2, 2));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS du2dt2_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',vecnorm(du2dt2, 2, 2));
    fprintf(fid,'\n');
elseif nD==1
    fprintf(fid,'SCALARS u_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',uvw);
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS dudt_x FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',dudt(:,1));
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS du2dt2_vec FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',du2dt2(:,1));
    fprintf(fid,'\n');
end

fprintf(fid,'SCALARS Modulus FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',mCst(:,1));
fprintf(fid,'\n');

% fprintf(fid,'SCALARS Yield_Stress FLOAT %i\n',1);
% fprintf(fid,'LOOKUP_TABLE default\n');
% fprintf(fid,'%f\n',mCst(:,3)/1e3);
% fprintf(fid,'\n');

if sum(suction)~=0
    fprintf(fid,'SCALARS Suction FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',suction);
    fprintf(fid,'\n');
end

if sum(WC)~=0
    fprintf(fid,'SCALARS Water_Content FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',WC);
    fprintf(fid,'\n');
end
% H_Change_Seepage,H_Change_Stress

if sum(H_Change_Seepage)~=0
    fprintf(fid,'SCALARS H_Change_Seepage FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',H_Change_Seepage);
    fprintf(fid,'\n');
end

if sum(H_Change_Stress)~=0
    fprintf(fid,'SCALARS H_Change_Stress FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%f\n',H_Change_Stress);
    fprintf(fid,'\n');
end

Plastic_TZ1 = ( 200*(mean(abs(plastic_strain(:,1)))) < abs(plastic_strain(:,1)) );
Plastic_TZ2 = ( 200*(mean(abs(plastic_strain(:,2)))) < abs(plastic_strain(:,2)) );
if any(Plastic_TZ1) || any(Plastic_TZ2)
    plastic_strain(Plastic_TZ1 | Plastic_TZ2,1) = 0.0;
    plastic_strain(Plastic_TZ1 | Plastic_TZ2,2) = 0.0;
end

% plastic_strain(mpC(:,2)>8*1.5-0.78,:) = 0.0;

fprintf(fid,'SCALARS Plastic_Shear FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',plastic_strain(:,1));
fprintf(fid,'\n');

fprintf(fid,'SCALARS Plastic_Volume FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',plastic_strain(:,2));
fprintf(fid,'\n');

fprintf(fid,'SCALARS Add_s2 FLOAT %i\n',1);
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%f\n',Add_s2);
fprintf(fid,'\n');

fclose('all');
end