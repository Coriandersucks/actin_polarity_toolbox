%% categorize actin to membrane orientation
% Load filaments and BasePath
BasePath = {};

for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step7/filaments_160_tomo_' num2str(t) '.mat'];
end
mkdir ./mapping3d/filaments_step7/
cd ./mapping3d/filaments_step7/

% define pixel

pixelR = 0.220651*4;% conver to nm
sphereR = 40/pixelR;


for k = 1:size(BasePath,2)
    % load filament structure
    load(FilamentStruct{k},'filament_struct_expanded');
    
    % load membrane and fit
    load([BasePath{k} '/cor/particle_list_manual_seg_mem_rot.mat'], 'mPlist3');
    [sf,gof] = fit([mPlist3(:,1),mPlist3(:,3)],mPlist3(:,2),'poly55');
    
    % load arp23
    try
        load([BasePath{k} '/cor/particle_list_manual_seg_arp_rot_D.mat'], 'alg_inv');
        arpplist = alg_inv.pick_position;
        
        arp23_selected = 1;
    catch
        arp23_selected = 0;
    end
    
    %%% categorze filaments goining to membrane==1; 
    %%%                     going away from membrane ==-1 
    %%%                     parallel to the membrane == 2
    %%%                     not processd == 0
    %%%
    %%% write normal vector and angle difference in stucture
    %%% psi+ to right psi- to left after rotation
    
    for i = 1:size(filament_struct_expanded,2)
        for j = 1: size(filament_struct_expanded(i).cor_filament_ext_psi_rot,1)
            mempoint = filament_struct_expanded(i).closest_mempoint(j,:);
            [fx,fz] = differentiate(sf,mempoint(1),mempoint(3));
            normalvector = -([fx -1 fz]./norm([fx -1 fz]));
            filament_struct_expanded(i).mem_NormV(j,:) = normalvector;
            if filament_struct_expanded(i).filament_processed_vec(j) ~= 0
                if filament_struct_expanded(i).filament_spin_cat(j)== filament_struct_expanded(i).filament_spin_indx
                    pointrot = filament_struct_expanded(i).rot_xyz_filament_ext_rot(j,:);
                    pointpsi = filament_struct_expanded(i).psi_filament_ext_rot(j);
                elseif filament_struct_expanded(i).filament_spin_cat(j)== -filament_struct_expanded(i).filament_spin_indx
                    pointrot = -filament_struct_expanded(i).rot_xyz_filament_ext_rot(j,:);
                    pointpsi = wrapTo180(filament_struct_expanded(i).psi_filament_ext_rot(j)-180);
                end
                crossdeg = acosd(dot(pointrot,normalvector)/(norm(pointrot)*norm(normalvector)));
                if pointpsi >= 0 
                    filament_struct_expanded(i).point2membrane_angle(j) = crossdeg;
                elseif pointpsi < 0 
                    filament_struct_expanded(i).point2membrane_angle(j) = -crossdeg;
                end
            else
                filament_struct_expanded(i).point2membrane_angle(j) = 0;
            end
        end
    end
    
    % compare angle difference by a range of membrane average normal angle
    % by sphere range
    
    for i = 1:size(filament_struct_expanded,2)
        for j = 1: size(filament_struct_expanded(i).cor_filament_ext_psi_rot,1)
            mempoint = filament_struct_expanded(i).closest_mempoint(j,:);
            memdistance = pdist2(mempoint,mPlist3);
            mempointrange = mPlist3((memdistance<= sphereR),:);
            [sf,gof] = fit([mempointrange(:,1),mempointrange(:,3)],mempointrange(:,2),'poly55');
            [fx,fz] = differentiate(sf,mempointrange(:,1),mempointrange(:,3));
            normalvector = [];
            normalvector = -([fx -ones(size(fx,1),1) fz]./vecnorm([fx -ones(size(fx,1),1) fz],2,2));
            mean_normalvector = mean(normalvector);
            filament_struct_expanded(i).mem_NormavgV(j,:) = mean_normalvector;
            if filament_struct_expanded(i).filament_processed_vec(j) ~= 0
                if filament_struct_expanded(i).filament_spin_cat(j)== filament_struct_expanded(i).filament_spin_indx
                    pointrot = filament_struct_expanded(i).rot_xyz_filament_ext_rot(j,:);
                    pointpsi = filament_struct_expanded(i).psi_filament_ext_rot(j);
                elseif filament_struct_expanded(i).filament_spin_cat(j)== -filament_struct_expanded(i).filament_spin_indx
                    pointrot = -filament_struct_expanded(i).rot_xyz_filament_ext_rot(j,:);
                    pointpsi = wrapTo180(filament_struct_expanded(i).psi_filament_ext_rot(j)-180);
                end
                crossdeg = acosd(dot(pointrot,mean_normalvector)/(norm(pointrot)*norm(mean_normalvector)));
                if pointpsi >= 0 
                    filament_struct_expanded(i).point2sphere_membrane_angle(j) = crossdeg;
                elseif pointpsi < 0 
                    filament_struct_expanded(i).point2sphere_membrane_angle(j) = -crossdeg;
                end
            else
                filament_struct_expanded(i).point2sphere_membrane_angle(j) = 0;
            end
        end
    end
    
    % Define filament category by 80 and 100 degree 
    for i = 1:size(filament_struct_expanded,2)
        if ~isempty(filament_struct_expanded(i).filament_spin_indx)
            for j = 1: size(filament_struct_expanded(i).psi_filament_ext_rot,1)
                crossdeg = abs(filament_struct_expanded(i).point2membrane_angle(j));
                if filament_struct_expanded(i).filament_spin_cat(j) ~= 0
                    
                    if crossdeg > 80 && crossdeg <100
                        filament_struct_expanded(i).filament_move_cat(j) = 2;
                    elseif crossdeg <= 80
                        filament_struct_expanded(i).filament_move_cat(j) = 1;
                    elseif crossdeg >= 100
                        filament_struct_expanded(i).filament_move_cat(j) = -1;
                        
                    end
                    
                else
                    filament_struct_expanded(i).filament_move_cat(j) = 0;
                end
                
            end
        end
    end
    
    % define barbed-end and average orientation
    for i = 1:size(filament_struct_expanded,2)
        if ~isempty(filament_struct_expanded(i).filament_spin_indx)
            Fvector = filament_struct_expanded(i).cor_filament_ext_psi_rot(end,:)-filament_struct_expanded(i).cor_filament_ext_psi_rot(1,:);
            filamentrealindex = find(filament_struct_expanded(i).filament_spin_cat~=0);
            filament_struct_expanded(i).filament_avg_point2memangle = rad2deg(atan2(mean(sind((filament_struct_expanded(i).point2sphere_membrane_angle(filamentrealindex)))),...
                mean(cosd((filament_struct_expanded(i).point2sphere_membrane_angle(filamentrealindex))))));
            if abs(filament_struct_expanded(i).filament_avg_point2memangle) > 80 && abs(filament_struct_expanded(i).filament_avg_point2memangle) <100
                filament_struct_expanded(i).filament_avg_move_cat= 2;
            elseif abs(filament_struct_expanded(i).filament_avg_point2memangle) <= 80
                 filament_struct_expanded(i).filament_avg_move_cat= 1;
            elseif abs(filament_struct_expanded(i).filament_avg_point2memangle) >= 100
                 filament_struct_expanded(i).filament_avg_move_cat= -1;
            end
            if filament_struct_expanded(i).filament_spin_cat(filamentrealindex(1)) == filament_struct_expanded(i).filament_spin_indx
                pointrot  = filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(1),:);
                crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                if crossdeg > 90
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                    
                else
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                    
                end
            elseif filament_struct_expanded(i).filament_spin_cat(filamentrealindex(1)) == -filament_struct_expanded(i).filament_spin_indx
                pointrot  = -filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(1),:);
                crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                if crossdeg > 90
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                else
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                end
            elseif filament_struct_expanded(i).filament_spin_cat(filamentrealindex(end)) == filament_struct_expanded(i).filament_spin_indx
                pointrot  = filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(end),:);
                crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                if crossdeg > 90
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                    
                else
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                end
            elseif filament_struct_expanded(i).filament_spin_cat(filamentrealindex(end)) == -filament_struct_expanded(i).filament_spin_indx
                pointrot  = -filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(1),:);
                crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                if crossdeg > 90
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                else
                    filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                    filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                end
                
            end
            
        end
    end
    
    
    
    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k)
end
clear




