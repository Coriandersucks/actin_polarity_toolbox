%% rotate filaments and membrane to psi and base

%Load filaments and BasePath
BasePath = {};

for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step4/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step5/filaments_160_tomo_' num2str(t) '.mat'];
end


mkdir ./mapping3d/filaments_step5/

% Rotate filaments arp23 and membrane

for k = 1:size(BasePath,2)
    
    % load filament structure and membrane and base
    load(FilamentStruct{k},'filament_struct_expanded');
    load([BasePath{k} '/cor/particle_list_manual_seg_mem.mat'], 'mPlist');
    load([BasePath{k} '/cor/particle_list_manual_seg_base.mat'], 'bPlist');
    % load Arp23
    try
        load([BasePath{k} '/cor/particle_list_manual_seg_arp.mat'], 'alg_inv');
        arpplist = alg_inv.pick_position;
        DV = alg_inv.DaughterV;
        MV = alg_inv.MotherV;
        
        arp23_selected = 1;
    catch
        arp23_selected = 0;
    end
    
    tomo_rot_psi = filament_struct_expanded(1).tomogram_psi_corr;
    tomo_rot_b = filament_struct_expanded(1).tomogram_base_corr;
    
    % rotate and move base to xy = 0
    corbPlist = [];
    corbPlist2 = [];
    corbPlist3 = [];
    for i = 1:size(bPlist,1)
        
        corbPlist(i,:) = actin_RotatePoints(bPlist(i,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),1024,1024,1024);
        corbPlist2(i,:) = RotatePoints(corbPlist(i,:),0,0,-tomo_rot_psi,1024,1024,1024);
        
        
    end
    
    tomo_z_shift = floor(mean(corbPlist3(:,3)));
    
    
    % rotate filaments
    for i = 1:size(filament_struct_expanded,2)
        cor = filament_struct_expanded(i).cor_filament_ext;
        psi_angle = filament_struct_expanded(i).psi_filament_ext;
        rot_xyz = filament_struct_expanded(i).rot_xyz_filament_ext;
        cor_rot = [];
        cor_rot2 = [];
        cor_rot3 = [];
        psi_angle_rot = [];
        rot_xyz_rot = [];
        rot_xyz_rot2 = [];
        rot_xyz_rot3 = [];
        for j = 1:size(cor,1)
            cor_rot(j,:) = actin_RotatePoints(cor(j,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),1024,1024,1024); %in rad
            cor_rot2(j,:) = RotatePoints(cor_rot(j,:),0,0,-tomo_rot_psi,1024,1024,1024); % in degree
            if ~ismember(k,[13  14  15  16])
                cor_rot3(j,:) = actin_RotatePoints(cor_rot2(j,:),0,pi,0,1024,1024,1024);
            else
                cor_rot3(j,:) = cor_rot2(j,:);
            end
            
            if filament_struct_expanded(i).filament_processed_vec(j) ==1
                rot_xyz_rot(j,:) = actin_RotatePoints(rot_xyz(j,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),0,0,0);
                rot_xyz_rot2(j,:) = RotatePoints(rot_xyz_rot(j,:),0,0,-tomo_rot_psi,0,0,0);
                psi_angle_rot(j,1) = psi_angle(j,1)-rad2deg(tomo_rot_b(1,3))-tomo_rot_psi;
                if ~ismember(k,[13  14  15  16])
                    psi_angle_rot(j,1) = -(psi_angle_rot(j,1)-180);
                    rot_xyz_rot3(j,:) = actin_RotatePoints(rot_xyz_rot2(j,:),0,pi,0,0,0,0);
                else
                    rot_xyz_rot3(j,:) = rot_xyz_rot2(j,:);
                end
                
            else
                psi_angle_rot(j,1) = 0;
                rot_xyz_rot3(j,:) = rot_xyz(j,:);
            end
            
        end
        cor_rot3(:,3) = cor_rot3(:,3)-tomo_z_shift;
        filament_struct_expanded(i).cor_filament_ext_psi_rot = cor_rot3;
        filament_struct_expanded(i).psi_filament_ext_rot = wrapTo180(psi_angle_rot);
        filament_struct_expanded(i).rot_xyz_filament_ext_rot = -rot_xyz_rot3 ;
    end
    
    
    % rotate membrane
    mPlist2 = [];
    mPlist3 = [];
    parfor i = 1:size(mPlist,1)
        mPlist(i,:) = actin_RotatePoints(mPlist(i,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),1024,1024,1024);
        mPlist2(i,:) = RotatePoints(mPlist(i,:),0,0,-tomo_rot_psi,1024,1024,1024);
        if ~ismember(k,[13  14  15  16])
            mPlist3(i,:) = actin_RotatePoints(mPlist2(i,:),0,pi,0,1024,1024,1024);
        else
            mPlist3(i,:) = mPlist2(i,:);
        end
        
    end
    mPlist3(:,3) = mPlist3(:,3)-tomo_z_shift;
    save([BasePath{k} '/cor/particle_list_manual_seg_mem_rot.mat'],'mPlist3');
    disp(k)
    
    
    % find the minZ in membrane and max Z of tomogram
    up = [];
    down = [];
    for i = 1:size(filament_struct_expanded,2)
        up(i,:) = max(filament_struct_expanded(i).cor_filament_ext_psi_rot(:,3));
    end
    down = min(mPlist3(:,3));
    upndown2 = [max(up) down max(up)-down];
    for i = 1:size(filament_struct_expanded,2)
        zzz = filament_struct_expanded(i).cor_filament_ext_psi_rot(:,3);
        zzzR = (zzz-upndown2(2))./upndown2(3);
        filament_struct_expanded(i).Z_thick = upndown2(3);
        filament_struct_expanded(i).Z_Ration = zzzR;
    end
    
    save(FilamentStruct2{k},'filament_struct_expanded');
    
    disp(k)
    
    % rotate Arp23
    if arp23_selected == 1
        arpplist2 = [];
        arpplist3 = [];
        DV2 = [];
        DV3 = [];
        MV2 = [];
        MV3 = [];
        for i = 1:size(arpplist,1)
            arpplist(i,:) = actin_RotatePoints(arpplist(i,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),1024,1024,1024); %in rad
            arpplist2(i,:) = RotatePoints(arpplist(i,:),0,0,-tomo_rot_psi,1024,1024,1024);
            DV(i,:) = actin_RotatePoints(DV(i,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),1024,1024,1024);
            DV2(i,:) = RotatePoints(DV(i,:),0,0,-tomo_rot_psi,1024,1024,1024);
            MV(i,:) = actin_RotatePoints(MV(i,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),1024,1024,1024);
            MV2(i,:) = RotatePoints(MV(i,:),0,0,-tomo_rot_psi,1024,1024,1024);
            if ~ismember(k,[13  14  15  16])
                arpplist3(i,:) = actin_RotatePoints(arpplist2(i,:),0,pi,0,1024,1024,1024);
                DV3(i,:) = actin_RotatePoints(DV2(i,:),0,pi,0,1024,1024,1024);
                MV3(i,:) = actin_RotatePoints(MV2(i,:),0,pi,0,1024,1024,1024);
            else
                arpplist3(i,:) = arpplist2(i,:);
                DV3(i,:) = DV2(i,:);
                MV3(i,:) = MV2(i,:);
            end
        end
        arpplist3(:,3) = arpplist3(:,3)-tomo_z_shift;
        DV3(:,3) = DV3(:,3)-tomo_z_shift;
        MV3(:,3) = MV3(:,3)-tomo_z_shift;
        alg_inv.pick_position_rot = arpplist3;
        alg_inv.DaughterV_rot = DV3;
        alg_inv.MotherV_rot = MV3;
        save([BasePath{k} '/cor/particle_list_manual_seg_arp_rot.mat'],'alg_inv');
        
    end
    disp(k)
    
    
end

clear
    







