%% find the shortest distance to the segs and membrane of each arp23 


BasePath = {};

%define pixel and membrane range
pixelR = 0.220651*4;
sphereR = 50/pixelR;

% load filaments
for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end

for k =1:size(BasePath,2)
    
    % load arp23 and filaments
    try
        load([BasePath{k} '/cor/particle_list_manual_seg_arp_rot.mat'], 'alg_inv');
        arpplist = alg_inv.pick_position_rot;
        
        arp23_selected = 1;
    catch
        arp23_selected = 0;
    end
    
    
    
    % find distance and closest filament 
    if arp23_selected == 1
        load([BasePath{k} '/cor/particle_list_manual_seg_mem_rot.mat'], 'mPlist3');
        load(FilamentStruct{k},'filament_struct_expanded');
        counter = 1;
        filamentpoint = [];
        mempoint = [];
        for i = 1:size(filament_struct_expanded,2)
            for j = 1:size(filament_struct_expanded(i).cor_filament_ext)
                filamentpoint{counter,1} = filament_struct_expanded(i).cor_filament_ext_psi_rot(j,:);
                filamentpoint{counter,2} = filament_struct_expanded(i).rot_xyz_filament_ext_rot(j,:);
                filamentpoint{counter,3} = filament_struct_expanded(i).filament_spin_cat(j);
                filamentpoint{counter,4} = filament_struct_expanded(i).filament_spin_indx;
                counter = counter + 1;
            end
        end
        for i = 1:size(arpplist,1)
            minD = {};
            pointtoarp23 = [];
            arp23point = arpplist(i,1:3);
            arp23V = (alg_inv.DaughterV_rot(i,:)+alg_inv.MotherV_rot(i,:))/2 -arp23point;
            pointtoarp23 = pdist2(arp23point,cell2mat(filamentpoint(:,1)));
            pointtomem = pdist2(arp23point,mPlist3);
            
            minD{1} = min(pointtoarp23);
            minDm{1} = min(pointtomem);
            filamentpp = cell2mat(filamentpoint(find(minD{1}==pointtoarp23,1),1));
            
            mempoint = mPlist3(find(minDm{1}==pointtomem,1),1:3);
            memdistance = pdist2(mempoint,mPlist3);
            mempointrange = mPlist3((memdistance<= sphereR),:);
            [sf,gof] = fit([mempointrange(:,1),mempointrange(:,3)],mempointrange(:,2),'poly55');
            
            minD{2} = filamentpp;
            if filamentpoint{find(minD{1}==pointtoarp23,1),3} ~= 0
                if filamentpoint{i,3} == filamentpoint{i,4}
                    crossdeg_arp2filament = acosd(dot(filamentpoint{find(minD{1}==pointtoarp23,1),2},arp23V)/(norm(filamentpoint{find(minD{1}==pointtoarp23,1),2})*norm(arp23V)));
                    minD{3} = crossdeg_arp2filament;
                elseif filamentpoint{i,3} == -filamentpoint{i,4}
                    crossdeg_arp2filament = acosd(dot(-filamentpoint{find(minD{1}==pointtoarp23,1),2},arp23V)/(norm(-filamentpoint{find(minD{1}==pointtoarp23,1),2})*norm(arp23V)));
                    minD{3} = crossdeg_arp2filament;
                else
                    minD{3} = 1000;
                end
            else
                minD{3} = 1000;
            end
            minDm{2} = mempoint;

            alg_inv.closest_filamentpoint(i,:) = minD{2};
            alg_inv.closest_filament_D(i,:)= minD{1};
            alg_inv.closest_mempoint(i,:) = minDm{2};
            alg_inv.arp2filament_angle(i,:) = minD{3};
            [fx,fz] = differentiate(sf,mempoint(1),mempoint(3));
            normalvector = -([fx -1 fz]./norm([fx -1 fz]));
            alg_inv.closest_mem_normV(i,:)= normalvector;
            normalvector = [];
            [fx,fz] = differentiate(sf,mempointrange(:,1),mempointrange(:,3));
            normalvector = -([fx -ones(size(fx,1),1) fz]./vecnorm([fx -ones(size(fx,1),1) fz],2,2));
            mean_normalvector = mean(normalvector);
            alg_inv.closest__mean_mem_normV(i,:)= mean_normalvector;
            alg_inv.closest_mem_D(i,:)= minDm{1};
            dau_rot = alg_inv.DaughterV_rot(i,:)-alg_inv.pick_position_rot(i,:);
            momrot = alg_inv.MotherV_rot(i,:)-alg_inv.pick_position_rot(i,:);
            VVV = (dau_rot+momrot)/2;
            crossdeg_dau = acosd(dot(dau_rot,mean_normalvector)/(norm(dau_rot)*norm(mean_normalvector)));
            crossdeg_mom = acosd(dot(momrot,mean_normalvector)/(norm(momrot)*norm(mean_normalvector)));
            crossdeg_VVV = acosd(dot(VVV,mean_normalvector)/(norm(VVV)*norm(mean_normalvector)));
            alg_inv.dau2membrane_angle(i,1) = crossdeg_dau;
            alg_inv.mom2membrane_angle(i,1) = crossdeg_mom;
            alg_inv.arp2membrane_angle(i,1) = crossdeg_VVV;

        end
        
        
        
    
    
    

    save([BasePath{k} '/cor/particle_list_manual_seg_arp_rot_D.mat'],'alg_inv');
    disp(k)
    end

    
end
    

clear

