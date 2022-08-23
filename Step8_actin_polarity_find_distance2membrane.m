%% find the shortest distance to the membrane of each segs


BasePath = {};

% load filaments
for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step5/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir ./mapping3d/filaments_step6

for k =1:size(BasePath,2)
    
    % loadmembrane and filaments
    load([BasePath{k} '/cor/particle_list_manual_seg_mem_rot.mat'], 'mPlist3');
    load(FilamentStruct{k},'filament_struct_expanded');
    
    % find distance
    for i = 1:size(filament_struct_expanded,2)
        minD = {};
        pointtomem = [];
        mempoint = [];
        filamentpoints = filament_struct_expanded(i).cor_filament_ext_psi_rot(:,1:3);
        for j = 1: size(filamentpoints,1)

            
            pointtomem = pdist2(filamentpoints(j,1:3),mPlist3);
            minD{1,j} = min(pointtomem);
            mempoint(j,:) = mPlist3(find(minD{1,j}==pointtomem,1),1:3);
            minD{2,j} = mempoint(j,:);

            filament_struct_expanded(i).closest_mempoint = mempoint;
            filament_struct_expanded(i).closest_mempoint_D = [minD{1,:}];
        



         end

     end

    
    

    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k)

    
end
    

clear

