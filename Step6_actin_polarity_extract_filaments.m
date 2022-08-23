%% Break filament assignment into smaller pieces to reduce segementation artifect

% Load data matrix
load('./mapping/actin_polarity_mapping_160_step_2_new_more_clear3.mat','actin_data_matrix');

% Define filament length min threshold
filament_length_min_thres = 30;

% Define filament length max threshold
filament_length_max_thres = 30;

% Process filament
for k=1:58

    % Select only those segments that were used in the respective tomogram and in Relion 3D refine
    indx_segments = find(actin_data_matrix(:,1) == k & actin_data_matrix(:,4) == 1);
    
    % Reduce actin data matrix
    actin_data_matrix_red = actin_data_matrix(indx_segments,:);
    
    % Extract filament assignment of each segment in the respective tomogram
    filament_vect = actin_data_matrix_red(:,38);
    
    % Initialize refined filament assignment
    filament_vect_ref = filament_vect;
    
    % Initalize number of segment chains longer than filament length max threshold
    num_of_filament_length_max_thres = 1;
    
    for i=1:max(filament_vect(:))
        
         % Extract filament indices
         indx_filament_ext = find(filament_vect==i);
         
         if size(indx_filament_ext,1) > filament_length_max_thres
              
              correction_constant = 0;
              
              for j=1:size(indx_filament_ext,1)
              
                   filament_vect_ref(indx_filament_ext(j,1)) = filament_vect_ref(indx_filament_ext(j,1)) + num_of_filament_length_max_thres .* max(filament_vect(:)) + correction_constant;
                   
                   if mod(j,filament_length_min_thres) == 0
                       correction_constant = correction_constant + 1;
                   end
                   
              end
              
              num_of_filament_length_max_thres = num_of_filament_length_max_thres + 1;
              
         end
        
    end
    
    % Replace refined filament assignment
    actin_data_matrix_red(:,38) = filament_vect_ref;
    actin_data_matrix(indx_segments,:) = actin_data_matrix_red;
    
    disp(k);
    
    clear indx_segments actin_data_matrix_red filament_vect filament_vect_ref num_of_filament_length_max_thres indx_filament_ext correction_constant;
    
end

% Save workspace
save('./mapping/actin_polarity_mapping_160_step_3_new.mat','actin_data_matrix');
clear all;


%% Build filaments structure 


% Load actin data matrix
load('./mapping/actin_polarity_mapping_160_step_3_new.mat','actin_data_matrix');

% Define filament length in terms of segments to analyze threshold
filament_length_thres = 3;

% Define filament processed fraction threshold
filament_processed_frac_thres = 0.5;

% construct data path
for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step0/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir ./mapping3d/filaments_step0

% Process filament
for k=1:58

    % Select filament from tomogram
    indx_segments = find(actin_data_matrix(:,1) == k);
    
    % Reduce actin data matrix
    actin_data_matrix_red = actin_data_matrix(indx_segments,:);
    
    % Extract filament assignment of each segment in the respective tomogram
    filament_vect = actin_data_matrix_red(:,38);
    
    % Extract filaments and filament distributions
    filament_struct = [];
    rot_distribution = [];
    ccc_distribution = [];
    tomo_distribution = [];
    segment_sensitivity_distribution = [];
    zaehler_f = 1;
    
    
    
    for i=1:max(filament_vect(:))
        
        % Extract filament indices
        indx_filament_ext = find(filament_vect==i);
        
        % Extract real_processed_filament indices
        
        indx_filament_ext_real = find(actin_data_matrix_red(indx_filament_ext,4) == 1);
                        
        % Measure filament length
        filament_length = size(indx_filament_ext,1);
        
        % Measure filament processed fraction
        filament_processed_fraction = size(indx_filament_ext_real,1)/filament_length;
        
        if ~isempty(indx_filament_ext) && size(indx_filament_ext_real,1) >= filament_length_thres && (filament_processed_fraction > filament_processed_frac_thres)
        
              filament_struct(zaehler_f).filament_indx = zaehler_f;
              
              % Extract coordinates of alignment corrected segments
              cor_filament_ext = [actin_data_matrix_red(indx_filament_ext,8)./4 actin_data_matrix_red(indx_filament_ext,9)./4 actin_data_matrix_red(indx_filament_ext,10)./4];
              cor_filament_delta = [actin_data_matrix_red(indx_filament_ext,14)./4 actin_data_matrix_red(indx_filament_ext,15)./4 actin_data_matrix_red(indx_filament_ext,16)./4];
              cor_filament_inv_delta = [actin_data_matrix_red(indx_filament_ext,21)./4 actin_data_matrix_red(indx_filament_ext,22)./4 actin_data_matrix_red(indx_filament_ext,23)./4];
              % Extract actin directions
              phi_filament_ext = actin_data_matrix_red(indx_filament_ext,18);
              psi_filament_ext = actin_data_matrix_red(indx_filament_ext,19);
              theta_filament_ext = actin_data_matrix_red(indx_filament_ext,20);
              
              % Extract actin directions
              rot_filament_ext = actin_data_matrix_red(indx_filament_ext,24:26);
              rot_distribution = [rot_distribution; rot_filament_ext];
              
              % Extract ccc distribution
              ccc_filament_ext = actin_data_matrix_red(indx_filament_ext,35);
              ccc_distribution = [ccc_distribution; ccc_filament_ext];
              
              % Extract tomo distribution
              tomo_filament_ext = actin_data_matrix_red(indx_filament_ext,1);
              tomo_distribution = [tomo_distribution; tomo_filament_ext];
              
              % Extract filament sensitivity distribution
              segment_sensitivity_filament_ext = actin_data_matrix_red(indx_filament_ext,37);
              segment_sensitivity_distribution = [segment_sensitivity_distribution; segment_sensitivity_filament_ext];
              filament_processed_vec = actin_data_matrix_red(indx_filament_ext,4);
              
              % Extract actin vector
              vect_filament_ext = [actin_data_matrix_red(indx_filament_ext,24) actin_data_matrix_red(indx_filament_ext,25) actin_data_matrix_red(indx_filament_ext,26)];
              
              filament_struct(zaehler_f).cor_filament_ext = cor_filament_ext;
              filament_struct(zaehler_f).cor_filament_delta = cor_filament_delta;
              filament_struct(zaehler_f).cor_filament_inv_delta = cor_filament_inv_delta;
              filament_struct(zaehler_f).phi_filament_ext = phi_filament_ext;
              filament_struct(zaehler_f).psi_filament_ext = psi_filament_ext;
              filament_struct(zaehler_f).theta_filament_ext = theta_filament_ext;
              filament_struct(zaehler_f).rot_xyz_filament_ext = rot_filament_ext;
              filament_struct(zaehler_f).ccc_filament_ext = ccc_filament_ext;
              filament_struct(zaehler_f).segment_sensitivity_filament_ext = segment_sensitivity_filament_ext;
              filament_struct(zaehler_f).vect_filament_ext = vect_filament_ext;
              filament_struct(zaehler_f).filament_processed_vec = filament_processed_vec;
              filament_struct(zaehler_f).filament_processed_fraction = filament_processed_fraction;
              filament_struct(zaehler_f).filament_length = size(indx_filament_ext,1);
              
              zaehler_f = zaehler_f + 1;
              
        end
        
    end
    
save(FilamentStruct{k},'filament_struct');

disp(k)

end


clear
    

    
    

%% break V turn filaments

% Define filament length in terms of segments to analyze threshold
filament_length_thres = 3;

% Define filament processed fraction threshold
filament_processed_frac_thres = 0.5;

% angle of two filament
angle_turn = 90;


for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step0/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step1/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir ./mapping3d/filaments_step1/

for k = 1:58
    
    % load filament structure
    load(FilamentStruct{k},'filament_struct');
           
    % new filament structure
    filament_struct_expanded = [];
    filament_struct_N = [];
     
    zaehler_f = 1;
    

    
    
    for i=1:size(filament_struct,2)
     
        %fix the u-turn (over 90 degree turn)
        
        filament_cor_delta = [];
        for j = 2:size(filament_struct(i).cor_filament_ext,1)-1
            V1 = [];
            V2 = [];
            V1 = filament_struct(i).cor_filament_ext(j,:)-filament_struct(i).cor_filament_ext(j-1,:);
            V2 = filament_struct(i).cor_filament_ext(j+1,:)-filament_struct(i).cor_filament_ext(j,:);
            filament_cor_delta(j-1,1) = acosd(dot(V1,V2)/(norm(V1)*norm(V2)));
                
        end
            
        breakpoint = find(filament_cor_delta > angle_turn);
        startpoint = 0;
        ultrafilament = {};
        counter = 1;
        if ~isempty(breakpoint)
            for j = 1:size(breakpoint,1)
                ultrafilament{counter}= (startpoint+1):1:(breakpoint(j)+1);
                startpoint = breakpoint(j)+1;
                counter = counter+1;
            end
        end
        ultrafilament{counter} = (startpoint+1):1:size(filament_struct(i).cor_filament_ext,1);
        
        for j =1: size(ultrafilament,2)
                filament_processed_fraction =...
                    size(find(filament_struct(i).filament_processed_vec(ultrafilament{j},:)==1),1)/...
                    size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
                if (size(ultrafilament{j},2)>= filament_length_thres)  && (filament_processed_fraction >= filament_processed_frac_thres)

                    filament_struct_expanded(zaehler_f).filament_indx = zaehler_f;
                    filament_struct_expanded(zaehler_f).cor_filament_ext = filament_struct(i).cor_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).cor_filament_delta = filament_struct(i).cor_filament_delta(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).cor_filament_inv_delta = filament_struct(i).cor_filament_inv_delta(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).phi_filament_ext = filament_struct(i).phi_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).psi_filament_ext = filament_struct(i).psi_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).theta_filament_ext = filament_struct(i).theta_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).rot_xyz_filament_ext = filament_struct(i).rot_xyz_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).ccc_filament_ext = filament_struct(i).ccc_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).segment_sensitivity_filament_ext = filament_struct(i).segment_sensitivity_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).vect_filament_ext = filament_struct(i).vect_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).filament_processed_vec = filament_struct(i).filament_processed_vec(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).filament_processed_fraction = filament_processed_fraction;
                    filament_struct_expanded(zaehler_f).filament_length = size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);

                    zaehler_f = zaehler_f + 1;
                end
                
        end
    end
    filament_struct = filament_struct_expanded;
    save(FilamentStruct2{k},'filament_struct');
    disp(k)
end

clear

        
%% checking filament independence by the distance and rot angle of adjacent segement

% Define filament length in terms of segments to analyze threshold
filament_length_thres = 3;

% Define filament processed fraction threshold
filament_processed_frac_thres = 0.5;

% Define the longest filament to check
filament_length_check = 3;

% Define the angle range
group_angle = 30;

% Define the distance to the next point in pix
filament_point_D = 60;

for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step1/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step2/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir ./mapping3d/filaments_step2

for k = 1:58
    
    % load filament structure
    load(FilamentStruct{k},'filament_struct');
    
    % new filament structure
    filament_struct_expanded = [];
    zaehler_f = 1;
    
    for i=1:size(filament_struct,2)
        
        % grouping with rot vector
        
        filament_realindx = find(filament_struct(i).filament_processed_vec==1);
        filament_realV = filament_struct(i).rot_xyz_filament_ext(filament_realindx,:);
        
        if  size(filament_realindx,1)> filament_length_check
            psiindex = [1];
            filamentindex = [1:size(filament_realindx,1)];
            subfilament = {};
            fcounter =1;
            while ~isempty(psiindex)
                Vcandidate = filament_realV(filamentindex(1),:);
                counter =1;
                
                psiingroup = [];
                for j = 1:size(filamentindex,2)
                    if abs(acosd(dot(filament_realV(filamentindex(j),:),Vcandidate)/(norm(filament_realV(filamentindex(j),:))*norm(Vcandidate))))<=group_angle
                        psiingroup(counter) = filamentindex(j);
                        counter =counter+1;
                    elseif abs(acosd(dot(filament_realV(filamentindex(j),:),-Vcandidate)/(norm(filament_realV(filamentindex(j),:))*norm(-Vcandidate))))<=group_angle
                        psiingroup(counter) = filamentindex(j);
                        counter =counter+1;
                    end
                end
                subfilament{fcounter} = filament_realindx(psiingroup);
                fcounter = fcounter+1;
                overlap = ismember(filamentindex,psiingroup);
                psiindex = find(overlap == 0);
                filamentindex = filamentindex(psiindex);
                
            end
            % write the non-process filaments into the last group
            subfilament{fcounter} = find(filament_struct(i).filament_processed_vec==0);
            
            % split according to distance
            
            ultrafilament = {};
            counter = 1;
            
            % find the closest next point
            for j = 1:(size(subfilament,2)-1)
                subsubfilament = [];
                filamentpointsD = [];
                counterS = 1;
                counterF = 1;
                pointremain = subfilament{j}';
                subsubfilament(counterF,counterS)=pointremain(1);
                pointremain(1) = [];
                pointoutside = subfilament{size(subfilament,2)};
                while ~isempty(pointremain)
                    filamentpointsD = pdist2(filament_struct(i).cor_filament_ext(subsubfilament(counterF,counterS),1:3),filament_struct(i).cor_filament_ext(pointremain(1,:),1:3));
                    if ~isempty(pointoutside)
                        filamentpointsD2outside = pdist2(filament_struct(i).cor_filament_ext(subsubfilament(counterF,counterS),1:3),filament_struct(i).cor_filament_ext(pointoutside(:,1),1:3));
                    end
                    pointremain(2,:) = filamentpointsD;
                    if min(pointremain(2,:)) <= filament_point_D
                        counterS = counterS + 1;
                        nextpointindx = find(min(pointremain(2,:))==pointremain(2,:));
                        nextpointindx = pointremain(1,nextpointindx);
                        subsubfilament(counterF,counterS) = nextpointindx(1);
                        pointremain(:,find(nextpointindx(1)==pointremain(1,:))) = [];
                    elseif ~isempty(pointoutside) && min(filamentpointsD2outside) <= filament_point_D
                        nextpointindx = find(min(filamentpointsD2outside)==filamentpointsD2outside);
                        nextpointindx = pointoutside(nextpointindx);
                        nextpointindx = nextpointindx(1);
                        if min(pdist2(filament_struct(i).cor_filament_ext(nextpointindx,1:3),filament_struct(i).cor_filament_ext(pointremain(1,:),1:3)))<=20
                            counterS = counterS + 1;
                            subsubfilament(counterF,counterS) = nextpointindx(1);
                            pointoutside(find(pointoutside==nextpointindx)) = [];
                        else
                            counterS = 1;
                            counterF = counterF + 1;
                            subsubfilament(counterF,counterS)=pointremain(1);
                            pointremain(:,1) = [];
                        end
                    else
                        counterS = 1;
                        counterF = counterF + 1;
                        subsubfilament(counterF,counterS)=pointremain(1);
                        pointremain(:,1) = [];
                    end
                end
                
                for t = 1:size(subsubfilament,1)
                    ultrafilament{counter} = subsubfilament(t,find(subsubfilament(t,:)>0));
                    counter = counter + 1;
                end
            end
            pointoutsideR = ismember(pointoutside,cell2mat(ultrafilament));
            ultrafilament{counter} = pointoutside(find(pointoutsideR==0))';
            
            % assign the filaments again
            for j =1: size(ultrafilament,2)
                if ~isempty(ultrafilament{j})
                    filament_processed_fraction =...
                        size(find(filament_struct(i).filament_processed_vec(ultrafilament{j},:)==1),1)/...
                        size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
                    if (size(ultrafilament{j},2)>= filament_length_thres)  && (filament_processed_fraction >= filament_processed_frac_thres)
                        
                        filament_struct_expanded(zaehler_f).filament_indx = zaehler_f;
                        filament_struct_expanded(zaehler_f).cor_filament_ext = filament_struct(i).cor_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).cor_filament_delta = filament_struct(i).cor_filament_delta(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).cor_filament_inv_delta = filament_struct(i).cor_filament_inv_delta(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).phi_filament_ext = filament_struct(i).phi_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).psi_filament_ext = filament_struct(i).psi_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).theta_filament_ext = filament_struct(i).theta_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).rot_xyz_filament_ext = filament_struct(i).rot_xyz_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).ccc_filament_ext = filament_struct(i).ccc_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).segment_sensitivity_filament_ext = filament_struct(i).segment_sensitivity_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).vect_filament_ext = filament_struct(i).vect_filament_ext(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).filament_processed_vec = filament_struct(i).filament_processed_vec(ultrafilament{j},:);
                        filament_struct_expanded(zaehler_f).filament_processed_fraction = filament_processed_fraction;
                        filament_struct_expanded(zaehler_f).filament_length = size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
                        
                        zaehler_f = zaehler_f + 1;
                        
                        
                    end
                end
            end
            
        elseif size(filament_realindx,1) <= filament_length_check
            
            inter_dot_D = [];
            for j=2:size(filament_struct(i).cor_filament_ext,1)
                V1 = filament_struct(i).cor_filament_ext(j,:)-filament_struct(i).cor_filament_ext(j-1,:);
                inter_dot_D(j-1,1) = norm(V1);
            end
            
            breakpoint  = find(inter_dot_D > filament_point_D);
            
            startpoint = 1;
            ultrafilament = {};
            counter = 1;
            if ~isempty(breakpoint)
                for j = 1:size(breakpoint,1)
                    ultrafilament{counter}= startpoint:1:(breakpoint(j));
                    startpoint = breakpoint(j)+1;
                    counter = counter+1;
                end
            end
            ultrafilament{counter} = startpoint:1:size(filament_struct(i).cor_filament_ext,1);
            
            for j =1:size(ultrafilament,2)
                filament_processed_fraction =...
                    size(find(filament_struct(i).filament_processed_vec(ultrafilament{j},:)==1),1)/...
                    size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
                if (size(ultrafilament{j},2)>= filament_length_thres)  && (filament_processed_fraction >= filament_processed_frac_thres)
                    
                    filament_struct_expanded(zaehler_f).filament_indx = zaehler_f;
                    filament_struct_expanded(zaehler_f).cor_filament_ext = filament_struct(i).cor_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).cor_filament_delta = filament_struct(i).cor_filament_delta(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).cor_filament_inv_delta = filament_struct(i).cor_filament_inv_delta(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).phi_filament_ext = filament_struct(i).phi_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).psi_filament_ext = filament_struct(i).psi_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).theta_filament_ext = filament_struct(i).theta_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).rot_xyz_filament_ext = filament_struct(i).rot_xyz_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).ccc_filament_ext = filament_struct(i).ccc_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).segment_sensitivity_filament_ext = filament_struct(i).segment_sensitivity_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).vect_filament_ext = filament_struct(i).vect_filament_ext(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).filament_processed_vec = filament_struct(i).filament_processed_vec(ultrafilament{j},:);
                    filament_struct_expanded(zaehler_f).filament_processed_fraction = filament_processed_fraction;
                    filament_struct_expanded(zaehler_f).filament_length = size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
                    
                    zaehler_f = zaehler_f + 1;
                end
                
            end
            
            
        end
        
        
    end
    
    % reorder segements of each filament again
    
    for i=1:size(filament_struct_expanded,2)
        
        D = pdist2(mean(filament_struct_expanded(i).cor_filament_ext),filament_struct_expanded(i).cor_filament_ext);
        [~,indx_cm] = max(D);
        D = pdist2(filament_struct_expanded(i).cor_filament_ext(indx_cm,:),filament_struct_expanded(i).cor_filament_ext);
        [~,indx_sort] = sort(D);
        filament_struct_expanded(i).cor_filament_ext  = filament_struct_expanded(i).cor_filament_ext(indx_sort,:);
        filament_struct_expanded(i).cor_filament_delta  = filament_struct_expanded(i).cor_filament_delta(indx_sort,:);
        filament_struct_expanded(i).cor_filament_inv_delta  = filament_struct_expanded(i).cor_filament_inv_delta(indx_sort,:);
        filament_struct_expanded(i).phi_filament_ext  = filament_struct_expanded(i).phi_filament_ext(indx_sort,:);
        filament_struct_expanded(i).psi_filament_ext  = filament_struct_expanded(i).psi_filament_ext(indx_sort,:);
        filament_struct_expanded(i).theta_filament_ext  = filament_struct_expanded(i).theta_filament_ext(indx_sort,:);
        filament_struct_expanded(i).rot_xyz_filament_ext  = filament_struct_expanded(i).rot_xyz_filament_ext(indx_sort,:);
        filament_struct_expanded(i).ccc_filament_ext  = filament_struct_expanded(i).ccc_filament_ext(indx_sort,:);
        filament_struct_expanded(i).segment_sensitivity_filament_ext = filament_struct_expanded(i).segment_sensitivity_filament_ext(indx_sort,:);
        filament_struct_expanded(i).vect_filament_ext = filament_struct_expanded(i).vect_filament_ext(indx_sort,:);
        filament_struct_expanded(i).filament_processed_vec = filament_struct_expanded(i).filament_processed_vec(indx_sort,:);
        
    end
    
    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k)
    
end

clear




%% find vector for each filament


% Define psi angle range in degree
filament_psi_range = 30;

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step2/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct3{t} = ['./mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir ./mapping3d/filaments_step3/

for k =1:58
    
    % load filament structure
    load(FilamentStruct2{k},'filament_struct_expanded');
    
    for i=1:size(filament_struct_expanded,2)
        filament_realindx = find(filament_struct_expanded(i).filament_processed_vec==1);
        filament_struct_expanded(i).filament_spin_cat(1:size(filament_struct_expanded(i).cor_filament_ext,1))= 0;
        
        % extract end points of the filament to assign vector
        filamentpoints = filament_struct_expanded(i).cor_filament_ext(:,1:3);
        b = [filamentpoints(end,:)-filamentpoints(1,:)]/norm([filamentpoints(end,:)-filamentpoints(1,:)]);
        
        Regresspsi = b;
        RegresspsiR = -b;
        filament_struct_expanded(i).slope = b;
        
        
        for j = 1:size(filament_realindx,1)
            seg_V = [];
            seg_V = filament_struct_expanded(i).rot_xyz_filament_ext(filament_realindx(j),:);
            if abs(acosd(dot(seg_V,Regresspsi)/(norm(seg_V)*norm(Regresspsi))))<=(filament_psi_range)
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = 1;
            elseif abs(acosd(dot(seg_V,RegresspsiR)/(norm(seg_V)*norm(RegresspsiR))))<=(filament_psi_range)
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = -1;
            else
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = 3;
            end
        end
        
        
        
        
        
        
        portA = size(find(filament_struct_expanded(i).filament_spin_cat==1),2);
        portB = size(find(filament_struct_expanded(i).filament_spin_cat==-1),2);
        portmax = max([portA portB]);
        portconfi = portmax/size(filament_realindx,1);
        if portA>portB
            filament_struct_expanded(i).filament_spin_indx = 1;
            filament_struct_expanded(i).spin_filament_conf = portconfi;
            filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
            filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
        elseif portA<portB
            filament_struct_expanded(i).filament_spin_indx = -1;
            filament_struct_expanded(i).spin_filament_conf = portconfi;
            filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
            filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
        else
            filament_struct_expanded(i).filament_spin_indx = [];
            filament_struct_expanded(i).spin_filament_conf = portconfi;
            filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
            filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
        end
        
    end
    
    
    save(FilamentStruct3{k},'filament_struct_expanded');
    
    disp(k)
end

clear

        
%% confidence score calculation

filament_com_conf_thres = 0.6;

for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for k = 1:58
    load(FilamentStruct{k},'filament_struct_expanded');
    
    for i =1:size(filament_struct_expanded,2)
        
        % Calculate segment sensitivity confidence
        filament_struct_expanded(i).frac_of_pos_sense = size(find(filament_struct_expanded(i).segment_sensitivity_filament_ext==1),1)./size(find(filament_struct_expanded(i).filament_spin_cat~=0),2);
        
        % Calculate mean filament ccc
        filament_struct_expanded(i).mean_filament_ccc = mean(filament_struct_expanded(i).ccc_filament_ext);
        
        % Calculate combined confidence score / majority and sensitivity
        filament_struct_expanded(i).filament_conf_com_score = filament_struct_expanded(i).spin_filament_conf .* filament_struct_expanded(i).frac_of_pos_sense;
    end
    
    
    % Select filaments based on combined confidence score
    filament_struct_ref_conf = [filament_struct_expanded.filament_conf_com_score];
    indx_bad = find(filament_struct_ref_conf < filament_com_conf_thres);
    disp(['Bundle: ' num2str(k) ' | Out of ' num2str(size(filament_struct_ref_conf,2)) ' filaments in this bundle ' num2str(size(indx_bad,2)) ' will be excluded | Number of filaments: ' num2str(size(filament_struct_ref_conf,2)-size(indx_bad,2))]);
    clear filament_struct_ref_conf;
    clear indx_bad;
    save(FilamentStruct{k},'filament_struct_expanded');
    
    disp(k)
end

clear




                        
  
            

        
        
