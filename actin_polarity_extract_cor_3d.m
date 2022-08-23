function [actin_3d_cor_all,filament_list] = actin_polarity_extract_cor_3d(actin_seg,distance_to_next_point)
%%%%%%%%%%%
%%%%%%%%%%%
% This program extracts coordinates along a actin filament.
%
% INPUT
% actin_seg             --- segmentation volume
% distance_to_next_point --- distance to the next point in pixel
%
% plot_flag             --- if plot_flag == 0 no plotting
%                           if plot_flag > 0 plotting in figure(plot_flag)
%
% OUTPUT
% plist_out              --- actin coordinates
% filament_list = filaments identity for each segment

% Skeletonize segmentation
actin_seg = bwskel(logical(actin_seg));

% Exclude dust regions
cc = bwconncomp(actin_seg);
stats = regionprops(cc, 'Area');
idx = find([stats.Area] >= 2.*distance_to_next_point+1);
actin_seg = ismember(labelmatrix(cc), idx);

% Initialize filamnet list
filament_list = [];

% Label and count actin filaments
[actin_seg_label,num_of_actin_filaments] = bwlabeln(actin_seg, 26);

% Main loop over all filaments
actin_3d_cor_all = [];
for k=1:num_of_actin_filaments
   
    % Preprocess single actin filament
    [actin_filament_2d,actin_filament_3d] = actin_preprocess(actin_seg_label,k);
    
    % Extract 2d coordinates of actin segments
    [actin_2d_cor] = actin_extract_points(actin_filament_2d,distance_to_next_point);
    
    % Complement 2d coordinates with missing z-coordinate
    [actin_3d_cor] = Find3dCoordinate(actin_2d_cor,actin_filament_3d);
    
    % % 
    filament_list_pre(1:size(actin_3d_cor,1)) = k;
    % %
    
    % Put together all the coordinates into one matrix
    actin_3d_cor_all = [actin_3d_cor_all;actin_3d_cor];
    
    % % 
    filament_list = [filament_list; filament_list_pre'];
    
    clear filament_list_pre;
    
end


function [actin_filament_2d,actin_filament_3d] = actin_preprocess(actin_seg_label,actin_filament_indx)

% Extract actin filament 3d
actin_filament_3d = (actin_seg_label==actin_filament_indx);

% Extract actin filament 2d
actin_filament_2d = sum(actin_filament_3d,3) > 0;

% Skeletonization
actin_filament_2d = bwmorph(actin_filament_2d,'skel',Inf);


function [plist_out,actin_cor_2d_stats] = actin_extract_points(actin_proj_skel,distance_to_next_point,plot_flag)
%%%%%%%%%%%
%%%%%%%%%%%
% This program extracts coordinates along an actin filament.
%
% INPUT
% actin_proj             --- projection (binary image) of the actin segmentation of one actin filament
% distance_to_next_point --- distance to the next point (eg. the axial rise of F-actin is 27.3 A, and the
%                            optimal box size is ~4 times the axial rise, thus distance_to_next_point
%                            is ~8 pixels)
%
% plot_flag              --- if plot_flag == 0 no plotting
%                            if plot_flag > 0 plotting in figure(plot_flag)
%
% OUTPUT
% plist_out              --- actin coordinates

if nargin == 2
    plot_flag = 0;
end

% Evaluate region properties
STATS = regionprops(actin_proj_skel,'PixelList');

% Extract actin coordinates and switch x and y
xy = STATS.PixelList;
y = xy(:,1);
x = xy(:,2);
plist_all = [x,y];
plist_current = [x,y];

% Calculate the center of mass of the actin filament
actin_center_of_mass = round(mean(plist_current,1));

% Find the starting point
point_current = CalculateMaxDistPoint(plist_current,actin_center_of_mass);

% Display all points
if plot_flag > 0
    figure(plot_flag);tom_imagesc(actin_proj_skel);
    figure(plot_flag);hold on;plot(plist_current(:,1),plist_current(:,2),'g+');
end

% Extract actin coordinates with constant spacing
iteration = 1;

if  size(STATS.PixelList,1) < distance_to_next_point == 1
    
    plist_out = actin_center_of_mass;
    
elseif size(STATS.PixelList,1) >= distance_to_next_point ==1

    % Find the starting point
    point_current = CalculateMaxDistPoint(plist_current,actin_center_of_mass);

    % Display all points
    if plot_flag > 0
        figure(plot_flag);tom_imagesc(actin_proj_skel);
        figure(plot_flag);hold on;plot(plist_current(:,1),plist_current(:,2),'g+');
    end
    
    % Extract actin coordinates with constant spacing
max_iteration = 15000;

    while iteration < max_iteration

        % Stop if all points were visited
        if isempty(plist_current)
            break;
        end

        % Fill output list
        plist_out(iteration,:) = point_current;

        % Display output points
        if plot_flag > 0
             figure(plot_flag);hold on;plot(plist_out(iteration,1),plist_out(iteration,2),'r+');
        end

        % Choose neighbourhood around current point
        [neighbourhood_current] = ChooseNeighbourhood(point_current,plist_all,10);

        % Calculate next point estimation
        [point_est_neg,point_est_pos] = CalculateNextPointEst(neighbourhood_current,point_current,distance_to_next_point);

        % Find point closest to estimated point in negative direction
        [point_est_neg] = CalculateMinDistPoint(plist_current,point_est_neg);

        % Find point closest to estimated point in positive direction
        [point_est_pos] = CalculateMinDistPoint(plist_current,point_est_pos);

        % Maximize the distance between all points in the output list
        % Check point from negative direction
        plist_out_test_est_neg = plist_out;
        plist_out_test_est_neg(end+1,:) = point_est_neg;
        D_est_neg = sum(sum(squareform(pdist(plist_out_test_est_neg)),1),2);

        % Maximize the distance between all points in the output list
        % Check point from positive direction
        plist_out_test_est_pos = plist_out;
        plist_out_test_est_pos(end+1,:) = point_est_pos;
        D_est_pos = sum(sum(squareform(pdist(plist_out_test_est_pos)),1),2);

        % Remove already visited points
        for k=1:size(plist_current,1)
             D_plist_current_point_current(k,1) = sqrt((plist_current(k,1)-point_current(1,1)).^2+(plist_current(k,2)-point_current(1,2)).^2);
        end
        plist_current(D_plist_current_point_current <= distance_to_next_point,:) = [];
        clear D_plist_current_point_current;

        % Select new current point
        if D_est_neg > D_est_pos
            point_current = point_est_neg;
        end

        if D_est_pos > D_est_neg
            point_current = point_est_pos;
        end

        if D_est_pos == D_est_neg
            point_current = point_est_pos;
        end

        % Increase one iteration
        iteration = iteration + 1;

    end
end

% Calculate mean distance between the selected points
if size(plist_out,1) > 1
    for k=1:size(plist_out,1)-1
         D_mean_distance(k,1) = sqrt((plist_out(k,1)-plist_out(k+1,1)).^2+(plist_out(k,2)-plist_out(k+1,2)).^2);
    end
end

% Display number of selected points and mean distance between the selected points
disp(['Number of actin points selected: ' num2str(size(plist_out,1))]);
if size(plist_out,1) > 1
    disp(['Mean distance between the selected points: ' num2str(mean(D_mean_distance,1))]);
    disp(['Standard deviation between the selected points: ' num2str(std(D_mean_distance,1))]);
    actin_cor_2d_stats = mean(D_mean_distance,1);
elseif size(plist_out,1) == 1 
    actin_cor_2d_stats = 0;
end




function [neighbours,idx,neighbor_distance] = ChooseNeighbourhood(plist_point,elist,number_of_neighbours)

iteration_zaehler = 0;

search_switch = 0;
neighbor_distance = 0;
while search_switch == 0
    idx = rangesearch(elist,plist_point,neighbor_distance);
    idx = cell2mat(idx);
    if isempty(idx)
         search_switch = 0;
         neighbor_distance = neighbor_distance + 1;
    else
         if neighbor_distance <= 1
              elist(idx,:) = [];
              search_switch = 0;
              neighbor_distance = neighbor_distance + 1;
         else
              if size(idx,2) >= number_of_neighbours
                   search_switch = 1;
              else
                   search_switch = 0;
                   neighbor_distance = neighbor_distance + 1;
              end
         end
    end
    
    iteration_zaehler = iteration_zaehler + 1;
    
%     if iteration_zaehler == 100
%         search_switch = 1;
%     end
    
end
% if iteration_zaehler == 100
%     neighbours = elist(idx,:);
%     neighbours = neighbours;
% else
    neighbours = elist(idx,:);
    neighbours = neighbours(1:number_of_neighbours,:);
%end




function v = NormalizeVector(vector)
n = sqrt(vector(1).^2 + vector(2).^2);
if n == 0
    v = vector;
else
    v = vector./n;
end


function [point_est_neg,point_est_pos] = CalculateNextPointEst(selected_neighbourhood,plist_point,distance_to_next_point)

% Compute the covariance matrix C
C = cov(selected_neighbourhood);

% Compute the eigenvector of C
[v, lambda] = eig(C);

% Find the eigenvector corresponding to the minimum eigenvalue in C
[~, i] = max(diag(lambda));

% Normalize and scale
point_est_neg(1,:) = (-1) .* distance_to_next_point .* NormalizeVector(v(:,i)') + plist_point;

% Normalize and scale
point_est_pos(1,:) = (+1) .* distance_to_next_point .* NormalizeVector(v(:,i)') + plist_point;


function [point_max_dist] = CalculateMaxDistPoint(plist,p_point)

% Search for pixel that has the highest distance from all of the points
for k=1:size(plist,1)
    dist_plist_to_p_point(1,k) = sqrt((plist(k,1) - p_point(1,1)).^2 + (plist(k,2) - p_point(1,2)).^2);
end
[~,indx_max_pixel] = max(dist_plist_to_p_point);
point_max_dist = plist(indx_max_pixel,:);


function [point_min_dist] = CalculateMinDistPoint(plist,p_point)

% Search for pixel that has the smallest distance from all of the points
for k=1:size(plist,1)
    dist_plist_to_p_point(1,k) = sqrt((plist(k,1) - p_point(1,1)).^2 + (plist(k,2) - p_point(1,2)).^2);
end
[~,indx_min_pixel] = min(dist_plist_to_p_point);
point_min_dist = plist(indx_min_pixel,:);


function [actin_3d_cor] = Find3dCoordinate(actin_2d_cor,actin_3d_seg_single_filament)

actin_3d_cor = ones(size(actin_2d_cor,1),3);
actin_3d_cor(1:end,1:2) = actin_2d_cor;
for k=1:size(actin_2d_cor,1)
   actin_3d_cor(k,3) = round(mean(find(squeeze(actin_3d_seg_single_filament(actin_2d_cor(k,1),actin_2d_cor(k,2),:)) > 0),1));
end
