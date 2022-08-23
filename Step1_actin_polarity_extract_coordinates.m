
% Generate BasePath
BasePath = {};

% Extract Actin coordinates
for k=1:size(BasePath,2)
   
    % Load tomogram
    overview_wbp = tom_emreadc3f([BasePath{k} '/overview_wbp.em']);
    
    % Load segmentation
    seg = tom_mrcreadf([BasePath{k} '/tomosegresult.mrc']);
    
    % Create segmentation projection
    sl = tom_bin(mean(overview_wbp.*seg,3),1);
    
    % Extract coordinates
    [actin_3d_cor_all,actin_filament_list] = actin_polarity_extract_cor_3d(seg,12);%The distance between the segments is ~10.5 nm
    
    % Save coordinates
    save([BasePath{k} '/cor/actin_polarity_cor.mat'],'actin_3d_cor_all','actin_filament_list');
    
    % Plot figure without coordinates and save
    figure(1);tom_imagesc(sl);title(['Actin polarity model segmentation: ' num2str(k)]);
    saveas(gcf,['/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/doc/vis_segmentations/tomo_nocor_' num2str(k) '.tif']);
    delete(gcf);
    
    % Plot figure with coordinates and save
    figure(1);tom_imagesc(sl);title(['Actin polarity model segmentation: ' num2str(k)]);
    hold on;plot(actin_3d_cor_all(:,1)./2,actin_3d_cor_all(:,2)./2,'r+');
    saveas(gcf,['/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/doc/vis_segmentations/tomo_cor_' num2str(k) '.tif']);
    delete(gcf);
    
    disp(k);
    
end

% Count all segments
num_of_particles = 0;
all_particles = zeros(1,size(BasePath,2));
for k=1:size(BasePath,2)
    load([BasePath{k} '/cor/actin_polarity_cor.mat'],'actin_3d_cor_all');
    num_of_particles = num_of_particles + size(actin_3d_cor_all,1);
    all_particles(1,k) = size(actin_3d_cor_all,1);
    disp(k);
end
figure(1);plot(all_particles,'LineWidth',1);xlim([1 58]);
ylim([0 10000]);hold on;plot(mean(all_particles).*ones(1,58),'k--');
xlabel('Tomogram index');ylabel('Number of segments');box off;

%228848