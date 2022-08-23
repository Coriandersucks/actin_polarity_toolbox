
% Generate BasePath
BasePath = {};

% Define microscope parameter
Objectpixelsize = 0.220651;%nm

% Check particle coordinates
StartPath = pwd;
particle_avg_check = zeros(40,40,40);

for k=1:size(BasePath,2)
    disp('-----------------------------------------------------------------------');
    disp(['tomogram number ' num2str(k)]);
    cd(BasePath{k});
    overview_wbp = tom_emreadc3f('overview_wbp.em');
    try
    
        load([BasePath{k} '/cor/actin_polarity_cor.mat']);
        
        % clear repeated particles
        [plist,index] = unique(actin_3d_cor_all,'rows','stable');
        plist_filaments = actin_filament_list(index);
        
        save([BasePath{k} '/cor/actin_polarity_cor_final_160.mat'],'plist','plist_filaments');
        
        disp(['number of particles in this list: ' num2str(size(plist,1))]);
    catch
        disp('contains no particles in this list');
        continue;
    end
    
    del_indx = [];
    zaehler = 1;
    for i=1:size(plist,1)
         try
              particle = overview_wbp(plist(i,1)-20+1:plist(i,1)+20,plist(i,2)-20+1:plist(i,2)+20,plist(i,3)-20+1:plist(i,3)+20);
              particle_avg_check = particle_avg_check + particle;
         catch
              disp(['tomogram ' num2str(k) ' / particle ' num2str(i) ' delete']);
              del_indx(zaehler) = i;
              zaehler = zaehler + 1;
         end
    end
    
    if ~isempty(del_indx)
        
        plist(del_indx,:) = [];
        plist_filaments(del_indx,:) = [];
        
        save([BasePath{k} '/cor/actin_polarity_cor_final_160.mat'],'plist','plist_filaments');
    end
    clear zaehler;
    clear del_indx;
    cd(StartPath);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
end
clear StartPath i k overview_wbp particle plist;

