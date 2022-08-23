% Actin polarity / averaging / 

%% Create particle list and alignment structure

% Generate BasePath
BasePath = {};

% Define microscope parameter
Objectpixelsize = 0.220651;%nm

% Define particle path / bin0_104
ParticlesPath = '';

% Loop over actin particles
zaehler = 1;
zaehler_particles = 0;

for k=1:size(BasePath,2)
    
    % Check if particles are selected
    try
         load([BasePath{k} '/cor/actin_polarity_cor_final_160.mat'],'plist');
         particles_selected = 1;
         
    catch
         particles_selected = 0;
    end
    
    % Perform task
    if particles_selected
        
         for p=1:size(plist,1)
             
              % Identify particle
              actin_particle_pre = [ParticlesPath '/actin_particle_t_' num2str(k) '_p_' num2str(p) '.em'];
              disp(actin_particle_pre);
              
              % Create new particle identifier
              actin_particle = [ParticlesPath '/actin_particle_' actin_polarity_construct_particle_identifier(k,p,zaehler) '.em'];
              disp(actin_particle);
              
              % Rename particle
              unix(['mv ' actin_particle_pre ' ' actin_particle]);
              
              particles_list_pre{zaehler} = actin_particle_pre;
              particles_list{zaehler} = actin_particle;
              
              %disp(zaehler);
              
              zaehler = zaehler + 1;
             
         end
         
         zaehler_particles = zaehler_particles + size(plist,1);
         
         clear plist;
         
    end
    
    clear particles_selected;
    
    %disp(k);
    
end

% Build particle list
zaehler = 1;
for k=1:size(particles_list,2)
    
    particles_list_combined{zaehler} = particles_list{k};
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = 'singleaxiswedge 30 80';
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = 'allpass';
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = ' ';
    zaehler = zaehler + 1;
    
end

% Write empty file
placeh = num2str(0);
save('./averaging/actin_polarity_particles_list_bin0_160_clear.txt', 'placeh', '-ASCII');

% Write particles into file
fid = fopen('./averaging/actin_polarity_particles_list_bin0_160_clear.txt', 'w');
for k=1:size(particles_list_combined,2)
    fprintf(fid, '%s\n', char(particles_list_combined{k}));        
end
fclose(fid);

% Save workspace
save('./averaging/workspace_actin_polarity_averaging_160_step_1.mat');
clear all;


%% Project particles

% Parse particles list
[particles_list] = ParseParticleList('./averaging/actin_polarity_particles_list_bin0_160_clear.txt');

% Make projection folder
mkdir ./Particles_Proj_160_z50_clear/

% projection thickness in pixel

pt = 50;

% boxsize in pixel

boxsize = 160;


% Define objectpixelsize
Objectpixelsize = 0.220651 .* 10;% in Angstrom

% Project particles
for k=1:size(particles_list,2)
    
    % Load particle
    particle = double(tom_emreadc3f(particles_list{k}));
    
    % Project actin filament
    proj = mean(particle(:,:,(boxsize/2)-(pt/2)+1:(boxsize/2)+(pt/2)),3);%---> Projection thickness = 11 nm (50 pixels)
    
    % Normalize particle
    proj = (proj - mean(proj(:)))./std(proj(:));
    
    % Dissect file name into parts
    [~,particle_name] = fileparts(particles_list{k});
    
    % Write to disk
    tom_mrcwrite(single(proj),'name',['./Particles_Proj_160_z50_clear/' particle_name '.mrc']);
    
                      % Prepare proj list
                      proj_list{k} = ['./Particles_Proj_160_z50_clear/' particle_name '.mrc' '   ' './Particles_Proj_160_z50_clear/'  particle_name(1:21) '.mrc'];
                      
%     disp(k);
    
end

% Write proj list
placeh = num2str(0);
       save('./averaging/actin_polarity_particles_list_bin0_160_proj_clear.star', 'placeh', '-ASCII');
fid = fopen('./averaging/actin_polarity_particles_list_bin0_160_proj_clear.star', 'w');

% Add the following header:
% data_
% loop_
% _rlnImageName
% _rlnMicrographName
str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];

fprintf(fid, '%s\n', str1, str2, str3, str4);
for k=1:size(proj_list,2)
    fprintf(fid, '%s\n', char(proj_list{k}));
end
fclose(fid);



% Normalize with Relion preprocess
unix(['mpirun -np 64 relion_preprocess_mpi --operate_on ./averaging/actin_polarity_particles_list_bin0_160_proj_clear.star --norm --bg_radius 80 --invert_contrast --operate_out ./averaging/actin_polarity_particles_list_bin0_160_proj_clear_norm_inv --set_angpix 2.20651']);

% Save workspace
save('./averaging/workspace_actin_polarity_averaging_160_step_2.mat');
clear all;


%% Prepare 2D classification with priors / Class2D / job001 / it010 / prealignment with template library

% Parse data file / "this is the prealignment " 
% IMPORTANT: Delete the relion header before parse 
[alg_struct] = actin_polarity_parse_data_file_class2d_prealg('./relion_proj/Class2D/job344/run_it007_data_mod.star');

% Create new starfile
zaehler = 1;
for k=1:size(alg_struct,2)
    
    relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   ' num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
    
    zaehler = zaehler + 1;
    
end

% Write starfile
placeh = num2str(0);
       save('./averaging/particles_job344_it007_prior.star', 'placeh', '-ASCII');
fid = fopen('./averaging/particles_job344_it007_prior.star', 'w');



% Add the following header:
% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 

str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];
str5=["_rlnAngleRotPrior #3"];
str6=["_rlnAngleTiltPrior #4"];
str7=["_rlnAnglePsiPrior #5"];
str8=["_rlnOriginXPrior #6"];
str9=["_rlnOriginYPrior #7"];

fprintf(fid, '%s\n', str1, str2, str3, str4, str5, str6, str7, str8, str9);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);




%% Prepare 3D reconstruction with priors / Class2D / job003 / it100 / Select / job004 / alignment with prior psi angle

% Parse data file / "this is the psi finealignment"
% IMPORTANT: Delete the relion header before parse 

[alg_struct] = actin_polarity_parse_data_file_class2d_finealg('./relion_proj/Select/job365/particles_mod.star');
[particles_list] = ParseParticleList('./averaging/actin_polarity_particles_list_bin0_160_clear.txt');

BasePath = {};

counter = 1;
particleN = [];
for k=1:size(particles_list,2)
    [~,particle_name] = fileparts(particles_list{k});
    name = particle_name;
    idx1 = regexp(name, '_t_', 'end');
    idx2 = regexp(name, '_p_', 'start');
    particleN(counter,1)= str2num(name(idx1+1 : idx2-1));
    counter = counter+1;
    disp(k)
end


counter = 1;
for k = 1:size(BasePath,2)

    try
         load([BasePath{k} '/cor/actin_polarity_cor_final_160.mat']);%plist and plist_filaments
         particles_selected = 1;
         
    catch
         particles_selected = 0;
    end
    
    % Perform task
    if particles_selected == 1
        for p=1:size(plist_filaments,1)
            particleN(counter,2) = counter;
            particleN(counter,3) = plist_filaments(p);
            counter = counter+1;
            disp(counter)
        end
        
        
        
    end
end


% write particles from the same tomogram into the same group to avoid
% overfitting during 3D refinment

zaehler = 1;
counter1 = 0;
counter2 = 0;
for k=1:size(alg_struct,2)
    filamentN = particleN(find(particleN(:,2)==alg_struct(k).particle_indx),3);
    tomoN = particleN(find(particleN(:,2)==alg_struct(k).particle_indx),1);
    
    if rem(filamentN,2)==1
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(tomoN) '   ' num2str(0) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(1)];
         counter1 = counter1 + 1;
    
    elseif rem(filamentN,2)==0
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(tomoN) '   ' num2str(0) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(2)];
         counter2 = counter2 + 1;
    end
    
    zaehler = zaehler + 1;
    disp(k)
end



% Write starfile and add header manually
placeh = num2str(0);
       save('./relion_proj/rsc/particles_select_job365_new.star', 'placeh', '-ASCII');
fid = fopen('./relion_proj/rsc/particles_select_job365_new.star', 'w');
str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];
str5=["_rlnAngleRotPrior #3"];
str6=["_rlnAngleTiltPrior #4"];
str7=["_rlnAnglePsiPrior #5"];
str8=["_rlnRandomSubset #6"];


fprintf(fid, '%s\n', str1, str2, str3, str4, str5, str6, str7, str8);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);


