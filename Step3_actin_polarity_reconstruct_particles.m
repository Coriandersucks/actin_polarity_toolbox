
% Generate BasePath
BasePath = {};

% Define microscope parameter
Objectpixelsize = 0.220651;%nm

% Define particle path / bin0_144
ParticlesPath = '';

cd (ParticlesPath)

% Perform particle reconstruction
for k= size(BasePath,2)
    
    % Check if particles are selected
    try
         load([BasePath{k} '/cor/actin_polarity_cor_final_160.mat']);
         particles_selected = 1;
         clear plist;
    catch
         particles_selected = 0;
    end
    
    % Perform reconstruction of selected particles
    if particles_selected == 1
          actin_polarity_perform_particle_rec(BasePath{k},ParticlesPath,['actin_particle_t_' num2str(k) '_p_']);
    end
    
    clear particles_selected;
    
    disp(k);
    
end

