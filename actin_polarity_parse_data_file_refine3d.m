function [image_name,particle_indx,angle_rot,angle_tilt,angle_psi,origin_x,origin_y,class_indx,rlnLogLikeliContribution,rlnMaxValueProbDistribution,rlnNrOfSignificantSamples] = actin_polarity_parse_data_file_refine3d(relion_3d_data_file_mod)
%%%%%%%%
%%%%%%%%%
%
% This function reads a modified Relion Refine3D data file.
%
% INPUT
% relion_3d_data_file_mod --- Name of the modified Relion Refine3D data file 
%
% OUTPUT
% particle_list --- List of particle images
% particle_indx --- Vector of particle indices
% angle_rot --- Relion rot angle
% angle_tilt --- Relion tilt angle
% angle_psi --- Relion psi angle
% origin_x --- Vector of in-plane translations in x-direction
% origin_y --- Vector of in-plane translations in y-direction



% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnRandomSubset #6 
% _rlnGroupNumber #7 
% _rlnAngleRot #8 
% _rlnAngleTilt #9 
% _rlnAnglePsi #10 
% _rlnOriginX #11 
% _rlnOriginY #12 
% _rlnClassNumber #13 
% _rlnNormCorrection #14 
% _rlnLogLikeliContribution #15 
% _rlnMaxValueProbDistribution #16 
% _rlnNrOfSignificantSamples #17 

% Open particles list
fid = fopen(relion_3d_data_file_mod);

% Parse particles list
zaehler = 1;
while 1
    
    tline = fgetl(fid);
    
    if(tline==-1)
         break;
    end
    
    if regexp(tline,'.mrc')
    
    values = regexp(tline,' *','split');
    
    image_name{zaehler} = values{1};
    
    angle_rot(zaehler) = str2double(values{8});    % angle_rot
    angle_tilt(zaehler) = str2double(values{9});   % angle_tilt
    angle_psi(zaehler) = str2double(values{10});    % angle_psi
    
    origin_x(zaehler) = str2double(values{11});
    origin_y(zaehler) = str2double(values{12});
    
    class_indx(zaehler) = str2double(values{13});
    
    rlnLogLikeliContribution(zaehler) = str2double(values{15});
    rlnMaxValueProbDistribution(zaehler) = str2double(values{16});
    rlnNrOfSignificantSamples(zaehler) = str2double(values{17});
    
    zaehler = zaehler + 1;
    
    end

end

% Close particles list
fclose(fid);

% Extract particle indx
for k=1:size(image_name,2)
    image_name_line = image_name{k};
    particle_indx(k) = str2num(image_name_line(1:6));
end

