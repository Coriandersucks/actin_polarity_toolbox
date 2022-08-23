function [alg_struct] = actin_polarity_parse_data_file_class2d_prealg(relion_data_file)
%%%%%%%%%
%%%%%%%%%
%

% data_images
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnGroupNumber #3 
% _rlnAngleRot #4 
% _rlnAngleTilt #5 
% _rlnAnglePsi #6 
% _rlnOriginX #7 
% _rlnOriginY #8 
% _rlnClassNumber #9 
% _rlnNormCorrection #10 
% _rlnLogLikeliContribution #11 
% _rlnMaxValueProbDistribution #12 
% _rlnNrOfSignificantSamples #13 

% Open particles list
fid = fopen(relion_data_file);

% Parse particles list
zaehler = 1;
while 1
    
    tline = fgetl(fid);
    
    values = regexp(tline,' *','split');
    
    if tline==-1 | size(values,2) < 14
         disp('End of file reached!');
         break;
    end
    
    alg_struct(zaehler).data_line = values;
    
    alg_struct(zaehler).image_name = values{1};
    image_name_pre = alg_struct(zaehler).image_name;
    alg_struct(zaehler).particle_indx = str2num(image_name_pre(1:6));
    alg_struct(zaehler).micrograph_name = values{2};
    micrograph_name_pre = alg_struct(zaehler).micrograph_name;
    alg_struct(zaehler).micrograph_indx = str2num(micrograph_name_pre(123:126));
    
    alg_struct(zaehler).angle_rot = str2double(values{4});
    alg_struct(zaehler).angle_tilt = str2double(values{5});
    alg_struct(zaehler).angle_psi = str2double(values{6});
    alg_struct(zaehler).origin_x = str2double(values{7});
    alg_struct(zaehler).origin_y = str2double(values{8});
    
    alg_struct(zaehler).class_indx = str2double(values{9});
    alg_struct(zaehler).norm_corr = str2double(values{10});
    alg_struct(zaehler).log_likeli_con = str2double(values{11});
    alg_struct(zaehler).max_val_prob_dis = str2double(values{12});
    alg_struct(zaehler).nr_sig_samples = str2double(values{13});
    
    zaehler = zaehler + 1;
    
end

% Close particles list
fclose(fid);

