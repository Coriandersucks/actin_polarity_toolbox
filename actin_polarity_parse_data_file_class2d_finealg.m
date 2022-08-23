function [alg_struct] = actin_polarity_parse_data_file_class2d_finealg(relion_data_file)
%%%%%%%%%
%%%%%%%%%
%

% # RELION; version 3.0.4
% 
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
% _rlnGroupNumber #8 
% _rlnAngleRot #9 
% _rlnAngleTilt #10 
% _rlnAnglePsi #11 
% _rlnOriginX #12 
% _rlnOriginY #13 
% _rlnClassNumber #14 
% _rlnNormCorrection #15 
% _rlnLogLikeliContribution #16 
% _rlnMaxValueProbDistribution #17 
% _rlnNrOfSignificantSamples #18 

% Open particles list
fid = fopen(relion_data_file);

% Parse particles list
zaehler = 1;
while 1
    
    tline = fgetl(fid);
    
    values = regexp(tline,' *','split');
    
    if tline==-1 | size(values,2) < 19
         disp('End of file reached!');
         break;
    end
    
    alg_struct(zaehler).data_line = values;
    
    alg_struct(zaehler).image_name = values{1};
    image_name_pre = alg_struct(zaehler).image_name;
    alg_struct(zaehler).particle_indx = str2num(image_name_pre(1:6));
    alg_struct(zaehler).micrograph_name = values{2};
    micrograph_name_pre = alg_struct(zaehler).micrograph_name;
    alg_struct(zaehler).micrograph_indx = str2num(micrograph_name_pre(96:100));
    
    alg_struct(zaehler).angle_rot = str2double(values{9});
    alg_struct(zaehler).angle_tilt = str2double(values{10});
    alg_struct(zaehler).angle_psi = str2double(values{11});
    alg_struct(zaehler).origin_x = str2double(values{12});
    alg_struct(zaehler).origin_y = str2double(values{13});
    
    alg_struct(zaehler).class_indx = str2double(values{14});
    alg_struct(zaehler).norm_corr = str2double(values{15});
    alg_struct(zaehler).log_likeli_con = str2double(values{16});
    alg_struct(zaehler).max_val_prob_dis = str2double(values{17});
    alg_struct(zaehler).nr_sig_samples = str2double(values{18});
    disp(zaehler)
    zaehler = zaehler + 1;
    
end

% Close particles list
fclose(fid);

