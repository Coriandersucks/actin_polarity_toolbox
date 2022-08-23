% use the segement of base to rotate tomogram to the same orientation


%% bin tomogram for fast process

for k = 1:size(BasePath,2)
    
    tomo = [BasePath{k} '/overview_wbp.em'];
    tomo = tom_emread(tomo);
    bin_tomo = tom_bin(tomo.Value, 3);
    cd (BasePath{k})
    tom_emwrite('128bin.em',bin_tomo);

    
end

%% return the psi tilt and y shift value

actin_construct_base_path_cells;

for k = 1:size(BasePath,2)
    tomo = [BasePath{k} '/128bin.em'];
    tomo = tom_emread(tomo);
    disp(k);
    tom_volxyz(tomo.Value);
    
    Prmopt = {'psi angle turn'};
    Prmopt2 = {'y_shift'};
    Prmopt3 = {'tilt_angle'}; 
    tilt_angle(k,:) = inputdlg(Prmopt3); 





end
save('./data_matrix_workspace/tilt_angle.mat');


%% write rotation psi and tilt angle into filament structure
load('./data_matrix_workspace/psi_value.mat','psi_delta');
load('./data_matrix_workspace/corb_angle.mat','corb_angle');

%Load filaments and BasePath
BasePath = {};

for t=1:58
    FilamentStruct{t} = ['./mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:58
    FilamentStruct2{t} = ['./mapping3d/filaments_step4/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir ./mapping3d/filaments_step4/

for k=1:58
    
    % load filament structure
    load(FilamentStruct{k},'filament_struct_expanded');
    
    for i = 1:size(filament_struct_expanded,2)
        filament_struct_expanded(i).tomogram_psi_corr = str2num(psi_delta{k,1});
        filament_struct_expanded(i).tomogram_base_corr =corb_angle(k,:);
                
    end
    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k) 
    

end
clear