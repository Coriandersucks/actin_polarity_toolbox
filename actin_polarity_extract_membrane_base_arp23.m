%% extract membrane 

% Construct BasePath
BasePath = {};

% Extract particles
for i=1:size(BasePath,2)
    

    % Load tomogram
%     rec = tom_emreadc3f([BasePath{k} '/overview_wbp.em']);
    
    % Load membrane segmentation
    seg = tom_mrcread([BasePath{i} '/overview_wbp2_mem_seg.rec']); % manual_segmentation(rec or mrc file)
    seg = single(seg.Value);
    indx = find(seg ==1);
    
    [x y z] = ind2sub(size(seg),indx);
    mPlist(:,1) = x;
    mPlist(:,2) = y;
    mPlist(:,3) = z;
    
    save([BasePath{i} '/cor/particle_list_manual_seg_mem.mat'],'mPlist');%1024 form
    clearvars mPlist x y z indx seg;
end

for i = 1:46
    load([BasePath{i} '/cor/particle_list_manual_seg_mem.mat'], 'mPlist');
    mPlist(:,3) = mPlist(:,3) + 384;
    save([BasePath{i} '/cor/particle_list_manual_seg_mem.mat'],'mPlist');
    clearvars mPlist
end

for i = 47:58
    load([BasePath{i} '/cor/particle_list_manual_seg_mem.mat'], 'mPlist');
    mPlist(:,3) = mPlist(:,3) + 256;
    save([BasePath{i} '/cor/particle_list_manual_seg_mem.mat'],'mPlist');
    clearvars mPlist
end


%%  extract base and the normal angle

BasePath = {};

for i=1:size(BasePath,2)
    
    % Load membrane segmentation
    try
        seg = tom_mrcread([BasePath{i} '/overview_wbp_b_seg2.rec']);% manual_segmentation(rec or mrc file)
        newbase = 1;
        
    catch
        newbase = 0;
    end
    if newbase == 0
        seg = tom_mrcread([BasePath{i} '/overview_wbp_b_seg.rec']);
        
    end
    seg = single(seg.Value);
    indx = find(seg ==1);
    
    [x y z] = ind2sub(size(seg),indx);
    bPlist(:,1) = x;
    bPlist(:,2) = y;
    bPlist(:,3) = z;
    
    save([BasePath{i} '/cor/particle_list_manual_seg_base.mat'],'bPlist');%1024 form
    clearvars bPlist x y z indx seg;
end

for i = 1:size(BasePath,2)
    load([BasePath{i} '/cor/particle_list_manual_seg_base.mat'], 'bPlist');
    bPlist(:,3) = bPlist(:,3) + 384;
    save([BasePath{i} '/cor/particle_list_manual_seg_base.mat'],'bPlist');
    clearvars bPlist
end


% get normal angle of the base
for i = 1:size(BasePath,2)
    load([BasePath{i} '/cor/particle_list_manual_seg_base.mat'], 'bPlist');
    corbPlist = bPlist;

    [n,V,p] = affine_fit(corbPlist);
    Zv = [0 0 1];
    n = n';
    corR = vrrotvec(n,Zv);
    rotm = axang2rotm(corR);
    eul = rotm2eul(rotm,'XYZ');
    corb_angle(i,:) = eul;
    clear bPlist;

end

save('./data_matrix_workspace/corb_angle.mat');


%% extract arp23
BasePath = {};

xmlfile = arp23_xml2struct('./19-ParticleList.xml');
daughterV = [0 -10 0];
rotm = eul2rotm([0 0 -70/180*pi],'XYZ');
motherV = daughterV*rotm;
acosd(dot(daughterV,motherV)/(norm(daughterV)*norm(motherV)));

for k =1:size(BasePath,2)
    
    text = BasePath{k};
    idx1 = regexp(text, '/', 'end');
    idx2 = size(BasePath{k},2);
    tomoname = text(idx1+1 : idx2);
    [plist1,alg,arpplist] = arp23_2_actin_pytom_pytom2tom(xmlfile,tomoname);
    
    if ~isempty(alg)
        % Number of ribosomes:
        num_of_arp23 = size(alg.indx,2);
        alg_inv = [];

        % Invert combined transformations 
        alg_inv.indx = 1:1:num_of_arp23;
        alg_inv.ccc = zeros(1,num_of_arp23);
        alg_inv.phi = zeros(1,num_of_arp23);
        alg_inv.psi = zeros(1,num_of_arp23);
        alg_inv.theta = zeros(1,num_of_arp23);
        alg_inv.tx = zeros(1,num_of_arp23);
        alg_inv.ty = zeros(1,num_of_arp23);
        alg_inv.tz = zeros(1,num_of_arp23);
        for i=1:num_of_arp23

            rotations = [alg.phi(i) alg.psi(i) alg.theta(i)];

            translations = [alg.tx(i)/4 alg.ty(i)/4 alg.tz(i)/4];

            [r_sum_inv,t_sum_inv] = InvertTransformations(rotations,translations);

            alg_inv.phi(i) = r_sum_inv(1);
            alg_inv.psi(i) = r_sum_inv(2);
            alg_inv.theta(i) = r_sum_inv(3);
            alg_inv.tx(i) = t_sum_inv(1);
            alg_inv.ty(i) = t_sum_inv(2);
            alg_inv.tz(i) = t_sum_inv(3);

        end
    end
    if ~isempty(arpplist)
        arpplist(:,3) = arpplist(:,3)+384;
        for i=1:num_of_arp23
            cor = [];
            cor2 = [];
            corD = [];
            corD2 = [];
            corM = [];
            corM2 = [];
            cor = arpplist(i,:);
            cor2 = cor+[alg_inv.tx(i) alg_inv.ty(i) alg_inv.tz(i)];
            corD = RotatePoints(daughterV,alg_inv.phi(i),alg_inv.theta(i),alg_inv.psi(i),0,0,0);
            corD2 = cor2 + corD;
            corM = RotatePoints(motherV,alg_inv.phi(i),alg_inv.theta(i),alg_inv.psi(i),0,0,0);
            corM2 = cor2 + corM;
            alg_inv.pick_position(i,:) = cor2;
            alg_inv.DaughterV(i,:) = corD2;
            alg_inv.MotherV(i,:) = corM2;
        end
        save([BasePath{k} '/cor/particle_list_manual_seg_arp.mat'],'alg_inv');
    else
        disp(k)
    end

end


clear

