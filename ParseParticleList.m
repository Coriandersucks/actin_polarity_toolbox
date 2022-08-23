function [particles_list_1_mod,particles_list_2_mod,particles_list_3_mod,particles_list_4_mod] = ParseParticleList(particles_list_name,indx,particles_list_mod_name)
%%%%%%%%%
%%%%%%%%%
%
% This program parses a particle list. Particles with index indx remain in the 
% modified particle list which is written to a text file. If indx is empty all 
% particles remain in the list. If particles_list_mod_name is empty the program 
% only parses the particle list.
%
% INPUT
% particles_list_name --- name of the particle list
% indx --- indx of particles which remain in list 
% particles_list_mod_name --- name of the modified particle list
%
% OUTPUT
% particles_list_1_mod --- list with particles (first line of the particle list)
% particles_list_2_mod --- list with particle wedges (second line of the particle list)
% particles_list_3_mod --- list with real space masks (third line of the particle list)
% particles_list_4_mod --- list with fourier space masks (fourth line of the particle list)

if nargin == 1
    indx = [];
    particles_list_mod_name = [];
end

if nargin == 2
    particles_list_mod_name = [];
end

% Open particle list
fid = fopen(particles_list_name);

% Parse particle list
zaehler = 1;
while 1
    
    tline = fgetl(fid);
    
    if(tline==-1)
         break;
    end
    
    particles_list{zaehler} = tline;
    
    zaehler = zaehler + 1;
end
fclose(fid);

% Find particles path / first line
zaehler = 1;

for k=1:4:size(particles_list,2)
    
    particles_list_1{zaehler} = particles_list{k};
    
    zaehler = zaehler + 1;
end
%disp(['Number of particles parsed: ' num2str(size(particles_list_1,2))]);

% Find particle wedges / second line
zaehler = 1;

for k=2:4:size(particles_list,2)
    
    particles_list_2{zaehler} = particles_list{k};
    
    zaehler = zaehler + 1;
end
%disp(['Number of particle wedges parsed: ' num2str(size(particles_list_2,2))]);

% Find real space masks / third line
zaehler = 1;

for k=3:4:size(particles_list,2)
    
    particles_list_3{zaehler} = particles_list{k};
    
    zaehler = zaehler + 1;
end
%disp(['Number of real space masks parsed: ' num2str(size(particles_list_3,2))]);

% Find fourier space masks / fourth line
zaehler = 1;

for k=4:4:size(particles_list,2)
    
    particles_list_4{zaehler} = particles_list{k};
    
    zaehler = zaehler + 1;
end
%disp(['Number of fourier space masks parsed: ' num2str(size(particles_list_4,2))]);

% Check parsed lines
if size(particles_list_1,2) ~= size(particles_list_2,2) || ...
   size(particles_list_1,2) ~= size(particles_list_3,2) || ...
   size(particles_list_1,2) ~= size(particles_list_4,2)
    error('Number of parsed lines is not equal.');
end

% % Display last parsed lines
% disp(' ');
% disp('Last parsed lines:');
% disp('---');
% disp(particles_list_1{end});
% disp(particles_list_2{end});
% disp(particles_list_3{end});
% disp(particles_list_4{end});
% disp('---');

% Keep particles with indx
if ~isempty(indx)
    
    zaehler = 1;
    
    for k=1:size(particles_list_1,2)
        
         if ~isempty(find(indx==k,1))
            
              particles_list_1_mod{zaehler} = particles_list_1{k};
              particles_list_2_mod{zaehler} = particles_list_2{k};
              particles_list_3_mod{zaehler} = particles_list_3{k};
              particles_list_4_mod{zaehler} = particles_list_4{k};
              
              zaehler = zaehler + 1;
            
         end
        
    end
    
else
    particles_list_1_mod = particles_list_1;
    particles_list_2_mod = particles_list_2;
    particles_list_3_mod = particles_list_3;
    particles_list_4_mod = particles_list_4;
end

% Build modified particle list
zaehler = 1;
for k=1:size(particles_list_1_mod,2)
    
    particles_list_mod{zaehler} = particles_list_1_mod{k};
    
    zaehler = zaehler + 1;
    
    particles_list_mod{zaehler} = particles_list_2_mod{k};
    
    zaehler = zaehler + 1;
    
    particles_list_mod{zaehler} = particles_list_3_mod{k};
    
    zaehler = zaehler + 1;
    
    particles_list_mod{zaehler} = particles_list_4_mod{k};
    
    zaehler = zaehler + 1;
    
end

% Save modified particle list
if ~isempty(particles_list_mod_name)
    
    placeh = num2str(0);
    save(particles_list_mod_name, 'placeh', '-ASCII');
    fid = fopen(particles_list_mod_name, 'w');
    for k=1:size(particles_list_mod,2)
        
        fprintf(fid, '%s\n', char(particles_list_mod{k}));
        
    end
    fclose(fid);

end

% Remove keyword 'emfile' in the list with the particle wedges
for k=1:size(particles_list_2_mod,2)
   
    line_without_keyword_emfile = particles_list_2_mod{k};
    line_without_keyword_emfile = line_without_keyword_emfile(8:end);
    
    particles_list_2_mod_mod{k} = line_without_keyword_emfile;
    
end
particles_list_2_mod = particles_list_2_mod_mod;

