function [plist_rot] = actin_RotatePoints(plist,phi,theta,psi,dim_x,dim_y,dim_z)
%%%%%%%%%%
%%%%%%%%%%

% Translate
plist(:,1) = plist(:,1) - (dim_x/2);
plist(:,2) = plist(:,2) - (dim_y/2);
plist(:,3) = plist(:,3) - (dim_z/2);

% Rotation
for k=1:size(plist,1)
   
    cor_rot = Rx(phi)*Ry(theta)*Rz(psi)*[plist(k,1);plist(k,2);plist(k,3)];
    
    plist_rot(k,1) = cor_rot(1);
    plist_rot(k,2) = cor_rot(2);
    plist_rot(k,3) = cor_rot(3);
    
end

% Translate
plist_rot(:,1) = plist_rot(:,1) + (dim_x/2);
plist_rot(:,2) = plist_rot(:,2) + (dim_y/2);
plist_rot(:,3) = plist_rot(:,3) + (dim_z/2);


% Create rotation matrix Rz
function Rz=Rz(angle)

      Rz = [       cos(angle)               -sin(angle)                0              ;...
                   sin(angle)               cos(angle)                 0              ;...
                        0                         0                      1             ];
      
    
% Create rotation matrix Rx
function Rx=Rx(angle)
     
       Rx = [           1                   0                 0               ;...
                        0               cos(angle)       -sin(angle)        ;...
                        0               sin(angle)       cos(angle)        ];  
                   
                    
% Create rotation matrix Rx
function Ry=Ry(angle)
     
       Ry = [         cos(angle)           0               sin(angle)        ;...
                        0                   1                  0               ;...
                      -sin(angle)          0               cos(angle)        ];  
                    
                    
                    
                    
                    
                    