function [vector_z_rot] = actin_polarity_prepare_z_vector(phi,psi,theta)

vector_z_rot = NormalizeVectorEMB(RotateVectorEMB([0 0 1],phi,psi,theta));

function v = NormalizeVectorEMB(vector)
n = sqrt(vector(1).^2 + vector(2).^2 + vector(3).^2);
if n == 0
    v = vector;
else
    v = vector./n;
end

function sp = DotProductEMB(vector1,vector2)

sp = vector1(1,1).*vector2(1,1) + ...
     vector1(1,2).*vector2(1,2) + ...
     vector1(1,3).*vector2(1,3);

function rvector =  RotateVectorEMB(vector,phi,psi,theta)

% Rotation around phi around z-axis
phi_m = [  cosd(phi)    -sind(phi)      0      0; ...
        sind(phi)      cosd(phi)      0      0; ...
        0              0              1      0; ...
        0              0              0      1];

% Rotation around theta around x-axis
theta_m = [   1       0                   0                0;...
         0       cosd(theta)        -sind(theta)      0; ...
         0       sind(theta)         cosd(theta)      0; ...
         0       0                   0                1];

% Rotation around psi around z-axis
psi_m = [cosd(psi)     -sind(psi)       0                0; ...
        sind(psi)      cosd(psi)       0                0; ...
        0              0               1                0; ...
        0              0               0                1]; 

T = psi_m * theta_m * phi_m;

vector = [vector 1];
rvector = T * vector';
rvector = rvector(1:3)';

function rvector =  RotateVectorInvEMB(vector,phi,psi,theta)

% Rotation around phi around z-axis
phi_m = [  cosd(phi)     sind(phi)      0      0; ...
        -sind(phi)      cosd(phi)      0      0; ...
        0              0              1      0; ...
        0              0              0      1];

% Rotation around theta around x-axis
theta_m = [   1       0                   0                0;...
         0       cosd(theta)        sind(theta)      0; ...
         0       -sind(theta)         cosd(theta)      0; ...
         0       0                   0                1];

% Rotation around psi around z-axis
psi_m = [cosd(psi)     sind(psi)       0                0; ...
        -sind(psi)      cosd(psi)       0                0; ...
        0              0               1                0; ...
        0              0               0                1]; 

T = phi_m * theta_m * psi_m;

vector = [vector 1];
rvector = T * vector';
rvector = rvector(1:3)';
