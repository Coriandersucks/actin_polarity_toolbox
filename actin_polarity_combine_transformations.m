function [r_sum_out,t_sum_out] = actin_polarity_combine_transformations(rotations_tom,translations_tom,angle_psi,angle_tilt,angle_rot,tx_r,ty_r,tz_r)
%%%%%%%%%%%%
%%%%%%%%%%%%
% This program combines rotations and translations.
%
% INPUT
% rotations --- N x 3 Euler angles, phi, psi, theta
% translations --- N x 3 translations, tx, ty, tz
%
% OUTPUT
% r_sum_out --- 1 x 3 combined Euler angles, phi, psi, theta
% t_sum_out --- 1 x 3 combined translations, tx, ty, tz

% Create rotation matrices / TOM
for k=1:size(rotations_tom,1)
    R_tom(:,:,k) = [cosd(rotations_tom(k,1))  sind(rotations_tom(k,1))       0               0;...
               -sind(rotations_tom(k,1))  cosd(rotations_tom(k,1))       0               0;...
                           0                   0                 1               0;...
                           0                   0                 0               1]*...
               [           1                   0                 0               0;...
                           0         cosd(rotations_tom(k,3))  sind(rotations_tom(k,3))  0;...
                           0        -sind(rotations_tom(k,3))  cosd(rotations_tom(k,3))  0;...
                           0                   0                 0               1]*...
               [cosd(rotations_tom(k,2)) sind(rotations_tom(k,2))        0               0;...
               -sind(rotations_tom(k,2)) cosd(rotations_tom(k,2))        0               0;...
                           0                   0                 1               0;...
                           0                   0                 0               1];
end

% Create translation matrices / TOM
for k=1:size(translations_tom,1)
    S_tom(:,:,k) = eye(4);
    S_tom(1,4,k) = translations_tom(k,1);
    S_tom(2,4,k) = translations_tom(k,2);
    S_tom(3,4,k) = translations_tom(k,3);
end

% Prepare Relion angles
rot = -angle_psi.*pi./180;
tilt = -angle_tilt.*pi./180;
psi = -angle_rot.*pi./180;

% Create rotations matrices (right handed coordinate rotations)

    R_rln = [cos(rot)    -sin(rot)       0            0; ...
         sin(rot)     cos(rot)       0            0; ...
            0            0           1            0; ...
            0            0           0            1] * ...
        [cos(tilt)       0        sin(tilt)       0; ...
            0            1           0            0; ...
        -sin(tilt)       0       cos(tilt)        0; ...
            0            0           0            1] * ...
        [cos(psi)    -sin(psi)       0            0; ...
         sin(psi)     cos(psi)       0            0; ...
            0            0           1            0; ...
            0            0           0            1];
     
% Create translation matrice
    S_rln(:,:) = eye(4);
    S_rln(1,4) = tx_r;
    S_rln(2,4) = ty_r;
    S_rln(3,4) = tz_r;

% Add transformations / TOM
T = eye(4);
for k=1:size(rotations_tom,1)
    T = T*S_tom(:,:,k)*R_tom(:,:,k);
end

% Add transformations / RLN
T = T*S_rln(:,:,1)*R_rln(:,:,1);

% Extract translations
t_sum = T*T';

t_sum_out(1) = t_sum(1,4);
t_sum_out(2) = t_sum(2,4);
t_sum_out(3) = t_sum(3,4);

% Extract Euler angles
r_sum_out(1) = (180./pi).*atan2(T(1,3),T(2,3));%phi
r_sum_out(2) = (180./pi).*atan2(T(3,1),-T(3,2));%psi
r_sum_out(3) = (180./pi).*acos(T(3,3));%theta

if -(T(3,3)-1)<10e-8
    r_sum_out(1) = (180./pi).*atan2(T(1,2),T(1,1));
    r_sum_out(2) = 0;
    r_sum_out(3) = 0;
end

