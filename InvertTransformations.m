function [r_sum_out,t_sum_out] = InvertTransformations(rotations,translations)
%%%%%%%%%%%%
%%%%%%%%%%%%
% This program inverts (and combines) rotations and translations.
%
% INPUT
% rotations --- N x 3 Euler angles, phi, psi, theta
% translations --- N x 3 translations, tx, ty, tz
%
% OUTPUT
% r_sum_out --- 1 x 3 combined Euler angles, phi, psi, theta
% t_sum_out --- 1 x 3 combined translations, tx, ty, tz

% Create rotation matrices
for k=1:size(rotations,1)
    R(:,:,k) = [cosd(rotations(k,1))  sind(rotations(k,1))       0               0;...
               -sind(rotations(k,1))  cosd(rotations(k,1))       0               0;...
                           0                   0                 1               0;...
                           0                   0                 0               1]*...
               [           1                   0                 0               0;...
                           0         cosd(rotations(k,3))  sind(rotations(k,3))  0;...
                           0        -sind(rotations(k,3))  cosd(rotations(k,3))  0;...
                           0                   0                 0               1]*...
               [cosd(rotations(k,2)) sind(rotations(k,2))        0               0;...
               -sind(rotations(k,2)) cosd(rotations(k,2))        0               0;...
                           0                   0                 1               0;...
                           0                   0                 0               1];
end

% Create translation matrices
for k=1:size(translations,1)
    S(:,:,k) = eye(4);
    S(1,4,k) = translations(k,1);
    S(2,4,k) = translations(k,2);
    S(3,4,k) = translations(k,3);
end

% Add transformations
T = eye(4);
for k=1:size(rotations,1)
    T = T*S(:,:,k)*R(:,:,k);
end

% Invert transformations
T = inv(T);

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

