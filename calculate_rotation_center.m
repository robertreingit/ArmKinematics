function [z1,z2] = calculate_rotation_center(marker1,marker2)
% calculate_rotation_center
% Calculates the center of rotation of a ball-joint
% following the routines described in:
% Biryukova et al. (2000)
% INPUT:
% marker1 = no_frames x 6 matrix
%           columns 1-3: position
%           columns 4-6: az-el-ro in radiants
% OUTPUT
% z1    = vector from marker1 to axis
% z2    = vector from marker2 to axis
% SIDEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(size(marker1,2) == 6 ,...
    'marker1 should have 6 columns');
assert(size(marker2,2) == 6, ...
    'marker2 should have 6 columns');
assert(all(size(marker1)==size(marker2)), ...
    'marker1 and marker2 should be of equal size');

no_frames = size(marker1,1);

b1 = zeros(3,1);
b2 = zeros(3,1);
M = zeros(3,3);

% Calculate integral Eq12 + Eq13
for frame = 1:no_frames
    a1 = euler2RM(marker1(frame,4:6));
    a2 = euler2RM(marker2(frame,4:6));
    dr = marker1(frame,1:3)' - marker2(frame,1:3)';
    b1 = b1 - a1 * dr;
    b2 = b2 + a2 * dr;
    M = M + a1*a2';
end

% scale according to number of frames
M = M/no_frames;
b1 = b1/no_frames;
b2 = b2/no_frames;

% Calculate D (Eq14) and
% solve Eq13
D = eye(6);
D(1:3,4:6) = -M;
D(4:6,1:3) = -M';
% (eq15)
B = [b1;b2];
Z = D\B;
z1 = Z(1:3);
z2 = Z(4:6);
