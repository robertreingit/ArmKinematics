function [z1,z2,w_1,w_2] = calculate_skew_axis(marker1,marker2,usepinv)
% calculate_skew_axis
% Calculates one of axis according to a skew joint
% following the routines described in:
% Biryukova et al. (2000)
% INPUT:
% marker1 = no_frames x 6 matrix
%           columns 1-3: position
%           columns 4-6: az-el-ro in radiants
% usepinv = switch whether to use pinv from Matlab
%           or calculate inverse according to Biryukova routines.
% OUTPUT
% z1    = vector from marker1 to axis
% z2    = vector from marker2 to axis
% w_1   = axis direction in LCS of marker1
% w_2   = axis direction in LCS of marker2
% SIDEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(size(marker1,2) == 6 ,...
    'marker1 should have 6 columns');
assert(size(marker2,2) == 6, ...
    'marker2 should have 6 columns');
assert(all(size(marker1)==size(marker2)), ...
    'marker1 and marker2 should be of equal size');

if nargin == 2
    usepinv = 0;
end

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
[V,EV] = eig(D);
[~,min_ev] = min(diag(EV));
w_1 = V(1:3,min_ev);
w_1 = w_1/norm(w_1);
w_2 = V(4:6,min_ev);
w_2 = w_2/norm(w_2);

% (eq15)
B = [b1;b2];
% solve Eq6
if usepinv
    pseudo_inv = diag(1./diag(EV));
    pseudo_inv(min_ev,min_ev) = 0;
    Z = V*pseudo_inv*V'*B;
else
    Z = pinv(D)*B;
end
z1 = Z(1:3);
z2 = Z(4:6);
