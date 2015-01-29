%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azw
% Calculates the local coordinate system
% attached to a skew axis according to
% Prokopenko et al. (2001) p.183
% INPUT:
% z = local vector from marker to axis
% omega = axis direction vector
% OUTPUT:
% A = orientation of coordinate system e RO(3)
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Azw(z,omega)

    if isrow(omega)
        omega = omega';
    end

    ox = uvec(omega);
    oy = fast_cross(z,ox);
    oz = fast_cross(ox,oy);
    A = [ox,oy,oz];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast_cross
% cross product function
% optimized for 3d vectors and normalizes
% s
% s = u x v
% INPUT:
% u = vector e R^3
% v = vector e R^3
% OUTPUT:
% s = vector e R^3 |s| = 1
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = fast_cross(u,v)

    s = zeros(3,1);
    s(1) = u(2)*v(3) - u(3)*v(2);
    s(2) = u(3)*v(1) - u(1)*v(3);
    s(3) = u(1)*v(2) - u(2)*v(1);
    s = s / sqrt(s'*s);
        
end
