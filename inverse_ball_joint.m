%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% inverse_ball_joint
% Calculates the rotation angles from m1 (proximal)
% to m2 (distal)
% INPUT:
% m1 = proximal marker e R^6
% m2 = distal marker e R^6
% OUTPUT:
% euler = Euler angles [azimuth,elevation,roll]
function euler = inverse_ball_joint(m1,m2)

euler = RM2euler(euler2RM(m1(4:6)) * euler2RM(m2(4:6))');
