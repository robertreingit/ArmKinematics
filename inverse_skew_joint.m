%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse_skew_joint
% Calculates the inverse angle for a given
% orientation of markers m1 and m2 using a skew
% model as defined by model1 and model2.
% Prokopenko et al. (2001) p.184
% INPUT:
% joint = joint structure form calculate arm axes
% m1 = proximal marker data
% m2 = distal marker data
% alpha = angle between axes
% OUTPUT:
% ang = angle data in radiants e R^(no_frames,2)
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang = inverse_skew_joint(joint,m1,m2)

    assert(all(size(m1)==size(m2)),'datasets must be equal size.\n');
    
    no_frames = size(m1,1);
    ang = NaN(no_frames,2);
    
    A1 = Azw(joint.proximal.z1,joint.proximal.w1);
    A2 = Azw(-joint.distal.z2,joint.distal.w2);
    alpha = joint.a.value;
    
    for f = 1:no_frames
       
        [~,R1] = marker_data(m1(f,:));
        [~,R2] = marker_data(m2(f,:));
        
        B = A1' * R1' * R2 * A2;
        
        sin_phi = sin(sign(alpha) * acos(B(1,1)));
        % proximal angle
        ang(f,1) = sign(B(3,1)/sin_phi) * acos(B(2,1)/sin_phi);
        % distal angle
        ang(f,2) = sign(B(1,3)/sin_phi) * acos(-B(1,2)/sin_phi); 
        
    end

end

