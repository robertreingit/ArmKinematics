%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caculate_arm_axes
% Calculates a 7DOF model based on the
% routines described in Biryukova.
% INPUT:
% calibration_data
% OUTPUT:
% model = structure with all axis specifications
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = calculate_arm_axes(cd,p)

    if nargin == 1
        p = false;
    end
    
    % extract marker
    wrist_abduction     = cd.wrist_abduction;
    wrist_flexion       = cd.wrist_flexion;
    elbow_pronation     = cd.elbow_pronation;
    elbow_flexion       = cd.elbow_flexion;
    shoulder            = cd.shoulder;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate shoulder
    m1 = get_marker(shoulder,4);
    m2 = get_marker(shoulder,3);
    [tmp_z1,tmp_z2] = calculate_rotation_center(m1,m2);
    model.shoulder.z1 = tmp_z1;
    model.shoulder.z2 = tmp_z2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate elbow
    m1 = get_marker(elbow_flexion,3);
    m2 = get_marker(elbow_flexion,2);
    model.elbow.proximal = collect_skew_axis(m1,m2);
    model.elbow.proximal.name = 'flexion';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate elbow pronation
    m1 = get_marker(elbow_pronation,3);
    m2 = get_marker(elbow_pronation,2);
    model.elbow.distal = collect_skew_axis(m1,m2);
    model.elbow.distal.name = 'pronation';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wrist flexion
    m1 = get_marker(wrist_flexion,2);
    m2 = get_marker(wrist_flexion,1);
    model.wrist.proximal = collect_skew_axis(m1,m2);
    model.wrist.proximal.name = 'flexion';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wrist abduction
    m1 = get_marker(wrist_abduction,2);
    m2 = get_marker(wrist_abduction,1);
    model.wrist.distal = collect_skew_axis(m1,m2);
    model.wrist.distal.name = 'abduction';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check axis directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model = direct_elbow_flexion_axis(model,elbow_flexion);
    model = direct_elbow_pronation_axis(model);
    model = direct_wrist_flexion_axis(model,wrist_flexion);
    model = direct_wrist_abduction_axis(model);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate joint parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.elbow = calculate_skew_params(model.elbow);    
    model.wrist = calculate_skew_params(model.wrist);
    
    if p
        plot_joint(1,shoulder,model,'Shoulder');
        plot_joint(2,elbow_flexion,model,'Elbow flexion');
        plot_joint(3,elbow_pronation,model,'Elbow pronation');
        plot_joint(4,wrist_flexion,model,'Wrist flexion');
        plot_joint(5,wrist_abduction,model,'Wrist abduction');
    end
    
    fprintf('\n*****************************\n');
    fprintf('Calculated arm model.\n');
    fprintf('Elbow: distance = %.3f, angle = %.3f\n', model.elbow.h.value, ...
        rad2deg(model.elbow.a.value));
    fprintf('Wrist: distance = %.3f, angle = %.3f\n', model.wrist.a.value, ...
        rad2deg(model.wrist.a.value));
    fprintf('******************************\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect_skew_axis
% Just a wrapper to calculate_skew_axis
% INPUT:
% m1 = marker 1
% m2 = marker 2
% OUTPUT:
% submodel = struct with w1, w2, z1, z2
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function submodel = collect_skew_axis(m1,m2)

    [tmp_z1,tmp_z2,tmp_w1,tmp_w2] = calculate_skew_axis(m1,m2);
    submodel.z1 = tmp_z1;
    submodel.z2 = tmp_z2;
    submodel.w1 = tmp_w1;
    submodel.w2 = tmp_w2;
    submodel.m1 = m1;
    submodel.m2 = m2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate_skew_params
% Calculates the average distance and angle (rad)
% between the skew axes.
% INPUT:
% m1 = proximal marker
% m2 = distal marker
% z1 = point on proximal axis
% w1 = orientation of proximal axis
% z2 = point on distal axis
% w2 = orientation of distal axis
% OUTPUT:
% joint = results struct with average and std
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function joint = calculate_skew_params( joint )

    m1 = joint.proximal.m1;
    m2 = joint.proximal.m2;
    w1 = joint.proximal.w1;
    z1 = joint.proximal.z1;
    w2 = joint.distal.w2;
    z2 = joint.distal.z2;
    
    no_frames = size(m1,1);
    h = zeros(no_frames,1);
    h_ort = h;
    a = zeros(no_frames,1);
    a_n_h = zeros(no_frames,1);
    
    for f = 1:no_frames
        
        % local to global matrices
        R1_l_g = euler2RM(m1(f,4:6))';
        R2_l_g = euler2RM(m2(f,4:6))';
        
        % determine joint points in global CS
        joint_1 = m1(f,1:3)' + R1_l_g * z1;
        joint_2 = m2(f,1:3)' + R2_l_g * z2;
        
        % calculate w2 in m1 LCS
        w2_local = R1_l_g' * R2_l_g * w2;
        
        % axes in global CS
        w1_global = R1_l_g * w1;
        w2_global = R2_l_g * w2;
        % calculate true orthogonal distance between axes
        [~,~,h_ort(f)] = distance_between_lines(joint_1,w1_global,joint_2,w2_global);
        
        % angle between h vector between axes
        % and normal between both axis
        normal = uvec(cross(w1_global,w2_global));
        h_global = joint_1-joint_2;        
        
        h(f) = norm(joint_1 - joint_2);
        a(f) = angle_between_vectors(w1,w2_local);
        a_n_h(f) = angle_between_vectors(normal,h_global);
        
    end
    
    dist_mean = mean(h);
    dist_std = std(h);
    if dist_mean > 5
        fprintf('Using corrected distance measure.\n');
        dist_mean = mean(h_ort);
        dist_std = std(h_ort);
    end
    
    joint.h.value    = dist_mean;
    joint.h.std      = dist_std;
    joint.a.value    = mean(a);
    joint.a.std      = std(a);
    
    % determine direction of angle between axes
    % use h vector between axes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct_elbow_flexion_axis
% Checks the direction of the elbow axis.
% Calculates the plane in which the elbow flexion
% extension movements are mode ==> approx. sagital
% Calculates the plane normal based on the vectors from
% v1: elbow-flex to shoulder  and v2: elbow-pron to wrist-flex.
% plane normal n = cross(v2,v1). Afterwards
% Checks the angle between elbow flex vector and plane normal
% if greater pi/2 than elbow flex points into the wrong direction.
% INPUT:
% model = model structure
% e_flex = calibration data for elbow flexion
% OUTPUT:
% model = updated model structure
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = direct_elbow_flexion_axis(model,e_flex)

    fprintf('Checking direction of elbow flexion axis\n');

    m4 = get_marker(e_flex,4);
    m3 = get_marker(e_flex,3);
    m2 = get_marker(e_flex,2);

    % Calculates the distances between m4 and m2
    d = matVecNorm(m4(:,1:3) - m2(:,1:3));

    % maximum distance is approximately when elbow is 90° flexed
    [~,max_idx] = max(d);
    
    % establish flexion-extension plane to determine lateral
    % check angle between norm and flexion axis
    shoulder_z2     = model.shoulder.z2;
    elbow_z1        = model.elbow.proximal.z1;
    elbow_z2        = model.elbow.distal.z2;
    wrist_z1        = model.wrist.proximal.z1;
    
    humerus = shoulder_z2 - elbow_z1;
    forearm = wrist_z1 - elbow_z2;
    
    R_m3 = euler2RM(m3(max_idx,4:6))';
    v1 = R_m3 * humerus;
    R_m2 = euler2RM(m2(max_idx,4:6))';
    v2 = R_m2 * forearm;
    % plane normal
    n_lat = uvec(cross(v2,v1));

    % check angle between elbow flexion axis and plane normal
    elbow_flex_axis = R_m3 * model.elbow.proximal.w2;
    angle = rad2deg(angle_between_vectors(elbow_flex_axis,n_lat));
    
    fprintf('Angle betwen plane normal and elbow flexion axis is %.3f°\n', ...
        angle);
    if abs(angle) > 90
        fprintf('Correcting.\n');
        model.elbow.proximal.w1 = -1 * model.elbow.proximal.w1;
    end
    fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct_elbow_pronation_axis
% Checks the angle between the directed forarm vector
% along the forearm from elbow pronation axis to
% wrist flexion and the elbow pronation axis.
% If the angle is <90° the axis vector is inverted.
% INPUT:
% model = model structure
% OUTPUT:
% model = updated model structure
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = direct_elbow_pronation_axis(model)

    % calculate direction from elbow-pronation to wrist-flexion
    % check angle between pronation axis and this vector
    elbow_pronation_axis = model.elbow.distal.w2;
    forearm_direction = -model.elbow.distal.z2 + model.wrist.proximal.z1;
    angle = rad2deg(angle_between_vectors(elbow_pronation_axis,forearm_direction));
    fprintf('Checking direction of elbow pronation axis\n');
    fprintf('Angle between pronation axis and forearm %.3f°\n', angle);
    if angle < 90
        fprintf('Correcting.\n');
        model.elbow.distal.w2 = -1 * model.elbow.distal.w2;
    end
    fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct_wrist_flexion_axis
% Checks the angle between the elbow flexion
% axis and wrist flexion axis during wrist flexion
% axis calibration. Should be in the same direction.
% INPUT:
% model = model structure
% wrist_flexion = wrist_flexion dataset
% OUTPUT:
% model = updated model structure
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = direct_wrist_flexion_axis(model, wrist_flexion)
  
    fprintf('Checking direction of wrist flexion axis.\n');
    m3 = get_marker(wrist_flexion,3);
    m2 = get_marker(wrist_flexion,2);
    no_frames = size(m3,1);
    ang = NaN(no_frames,1);

    for f = 1:no_frames
        axis_1 = euler2RM(m3(f,4:6))' * model.elbow.proximal.w1;
        axis_2 = euler2RM(m2(f,4:6))' * model.wrist.proximal.w1;
        ang(f) = angle_between_vectors(axis_1,axis_2);
    end
    
    angle = rad2deg(median(ang));
    angle_std = std(rad2deg(ang));
    fprintf('Median angle between flexion axis: %3.f°+-%.3f°.\n',...
        angle, angle_std);
    if angle > 90
        fprintf('Correcting.\n');
        model.wrist.proximal.w1 = -1 * model.wrist.proximal.w1;
    end
    fprintf('\n');
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct_wrist_abduction_axis
% Uses the orientation of the marker attached
% to the hand which will always be on the 
% dorsal side and accordingly the z-axis fill face
% palmwards.
% Direct the axis in same direction.
% INPUT:
% OUTPUT:
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = direct_wrist_abduction_axis(model)
    
    fprintf('Checking wrist abduction axis direction.\n');
    w2 = model.wrist.distal.w2;
    z = [0 0 1]';
    if angle_between_vectors(w2,z) > pi/2
        fprintf('Correcting wrist abduction angle.\n');
        model.wrist.distal.w2 = -w2;
    end
    fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_joint
% Simple plotting function
% INPUT:
% pos = index of subplot 1-5
% dataset = marker set of marker 1 - 4
% m1 = proximal marker
% m2 = distal marker
% model = model struct with axis specs
% title_ = string for subplot title
% OUTPUT:
% SIDEEFFECTS:
% Subplot generated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_joint(pos,dataset,model,title_)

    subplot(1,5,pos);
    
    m1 = get_marker(dataset,1);
    m2 = get_marker(dataset,2);   
    m3 = get_marker(dataset,3);
    m4 = get_marker(dataset,4);
    
    no_frames = size(dataset,1);
    frame_id = round(no_frames/2);
    
    m4_pos = m4(frame_id,1:3)';
    m3_pos = m3(frame_id,1:3)';
    m2_pos = m2(frame_id,1:3)';
    m1_pos = m1(frame_id,1:3)';
    
    R4 = euler2RM(m4(frame_id,4:6))';
    R3 = euler2RM(m3(frame_id,4:6))';
    R2 = euler2RM(m2(frame_id,4:6))';
    R1 = euler2RM(m1(frame_id,4:6))';
    
    joints = NaN(3,5);
    joints(:,1) = m4_pos + R4 * model.shoulder.z1;
    joints(:,2) = m3_pos + R3 * model.elbow.proximal.z1;
    joints(:,3) = m2_pos + R2 * model.elbow.distal.z2;
    joints(:,4) = m2_pos + R2 * model.wrist.proximal.z1;
    joints(:,5) = m1_pos + R1 * model.wrist.distal.z2;
    
    axes = NaN(3,5);
    axes(:,2) = R3 * model.elbow.proximal.w1;
    axes(:,3) = R2 * model.elbow.distal.w2;
    axes(:,4) = R2 * model.wrist.proximal.w1;
    axes(:,5) = R1 * model.wrist.distal.w2;
    
    % sub-function for plotting
    function plot_axes(joint,axis,col)
        plot3a(joint,'mo','markerface','magenta');
        if ~any(isnan(axis))
            h = plot_vector(joint,joint + 10*axis);
            set(h,'color',col,'linewidth',2);
        end
    end
    
    % Plotting main axis
    colors = {'','r','g','b','m'};
    plot_marker(dataset);
    set(gca,'nextplot','add');
    for j = 1:5
        plot_axes(joints(:,j),axes(:,j),colors{j});
    end
    set(gca,'nextplot','replace');
    
    % plotting other axes
    title(title_);
    xlabel('X');ylabel('Y');zlabel('Z');

end
