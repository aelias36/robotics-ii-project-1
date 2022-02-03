classdef mapping_kalman < handle

properties
    state % [x; y; z] for a single anchor
    cov
    p_tag_front
    p_tag_back
    UWB_cov
    k_anchor % anchor number (1~4)
end

methods
    function obj = mapping_kalman()
        obj.state = [0;0;10]; % Assume anchor above
        obj.cov = eye(3) * 1e6;
    end

    function [state_c, cov_c] = UWB_meas(obj, meas_UWB, robot_state)
        % There are no system dynamics, so only correction step is needed

        % Convert to 2 measurements: [d_front; d_back]
        N_anchors = 4;
        meas = meas_UWB(obj.k_anchor + [0 N_anchors]);

        dist_p_front = uwb_model(robot_state, obj.p_tag_front, obj.state);
        dist_p_back = uwb_model(robot_state,  obj.p_tag_back, obj.state);
        meas_p = [dist_p_front; dist_p_back];

        dh_dx_1 = diff_meas_UWB_wrt_anchor(robot_state, obj.p_tag_front, obj.state);
        dh_dx_2 = diff_meas_UWB_wrt_anchor(robot_state, obj.p_tag_back, obj.state);


        [state_c, cov_c] = EKF_correct([dh_dx_1;dh_dx_2], eye(2), obj.state, obj.cov, meas_p, meas, obj.UWB_cov, [], []);
        obj.state = state_c;
        obj.cov = cov_c;
    end
end
end

function [state_c, cov_c] = EKF_correct(dh_dx, dh_dw, state_p, cov_p, meas_p, meas, meas_cov, state_angle_ind, meas_angle_ind)
    S = dh_dx*cov_p*dh_dx' + dh_dw*meas_cov*dh_dw';
    K = cov_p*dh_dx'/S;
    
    meas_diff = meas - meas_p;
    meas_diff(meas_angle_ind) = wrapToPi(meas_diff(meas_angle_ind));
    
    state_c = state_p + K*meas_diff;
    state_c(state_angle_ind) = wrapToPi(state_c(state_angle_ind));
    
    cov_c = cov_p - K*S*K';
    cov_c = (cov_c + cov_c')/2;

end

function d_h_d_x = diff_meas_UWB_wrt_anchor(robot_state, tag_pos, anchor)
    tag_world = robot_state(1:2) + rot2(robot_state(3))*tag_pos;
    p_tag_anchor = anchor - [tag_world; 0];
    d_h_d_x = p_tag_anchor' / norm(p_tag_anchor);

end