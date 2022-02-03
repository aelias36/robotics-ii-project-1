classdef slam_kalman < handle
    properties
        state 
%         1 p_x
%         2 p_y
%         3 theta
%         4 v_x(R)
%         5 omega
%         6 alpha_x(R)
%         7  A_z
%         8  B_x
%         9  B_z
%         10 C_x
%         11 C_y
%         12 C_z
%         13 D_x
%         14 D_y
%         15 D_z

        cov
        p_tag_front
        p_tag_back
        T
        UWB_cov_mat
        IMU_cov_mat
        odom_cov_mat
        process_cov
    end

    methods
        function obj = slam_kalman()
            obj.state = zeros([15 1]);
            obj.state([7 9 12 15]) = 10; % Assume positive Z

            obj.cov = eye(15) * 1e6;
        end

        function [state_c, cov_c] = multi_meas(obj, meas_UWB, meas_IMU, meas_odom)
        % First, predict
        [state_p, cov_p] = EKF_predict(@transition, @diff_f_wrt_x, @diff_f_wrt_w, obj.process_cov, obj.state, obj.cov, obj.T, 3);
        
        if ~isempty(meas_UWB)
            anchor_mat_p = [0            obj.state(8) obj.state(10) obj.state(13)
                            0            0            obj.state(11) obj.state(14)
                            obj.state(7) obj.state(9) obj.state(12) obj.state(15)];
            dist_vec_p_front = uwb_model(state_p(1:3), obj.p_tag_front, anchor_mat_p);
            dist_vec_p_back = uwb_model(state_p(1:3), obj.p_tag_back, anchor_mat_p);
        end

        if ~isempty(meas_IMU)
            meas_IMU_p = state_p([6 5]); % alpha, omega
        end

        if ~isempty(meas_odom)
            meas_odom_p = state_p([4 5]); % v_x(R), omega
        end
        
        % Then, correct
        dh_dx = [];
        dh_dw = [];
        meas_p = [];
        meas = [];
        meas_cov = [];

        if ~isempty(meas_UWB)
            d_uwb_front_d_x = diff_meas_UWB_wrt_x(state_p, obj.p_tag_front, anchor_mat_p);
            d_uwb_back_d_x = diff_meas_UWB_wrt_x(state_p, obj.p_tag_back, anchor_mat_p);
            dh_dx = [dh_dx; d_uwb_front_d_x; d_uwb_back_d_x];
            dh_dw = blkdiag(dh_dw, eye(8));
            meas_p = [meas_p; dist_vec_p_front'; dist_vec_p_back'];
            meas = [meas; meas_UWB];
            meas_cov = blkdiag(meas_cov, obj.UWB_cov_mat);
        end

        if ~isempty(meas_IMU)
            dh_dx_IMU = [0 0 0 0 0 1
                         0 0 0 0 1 0];
            % extend for slam
            dh_dx_IMU = [dh_dx_IMU zeros([2 9])];
            dh_dx = [dh_dx; dh_dx_IMU];
            dh_dw = blkdiag(dh_dw, eye(2));
            meas_p = [meas_p; meas_IMU_p];
            meas = [meas; meas_IMU];
            meas_cov = blkdiag(meas_cov, obj.IMU_cov_mat);
        end

        if ~isempty(meas_odom)
            dh_dx_odom= [0 0 0 1 0 0
                         0 0 0 0 1 0];
            % extend for slam
            dh_dx_odom = [dh_dx_odom zeros([2 9])];
            dh_dx = [dh_dx; dh_dx_odom];
            dh_dw = blkdiag(dh_dw, eye(2));
            meas_p = [meas_p; meas_odom_p];
            meas = [meas; meas_odom];
            meas_cov = blkdiag(meas_cov, obj.odom_cov_mat);
        end

        [state_c, cov_c] = EKF_correct(dh_dx, dh_dw,state_p, cov_p, meas_p, meas, meas_cov, 3, []);

        obj.state = state_c;
        obj.cov = cov_c;
        end
    end
end

function df_dx = diff_f_wrt_x(state, T)
c = cos(state(3));
s = sin(state(3));

df_dx = [1 0  s*state(4)*T c*T 0 0  % p_x
         0 1 -c*state(4)*T s*T 0 0  % p_y
         0 0 1 0 T 0    % theta
         0 0 0 1 0 T    % v_x (R)
         0 0 0 0 1 0    % omega
         0 0 0 0 0 1];  % alpha_x (R)

% extend for slam
df_dx = blkdiag(df_dx, eye(9));

end

function df_dw = diff_f_wrt_w(state)
    df_dw = eye(6);

    % extend for slam
    df_dw = [df_dw; zeros([9 6])];
end

function next_state = transition(state, T)
c = cos(state(3));
s = sin(state(3));
A = [1 0 0 c*T 0 0  % p_x
         0 1 0 s*T 0 0  % p_y
         0 0 1 0 T 0    % theta
         0 0 0 1 0 T    % v_x (R)
         0 0 0 0 1 0    % omega
         0 0 0 0 0 1];  % alpha_x (R)
    % extend for slam
    A = blkdiag(A, eye(9));
    next_state = A*state;
    next_state(3) = wrapToPi(next_state(3));
end

function [state_p, cov_p] = EKF_predict(transitionfcn, df_dx_fcn, df_dw_fcn, process_cov, state, cov, dt, state_angle_ind)
    state_p = transitionfcn(state, dt);
    state_p(state_angle_ind) = wrapToPi(state_p(state_angle_ind));
    
    df_dx = df_dx_fcn(state, dt);
    df_dw = df_dw_fcn(state);
    cov_p = df_dx*cov*df_dx' + df_dw*process_cov*df_dw';
%         
%         meas_p = measurementfcn(state_p);
%         meas_p(meas_angle_ind) = wrapToPi(meas_p(meas_angle_ind));
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

function d_h_d_x = diff_meas_UWB_wrt_x(state, tag_pos, anchor_mat)
    % contribution from robot / tags
    robot_tag_world = rot2(state(3))*tag_pos;
    tag_world = state(1:2) + rot2(state(3))*tag_pos;
    p_anchor_tag = [tag_world; 0] - anchor_mat;
    d_h_d_tagxyz = NaN(size(anchor_mat));
    for i = 1:width(anchor_mat)
        d_h_d_tagxyz(:,i) = p_anchor_tag(:,i) / norm(p_anchor_tag(:,i));
    end
    d_h_d_tagxyz = d_h_d_tagxyz';

    
    d_tagxyz_d_xrob = [1 0 -robot_tag_world(2) 0 0 0
                    0 1  robot_tag_world(1) 0 0 0
                    0 0 0 0 0 0];

    d_h_d_xrob = d_h_d_tagxyz * d_tagxyz_d_xrob;

    % contribution from anchors
    
    % normal vectors pointing from tag to anchors
    n = -d_h_d_tagxyz; % (anchor, xyz)
%         1 A_z
%         2 B_x
%         3 B_z
%         4 C_x
%         5 C_y
%         6 C_z
%         7 D_x
%         8 D_y
%         9 D_z

    d_h_d_xanchor = [n(1,3) 0      0      0      0      0      0      0      0
                     0      n(2,1) n(2,3) 0      0      0      0      0      0
                     0      0      0      n(3,1) n(3,2) n(3,3) 0      0      0
                     0      0      0      0      0      0      n(4,1) n(4,2) n(4,3)];


    d_h_d_x = [d_h_d_xrob d_h_d_xanchor];
end