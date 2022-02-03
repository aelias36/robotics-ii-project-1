classdef room_cover_controller < handle
    properties
        theta_desired
        v_desired
        v_max

        THETA_DOT_MAX
        X_DOT_MAX

        THETA_DOT_DOT
        X_DOT_DOT

        critical_distance
        epsilon_theta

        T
    end

    methods
        function obj = room_cover_controller()
            obj.theta_desired = 0;
            obj.v_desired = 1;
    
            obj.THETA_DOT_MAX = 2;
            obj.X_DOT_MAX = 1;
            
            obj.THETA_DOT_DOT = 2;
            obj.X_DOT_DOT = 2;

            obj.critical_distance = 1;
            obj.epsilon_theta = deg2rad(10);

            obj.T = 1/15;

        end


        function [v_dot_goal, theta_dot_goal] = step(obj, robot_state, ranges, angles)
            % pick a new direction if close to a wall
            ranges_in_eps = ranges(abs(wrapToPi(angles)) < obj.epsilon_theta);
            if  any(ranges_in_eps < obj.critical_distance)
                ind = find(abs(wrapToPi(angles)) > obj.epsilon_theta...
                    & ranges > max(obj.critical_distance, median(ranges)));
                obj.theta_desired = robot_state(3) + angles(randsample(ind, 1));
                
            end

            % pick velocity goal based on convergence of theta
            if abs(wrapToPi(robot_state(3) - obj.theta_desired)) < obj.epsilon_theta
                obj.v_desired = obj.X_DOT_MAX;
            else
                obj.v_desired = 0;
            end
            
            
            v_dot_goal = obj.v_desired;
            theta_dot_goal = 3*wrapToPi(obj.theta_desired - robot_state(3));
            theta_dot_goal = min(theta_dot_goal, obj.THETA_DOT_MAX);
            theta_dot_goal = max(theta_dot_goal, -obj.THETA_DOT_MAX);

        end
    end
end