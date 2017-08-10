function [euler_angles] = RGete (rotm) 
    
    eul = zeros(3,1);

    % convention used by (*) and (**).
    % note: the final orientation is the same as in XYZ order about fixed axes ...
    if (rotm(1,3) > -1)
        if (rotm(1,3) < 1) % case 1: if r13 ~= ±1
            % Solution with positive sign. It limits the range of the values
            % of theta_y to (-pi/2, pi/2):
            eul(3,1) = atan2(-rotm(1,2), rotm(1,1)); % theta_z
            eul(2,1) = asin(rotm(1,3));            % theta_y
            eul(1,1) = atan2(-rotm(2,3), rotm(3,3)); % theta_x
        else % case 2: if r13 = 1
            % theta_x and theta_z are linked --> Gimbal lock:
            % There are infinity number of solutions for theta_x - theta_z = atan2(-r23, r22).
            % To find a solution, set theta_x = 0 by convention.
            eul(3,1) = atan2(rotm(2,3), rotm(2,2));
            eul(2,1) = pi/2;
            eul(1,1) = 0;
        end
    else % case 3: if r13 = -1
        % Gimbal lock: There is not a unique solution for
        %   theta_x + theta_z = atan2(-r23, r22), by convention, set theta_x = 0.
        eul(1,1) = -atan2(rotm(2,3), rotm(2,2));
        eul(2,1) = -pi/2;
        eul(3,1) = 0;
    end

    euler_angles = eul;
    return;

% 
%     R_x = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
%     R_y = [cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
%     R_z = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
%     
%     R_mat = R_x * R_y * R_z;
%     
%     euler_angles(1) = a;
%     euler_angles(2) = b;
%     euler_angles(3) = g;
end

