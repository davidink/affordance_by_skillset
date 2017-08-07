function [R_mat] = eGetR (euler_angles) % euler angles from v-rep
    a = euler_angles(1);
    b = euler_angles(2);
    g = euler_angles(3);
    
    R_x = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
    R_y = [cos(b) 0 sin(b);0 1 0;-sin(b) 0 cos(b)];
    R_z = [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
    
    R_mat = R_x * R_y * R_z;
end

