function rotm = quat2rotm(quatm)
    qx = quatm(1);
    qy = quatm(2);
    qz = quatm(3);
    qw = quatm(4);
    rotmat(1,1) = 1 - 2*qy^2-2*qz^2;
    rotmat(1,2) = 2*qx*qy-2*qz*qw;
    rotmat(1,3) = 2*qx*qz+2*qy*qw;
    rotmat(2,1) = 2*qx*qy+2*qz*qw;
    rotmat(2,2) = 1 - 2*qx^2-2*qz^2;
    rotmat(2,3) = 2*qy*qz-2*qx*qw;
    rotmat(3,1) = 2*qx*qz-2*qy*qw;
    rotmat(3,2) = 2*qy*qz+2*qx*qw;
    rotmat(3,3) = 1 - 2*qx^2-2*qy^2;
    
%     rotmat = [1-2*qy^2-2*qz^2 2*qx*qy-2*qz*qw 2*qx*qy+2*qy*qw; 
%                2*qx*qy+2*qz*qw 1-2*qx^2-2*qz^2 2*qy*qz-2*qx*qw;
%                2*qx*qz-2*qy*qw 2*qy*qz+2*qx*qw 1-2*qx^2-2*qy^2];
    rotm =rotmat(:,:);
    
    
end