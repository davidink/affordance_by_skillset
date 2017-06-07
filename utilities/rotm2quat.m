function quat = rotm2quat(rotm)
    qw = sqrt(1+rotm(1,1)+rotm(2,2)+rotm(3,3))/2;
    qx = (rotm(3,2)-rotm(2,3))/(4*qw);
    qy = (rotm(1,3)-rotm(3,1))/(4*qw);
    qz = (rotm(2,1)-rotm(1,2))/(4*qw);    
    quat = [qx qy qz qw];
end