granual = 0.01;
x_size = 0.2;
y_size = 0.05;
z_size = 0.05;
obj_cent = [0.05 -0.075 -0.025];
[modelpoints_l normpoints_l] = create_block_pcd(x_size,y_size,z_size,granual,obj_cent);
obj_cent = [0.05 0.075 -0.025];
[modelpoints_r normpoints_r] = create_block_pcd(x_size,y_size,z_size,granual,obj_cent);

modelpoints = [modelpoints_l;modelpoints_r];
normpoints = [normpoints_l;normpoints_r];

figure;
plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
hold on;
%plotCoord(crnt_obj_pose(1:3,4)',crnt_obj_pose(1:3,1:3),0.025);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal

save('PR2_gripper.mat','modelpoints','normpoints');

