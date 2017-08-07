granual = 0.005;
x_size = 0.025;
y_size = 0.025;
z_size = 0.1;
obj_cent = [-0.075 0.0 0];
[modelpoints_l normpoints_l] = create_block_pcd(x_size,y_size,z_size,granual,obj_cent);
obj_cent = [0.075 0.0 0];
[modelpoints_r normpoints_r] = create_block_pcd(x_size,y_size,z_size,granual,obj_cent);

modelpoints = [modelpoints_l;modelpoints_r];
normpoints = [normpoints_l;normpoints_r];

gripper_pose = [eGetR([1.57 0 1.57]) [0 0 0]'; 0 0 0 1];

modelPoints=gripper_pose*[modelpoints'; ones(1,size(modelpoints,1))];
normalPoints = [gripper_pose(1:3,1:3) [0 0 0]';0 0 0 1]*[normpoints'; ones(1,size(normpoints,1))];
modelpoints=modelPoints(1:3,:)';
normalpoints=normalPoints(1:3,:)';

figure;
plot3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),'Color',[1 0 0],'Marker','.','Linestyle','none');
hold on;
%plotCoord(crnt_obj_pose(1:3,4)',crnt_obj_pose(1:3,1:3),0.025);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
quiver3(modelpoints(:,1),modelpoints(:,2),modelpoints(:,3),normalpoints(:,1)/100,normalpoints(:,2)/100,normalpoints(:,3)/100,'Color',[0 0 1]);

save('PR2_gripper.mat','modelpoints','normalpoints');

