dataFolder='~/catkin_ws/data/';
table_dim = csvread([dataFolder,'start_table_dimension.txt']);

% Grasp!!
%update the state
disp('Now grasping the object..');
disp('Updating the scene..');
updated = false;
while ~updated
    [status,result] = system('rosrun affordance_prediction save_tabletop_pointcloud');
    err = strfind(result,'error');
    if size(err,1) ==0
        updated = true;
        disp('Scene updated');
    end
end
final_obj_pose = extract_obj_pose();
gripper_pose_above = [final_obj_pose(1) table_dim(3)-0.08 final_obj_pose(3)+0.065 0.7071 0 -0.7071 0];
gripper_pose_grasp = [final_obj_pose(1) table_dim(3)-0.03 final_obj_pose(3)+0.05 0.5 0.5 -0.5 -0.5];
in_above = in_polyhedron(fv_upright_right,gripper_pose_above(1:3));
in_grasp = in_polyhedron(fv_sidegrasp_right,gripper_pose_grasp(1:3));

%plot3(gripper_pose_grasp(1),gripper_pose_grasp(2),gripper_pose_grasp(3),'.r');
if in_above&in_grasp == 1
    %run the grasp
    csvwrite([dataFolder,'current_grasp_command.txt'],[gripper_pose_above gripper_pose_grasp]);
    system('rosrun affordance_prediction execute_simple_grasp');
end