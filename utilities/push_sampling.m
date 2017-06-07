clear all;
close all;

%Start
dataFolder='~/catkin_ws/data/';
dataSaveFolder=[dataFolder 'data_push_result/'];
%addpath(dataFolder);
load 'PR2_boundary.mat'

current_exp_num = size(dir(dataSaveFolder),1)-2 +1;
current_save_folder = [dataSaveFolder num2str(current_exp_num) '/'];
mkdir(current_save_folder);


% get the starting state
got_scene = false;
disp('Acquiring the initial scene..');
while ~got_scene
    [obj_status obj_result] = system('rosrun affordance_prediction acquire_obj_obs_pose');
    idx = strfind(obj_result,'Found');
    num_obj = str2num(obj_result(1,idx(1)+6));
    if num_obj == 2
        got_scene= true;
        disp('Obtained two object cluster on the table');
    end
end

obj_obs_pose = csvread([dataFolder,'start_obj_obs_pose.txt']);
obj_start_pose = obj_obs_pose(1,:);
obs_pose = obj_obs_pose(2,:);

obj_name = 'Book1';
obs_name = 'Battery1';

table_dim = csvread([dataFolder,'start_table_dimension.txt']);

manipulation_done = false;
crnt_obj_pose = obj_start_pose;

features_result = [];
push_data = [];
push_cnt = 0;

while ~manipulation_done
    push_cnt = push_cnt + 1;
    disp(['Testing no.' num2str(push_cnt) ' push..']);
    % sample a new pushing pose
    %max_r = comput_max_dim(obj_name);
    max_r = 0.15;
    push_dist = 0.15;
    cand_cnt = 20;
    smp_r = max_r + 0.05;
    push_cand = [];
    for i=1:cand_cnt
        %push_cand(i,:) = obj_pose;
        smp_theta(i) = rand * 3.141592 * 2;
        push_cand(i,1) = crnt_obj_pose(1) + smp_r * cos(smp_theta(i));
        push_cand(i,2) = crnt_obj_pose(2) + smp_r * sin(smp_theta(i));
        push_cand(i,3) = crnt_obj_pose(3);
        push_cand_end(i,1) = crnt_obj_pose(1) + (smp_r-max_r) * cos(smp_theta(i));
        push_cand_end(i,2) = crnt_obj_pose(2) + (smp_r-max_r) * sin(smp_theta(i));
        push_cand_end(i,3) = crnt_obj_pose(3);

        smp_theta_rot = smp_theta(i)+ 3.141592;
        m = [0 -sin(smp_theta_rot) cos(smp_theta_rot); 0 cos(smp_theta_rot) sin(smp_theta_rot); -1 0 0];
        qw = sqrt(1 + m(1,1)+m(2,2)+m(3,3))/2;
        push_cand(i,4) = (m(3,2) - m(2,3))/(4*qw);
        push_cand(i,5) = (m(1,3) - m(3,1))/(4*qw);
        push_cand(i,6) = (m(2,1) - m(1,2))/(4*qw);
        push_cand(i,7) = qw; 
        push_cand_end(i,4:7) = push_cand(i,4:7);
    end

    push_cand = [crnt_obj_pose;obs_pose;push_cand];
    csvwrite([dataFolder,'sample_push_candidates.txt'],push_cand);

    % request features of sample points
    disp('Extracting features..');
    [feat_status feat_result] = system('rosrun affordance_prediction compute_feature_from_cand');
    
    % load the features
    features_cand = csvread([dataFolder,'computed_features.txt']);
    
    push_success= false;
    i = 1;
    while ~push_success
        disp(['Trying no.' num2str(i) ' push candidate..']);
        if i == size(push_cand,1)-2
            disp('No push can be done at this pose');
            manipulation_done = true;
            break;
        end
        % check whether the push is feasible
        push_start_ok_right = in_polyhedron(fv_upright_right,push_cand(i+2,1:3));
        push_end_ok_right = in_polyhedron(fv_upright_left,push_cand_end(i,1:3));
        push_start_ok_left = in_polyhedron(fv_upright_right,push_cand(i+2,1:3));
        push_end_ok_left = in_polyhedron(fv_upright_left,push_cand_end(i,1:3));
        
        if (push_start_ok_right && push_end_ok_right) || (push_start_ok_left && push_end_ok_left)
            % execute the push
            current_push = [push_cand(i+2,1:7) push_cand_end(i,1:7)];
            csvwrite([dataFolder,'current_push_command.txt'],current_push);
            disp('Executing the push..');
            [push_status, push_result] = system('rosrun affordance_prediction execute_simple_push');
            if size(strfind(push_result,'successful'),1) == 0
                disp('Pushing was not successful..');
                i = i +1;
            else
                disp('Pushing successful!');
                 %update the state
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
                [got_obj rslt_obj_pose] = extract_obj_pose();
                if got_obj == true
                    features_result = [features_result; features_cand(i,:)];
                    push_data = [push_data; push_cand(i+2,:);crnt_obj_pose;rslt_obj_pose];
                    crnt_obj_pose = rslt_obj_pose;
                    csvwrite([current_save_folder,'features_collected.txt'],features_result);
                    csvwrite([current_save_folder,'results_collected.txt'],push_data);
                else 
                    disp('The object is gone!!');
                    manipulation_done = true;
                end
                push_success = true;
            end
        else i = i + 1;
        end
    end    
end