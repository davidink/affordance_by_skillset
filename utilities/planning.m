%clear all
%close all

load 'gp_a.mat'
load 'gp_b.mat'
load 'gp_nogroup.mat'
load 'GMM_4.mat'
load 'SVM_model.mat'

% get current state
system('rosrun affordance_prediction acquire_obj_obs_pose');
% pause(3);

dataFolder='~/catkin_ws/data/';
obj_obs_pose = csvread([dataFolder,'current_obj_obs_pose.txt']);
obj_start_pose = obj_obs_pose(1,:);
obs_pose = obj_obs_pose(2,:);

table_dim = csvread([dataFolder,'current_table_dimension.txt']);

% set the goal state
goal_state = [0.45 -0.48 0.74]; % edge of the table

% planning
reached_goal = false;
plan_step = 0;
state = [];
F_cand = [];
open_set = [];
closed_set = [];
state_idx = 1;

fig3d = figure;
%rectangle('Position',[table_dim(1) table_dim(3) table_dim(2)-table_dim(1) table_dim(4)-table_dim(3)])
plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
axis equal;
xlabel('x');
ylabel('y');
axis([0.15 1.0 -.6 .6]);
camroll(90);
obj_name = 'Book1';
obs_name = 'Battery1';
hold on;
visualize_model(obj_name,obj_start_pose,[0 1 0]);
visualize_model(obs_name,obs_pose,[1 0 0]);

fig_plan = figure;
hold on;
plot(goal_state(1),goal_state(2),'pg');
plot(obj_start_pose(1),obj_start_pose(2),'pg');
plot(obs_pose(1),obs_pose(2),'or');
axis equal;
axis([table_dim(1) table_dim(2) table_dim(3) table_dim(4)]);
camroll(90);

while reached_goal==false;
    
    % select the state of minimum F from current open set    
    if size(open_set,1)==0 % starting position
        idx_selected = 1;        
        closed_set = [1 1 0 0 obj_start_pose 0 0 0 0 0 0 0]; %start as current_state prev_state g_value h_value
        current_obj_pose = closed_set(1,5:11);           
    else
        [min_F idx_selected] = min(open_set(:,3)+open_set(:,4)); % f = g + h
        closed_set = [closed_set; open_set(idx_selected,:)];
        current_obj_pose = open_set(idx_selected,5:11);
        open_set(idx_selected,:) = [];
    end
        
    % execute the push and update status
    plan_step = plan_step +1;
    
    plot(closed_set(plan_step,5),closed_set(plan_step,6),'pg'); % current state
    if plan_step > 1
        plot(open_set(:,5),open_set(:,6),'.b'); % open sets
        prev_idx = closed_set(plan_step,2);
        prev_idx_in_closed_set = find(closed_set(:,1)==prev_idx);
        plot([closed_set(prev_idx_in_closed_set,5);closed_set(plan_step,5)],[closed_set(prev_idx_in_closed_set,6);closed_set(plan_step,6)],'-y');
    end
        
    %check whether the goal has been reached
    dist_goal = sqrt( (goal_state(1)-current_obj_pose(1))^2 + (goal_state(2)-current_obj_pose(2))^2);
    if dist_goal < 0.05
        reached_goal= true;
        plan_retrieved = false;
        plan_set = [];
        current_state = closed_set(end,1);
        current_closed_idx = find(closed_set(:,1)==current_state);
        while plan_retrieved == false    
            plan_set = [closed_set(current_closed_idx,:); plan_set];
            prev_step_idx = closed_set(current_closed_idx,2);
            if prev_step_idx ==1
                plan_set = [closed_set(1,:); plan_set];
                plan_retrieved = true;
            else current_closed_idx = find(closed_set(:,1)==prev_step_idx);            
            end            
        end
        break;
    end    
    
    % generate sample points of next state
    % sample a pushing position
    cand_cnt = 20;
    current_open_set_size = size(open_set,1);
    
    smp_r = 0.15;
    push_cand = [];
    for i=1:cand_cnt
        state_idx = state_idx + 1;
        state_ids(i) = state_idx;
        %push_cand(i,:) = obj_pose;
        smp_theta(i) = rand * 3.141592 * 2;
        push_cand(i,1) = current_obj_pose(1) + smp_r * cos(smp_theta(i));
        push_cand(i,2) = current_obj_pose(2) + smp_r * sin(smp_theta(i));
        push_cand(i,3) = current_obj_pose(3);

        smp_theta_rot = smp_theta(i)+ 3.141592;
        m = [0 -sin(smp_theta_rot) cos(smp_theta_rot); 0 cos(smp_theta_rot) sin(smp_theta_rot); -1 0 0];
        qw = sqrt(1 + m(1,1)+m(2,2)+m(3,3))/2;
        push_cand(i,4) = (m(3,2) - m(2,3))/(4*qw);
        push_cand(i,5) = (m(1,3) - m(3,1))/(4*qw);
        push_cand(i,6) = (m(2,1) - m(1,2))/(4*qw);
        push_cand(i,7) = qw;    
    end

    push_cand = [current_obj_pose;obs_pose;push_cand];
    csvwrite([dataFolder,'sample_push_candidates.txt'],push_cand);

    % request features of sample points
    system('rosrun affordance_prediction compute_feature_from_cand');

    % compute the output from model
    % load the features
    features_cand = csvread([dataFolder,'computed_features.txt']);

    valid_cand_cnt = 0;
    for i=1:size(features_cand,1)
        for j=1:3
            [y_cand_m_a(j) y_cand_v_a(j)] = predict(gprMdl_a{j},features_cand(i,:));
            [y_cand_m_b(j) y_cand_v_b(j)] = predict(gprMdl_b{j},features_cand(i,:)); 
            [y_cand_m_ng(j) y_cand_v_ng(j)] = predict(gprMdl_init{j},features_cand(i,:)); 
        end
        pred_label = predict(SVMModel,features_cand(i,:));
        if pred_label == 1
            y_cand_m = y_cand_m_a;
        else
            y_cand_m = y_cand_m_b;
        end
%        dist_a = sqrt((y_cand_m_a(1)-gmmodel.mu(1,1))^2+(y_cand_m_a(2)-gmmodel.mu(1,2))^2);
%        dist_b = sqrt((y_cand_m_b(1)-gmmodel.mu(1,1))^2+(y_cand_m_b(2)-gmmodel.mu(1,2))^2);
%         if dist_a < dist_b
%             y_cand_m = y_cand_m_a;
%         else
%             y_cand_m = y_cand_m_b;
%         end
%        y_cand_m = y_cand_m_b;
%        y_cand_m = y_cand_m_ng;
        
        %convert local frame to global frame
        push_transform = [quat2rotm(push_cand(i+2,4:7)) [push_cand(i+2,1) push_cand(i+2,2) push_cand(i+2,3)]';0 0 0 1];
        obj_transform_before_tf = [quat2rotm(current_obj_pose(4:7)) [current_obj_pose(1) current_obj_pose(2) current_obj_pose(3)]';0 0 0 1];
        obj_pose_push_frame = push_transform\ obj_transform_before_tf;
        obj_pose_push_frame_theta = atan2(obj_pose_push_frame(2,3),obj_pose_push_frame(3,3));
        output_pose_theta = obj_pose_push_frame_theta + y_cand_m(3);
        if output_pose_theta < -pi/2
            output_pose_theta = output_pose_theta + 3.141592;
        else if output_pose_theta > pi/2
                output_pose_theta = output_pose_theta - 3.141592;
            end
        end
        output_pose_push_frame(1) = obj_pose_push_frame(1,4); %x
        output_pose_push_frame(2) = obj_pose_push_frame(2,4) + y_cand_m(1); %y
        output_pose_push_frame(3) = obj_pose_push_frame(3,4) + y_cand_m(2); %z
        
        output_pose_push_frame_tf = [1 0 0 output_pose_push_frame(1);
                                      0 cos(output_pose_theta) -sin(output_pose_theta) output_pose_push_frame(2);
                                      0 sin(output_pose_theta) cos(output_pose_theta) output_pose_push_frame(3);
                                      0 0 0 1];
        
        output_pose_tf = push_transform * output_pose_push_frame_tf;
%         output_pose_fix_rot = [0 0 1 0; 0 1 0 0; -1 0 0 0; 0 0 0 1];
%         output_pose_tf = output_pose_tf *output_pose_fix_rot';
        
        output_pose(i,1) = output_pose_tf(1,4);
        output_pose(i,2) = output_pose_tf(2,4);
        output_pose(i,3) = output_pose_tf(3,4);
        output_pose(i,4:7) = rotm2quat(output_pose_tf(1:3,1:3));
        
        %collision_flag = collision_check(obj_name,output_pose(i,:),obs_name,obs_pose);
        collision_flag = false;
        
        if collision_flag == true
        else
            valid_cand_cnt = valid_cand_cnt +1;
            heu_val = heuristic_function(output_pose(i,:),goal_state);        

            open_set(current_open_set_size+valid_cand_cnt,1) = state_ids(i); % current state
            open_set(current_open_set_size+valid_cand_cnt,2) = closed_set(plan_step,1); % previous state
            open_set(current_open_set_size+valid_cand_cnt,3) = closed_set(plan_step,3) + 1; % g = g_prev + cost_to_travel
            open_set(current_open_set_size+valid_cand_cnt,4) = heu_val;
            open_set(current_open_set_size+valid_cand_cnt,5:11) = output_pose(i,:);
            open_set(current_open_set_size+valid_cand_cnt,12:18) = push_cand(i+2,:);
            %F_cand(i) = plan_step + heu_val;
        end
    end
    pause(0.05);    
%     for i=1:size(open_set,1)
%         visualize_model(obj_name, open_set(i,5:11), '.y');
%     end
end

%load 'plan_set.mat'
% for i=1:size(plan_set,1)-1
%     plot([plan_set(i:i+1,5)],[plan_set(i:i+1,6)],'-c');
% %    visualize_model(obj_name, plan_set(i,5:11), '.y')
% end

figure(fig3d);
for i=1:size(plan_set,1)
    visualize_model(obj_name, plan_set(i,5:11), [0 1-0.05*i 0])
end



