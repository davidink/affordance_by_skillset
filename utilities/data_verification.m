%function [feat rslt] = load_test_dataset()
clear all
close all

dataFolder='~/catkin_ws/data/';
dataSaveFolder=[dataFolder 'data_push_result/'];
load([dataFolder 'models/ARmodels/pair_map.mat']);

current_exp_num = size(dir(dataSaveFolder),1)-2;

for i=1:current_exp_num
    load([dataSaveFolder num2str(i) '/push_command.csv']);
    load([dataSaveFolder num2str(i) '/push_result.csv']);
    num_obj = size(push_result,1)/2;
    obj_poses_wID = push_result(1:num_obj,:);
    obj_poses_wID_ap = push_result(num_obj+1:end,:);
    
    if size(obj_poses_wID,1) ~= 2
        stop = 1;
    end
    %align ids
    if obj_poses_wID_ap(1,1) ~= obj_poses_wID(1,1)
        %swap lines
        temp = obj_poses_wID_ap(1,:);
        obj_poses_wID_ap(1,:) = obj_poses_wID_ap(2,:);
        obj_poses_wID_ap(2,:) = temp;
    end
end
    