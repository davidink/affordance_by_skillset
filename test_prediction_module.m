close all
clear all
%dataFolder = 'data/trajectories/push_1_object_4/';
dataFolder = 'data/reactive/';
num_test_scene = size(dir([dataFolder 'test*']),1);
for i=6:num_test_scene
    modular_prediction_push(dataFolder, i);
end