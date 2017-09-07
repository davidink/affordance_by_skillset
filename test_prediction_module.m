close all
clear all
%dataFolder = 'data/trajectories/push_1_object_4/';
testdataFolder = 'data/test/push_2_object/';

num_test_scene = size(dir([testdataFolder 'test*']),1);
for i=1:num_test_scene
    disp(['running scene ' num2str(i)]);
    test_trajectory = [testdataFolder 'test_trajectory' num2str(i)];
    fwd_gp_model = load('data/models/GP_models_feature_ard.mat');
    react_gp_model = load('data/models/GP_models_feature_ard_react.mat');
    modular_prediction_push(test_trajectory, fwd_gp_model, react_gp_model);
end