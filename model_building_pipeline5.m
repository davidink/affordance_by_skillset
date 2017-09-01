clear all
close all

dataFolder = 'data/trajectories/push_2_object/';

num_scene = size(dir([dataFolder 'trajectory*']),1);

%save features and results
opt_save_kernel_feat = false;
opt_save_feat_n_result = true;
for i=1:num_scene
    % load up the scene    
    load_trajectories_mono(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build gp models
build_GP_models_mono(dataFolder);

%test trajectories
testdataFolder = 'data/test/';
num_test_scene = size(dir([testdataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_mono(testdataFolder, i);
end

%dataFolder = 'data/reactive/';

%save features and results
opt_save_kernel_feat = false;
opt_save_feat_n_result = true;
for i=1:num_scene
    % load up the scene    
    load_trajectories_for_reactive_model(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build gp models
build_GP_models_kok(dataFolder);

testdataFolder = 'data/test/';
num_test_scene = size(dir([testdataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_with_react(testdataFolder, i);
end


