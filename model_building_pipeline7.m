clear all
close all

dataFolder = 'data/trajectories/place_1_object/';

num_scene = size(dir([dataFolder 'trajectory*']),1);

%save features and results
opt_save_kernel_feat = false;
opt_save_feat_n_result = true;
for i=1:num_scene
    % load up the scene    
    load_trajectories_place(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build gp models
build_GP_models_place(dataFolder);

%test trajectories
testdataFolder = 'data/test/place_1_object/';
num_test_scene = size(dir([testdataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_place(testdataFolder, i);
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

num_test_scene = size(dir([dataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_with_react(dataFolder, i);
end


