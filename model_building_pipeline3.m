clear all
close all

dataFolder = 'data/trajectories/push_1_object_4/';

num_scene = size(dir([dataFolder 'trajectory*']),1);

%extract kernel values
opt_save_kernel_feat = true;
opt_save_feat_n_result = false;
for i=1:num_scene
    % load up the scene    
    load_trajectories_kok(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build up kernel matrix
build_kernel_matrix_kok(dataFolder);

%save features and results
opt_save_kernel_feat = false;
opt_save_feat_n_result = true;
for i=1:num_scene
    % load up the scene    
    load_trajectories_kok(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build gp models
build_GP_models_kok(dataFolder);

%test trajectories
num_test_scene = size(dir([dataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_kok(dataFolder, i);
end

dataFolder = 'data/reactive/';
opt_save_kernel_feat = true;
opt_save_feat_n_result = false;
for i=1:num_scene
    load_trajectories_for_reactive_model(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build up kernel matrix
[kernel_global kernel_local] = build_kernel_matrix(dataFolder);

%save features and results
opt_save_kernel_feat = false;
opt_save_feat_n_result = true;
for i=1:num_scene
    % load up the scene    
    load_trajectories_for_reactive_model(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build gp models
build_GP_models(dataFolder);

num_test_scene = size(dir([dataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_with_react(dataFolder, i);
end


