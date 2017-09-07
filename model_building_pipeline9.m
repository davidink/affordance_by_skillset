clear all
close all

dataFolder = 'data/trajectories/push_2_object_3/';

num_scene = size(dir([dataFolder 'trajectory*']),1);

%save features and results
opt_save_kernel_feat = false;
opt_save_feat_n_result = true;
for i=1:num_scene
    % load up the scene    
    load_trajectories_mono_step(dataFolder,i,opt_save_kernel_feat,opt_save_feat_n_result);
end

%build gp models
build_GP_models_mono_step(dataFolder);

%test trajectories
testdataFolder = 'data/test/';
num_test_scene = size(dir([testdataFolder 'test*']),1);
for i=1:num_test_scene
    test_trajectory_prediction_mono_sam(testdataFolder, i);
end


