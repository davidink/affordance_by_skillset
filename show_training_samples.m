clear all
close all

addpath('utilities');

dataFolder = 'data/reactive/';
num_train_scene = size(dir([dataFolder 'trajec*']),1);
for i=1:num_train_scene
    figure;
    fig_gt = gcf;
    step_size = 10;
    disp_trajectory(fig_gt,[dataFolder 'trajectory' num2str(i)],step_size);
end



    