clear all
close all

dataFolder='~/catkin_ws/data/';
dataSaveFolder=[dataFolder 'data_push_result/'];

for i=1:90
    load([dataSaveFolder num2str(i) '/push_command.csv']);
    push_command(1,2) = -0.25;
    push_command(2,2) = -0.0;
    csvwrite([dataSaveFolder num2str(i) '/push_command.csv'],push_command);
end