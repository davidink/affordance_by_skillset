clear all
dataFolder='~/catkin_ws/data/';
dataSaveFolder = [dataFolder 'data_push_result/'];

num_folder = size(dir(dataSaveFolder),1) -2;

features_merged = [];
results_merged = [];

for i=1:num_folder
    current_folder = [dataSaveFolder num2str(i)];
    current_features = csvread([current_folder '/features_collected.txt']);
    current_results = csvread([current_folder '/results_collected.txt']);
    features_merged = [features_merged;current_features];
    results_merged = [results_merged;current_results];
end

features_set1 = csvread([dataFolder,'features_testset.txt']);
results_set1 = csvread([dataFolder,'results_testset.txt']);

features_merged = [features_merged;features_set1];
results_merged = [results_merged;results_set1];

csvwrite([dataFolder,'features_all.txt'],features_merged);
csvwrite([dataFolder,'results_all.txt'],results_merged);
