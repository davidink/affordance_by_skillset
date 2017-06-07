close all;
clear all;

dataFolder = '~/catkin_ws/data/';
graspFolder = [dataFolder 'data_grasp/'];
num_folder = size(dir(graspFolder),1)-2;

grasp_data = [];
for i=1:num_folder
    temp_data = csvread([graspFolder num2str(i) '/grasp_data.txt']);
    grasp_data = [grasp_data;temp_data];
end

figure;
hold on;
gscatter(grasp_data(:,3),grasp_data(:,4),grasp_data(:,1));
table_dim = mean(grasp_data(:,10:15));
plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
xlabel('x');
ylabel('y');
axis([table_dim(1)-0.1 table_dim(2)+0.1 table_dim(3)-0.1 table_dim(4)+0.1]);
camroll(90);
axis equal;

feat_glb_pose = grasp_data(1:10:end,3:9);
rslt_glb_pose = grasp_data(1:10:end,1);

SVMModel = fitcsvm(feat_glb_pose,rslt_glb_pose);
pred_glb_pose = predict(SVMModel,feat_glb_pose);

figure;
gscatter(feat_glb_pose(:,1),feat_glb_pose(:,2),pred_glb_pose,'rb');

%generate a grid over the config space

[X,Y] = meshgrid(table_dim(1):0.05:table_dim(2),table_dim(3):0.05:table_dim(4));
points = [X(:) Y(:)];
mean_z = mean(feat_glb_pose(:,3));
test_pose = [points repmat(mean_z,size(points,1),1)];
% figure;
% plot(test_pose(:,1),test_pose(:,2),'.r');

%generate random quaternion

for i=1:size(points,1)
    smp_theta(i) = rand * 3.141592 * 2;
    smp_theta_rot = smp_theta(i)+ 3.141592;
    m = [0 -sin(smp_theta_rot) cos(smp_theta_rot); 0 cos(smp_theta_rot) sin(smp_theta_rot); -1 0 0];
    qw = sqrt(1 + m(1,1)+m(2,2)+m(3,3))/2;
    test_pose(i,4) = (m(3,2) - m(2,3))/(4*qw);
    test_pose(i,5) = (m(1,3) - m(3,1))/(4*qw);
    test_pose(i,6) = (m(2,1) - m(1,2))/(4*qw);
    test_pose(i,7) = qw; 
end

pred_test_pose = predict(SVMModel,test_pose);
figure;
gscatter(test_pose(:,1),test_pose(:,2),pred_test_pose,'rg');
hold on;
plot3([table_dim(1) table_dim(2) table_dim(2) table_dim(1) table_dim(1)],[table_dim(3) table_dim(3) table_dim(4) table_dim(4) table_dim(3)],[table_dim(6) table_dim(6) table_dim(6) table_dim(6) table_dim(6)],'LineWidth',3);
xlabel('x');
ylabel('y');
axis equal;
axis([table_dim(1)-0.15 table_dim(2)+0.15 table_dim(3)-0.15 table_dim(4)+0.15]);
camroll(90);
legend('Not graspable','Graspable','Location','best');

cand_cnt = 0;
for i=1:size(test_pose,1)
    if pred_test_pose(i)==1
        cand_cnt = cand_cnt+1;
        visualize_model('Book1', test_pose(i,:), [0 1-0.015*cand_cnt 0]);
    end
end


    