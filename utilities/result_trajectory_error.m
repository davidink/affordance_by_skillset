load 'plan_execution_2.mat'

plan_pose = plan_set(:,5:7);
rslt_pose = rslt_obj_pose(:,1:3);
goal_pose = plan_pose(end,:);
start_pose = plan_pose(1,:);
plan_pose(1,:) = [];
plan_pose(end,:) = [];

diff = (plan_pose - rslt_pose).^2;
error = sqrt(diff(:,1) + diff(:,2));

for i=1:size(rslt_pose)
    goal_dist(i) = sqrt((rslt_pose(i,1)-goal_pose(1))^2+(rslt_pose(i,2)-goal_pose(2))^2);
end
