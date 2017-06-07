clear all
load 'plan_closed_set.mat'

plan_retrieved = false;
plan_set = [];
current_state = closed_set(end,1);
current_closed_idx = find(closed_set(:,1)==current_state);
while plan_retrieved == false    
    plan_set = [closed_set(current_closed_idx,:); plan_set];
    prev_step_idx = closed_set(current_closed_idx,2);
    if prev_step_idx ==1
        plan_set = [closed_set(1,:); plan_set];
        plan_retrieved = true;
    else current_closed_idx = find(closed_set(:,1)==prev_step_idx);            
    end            
end