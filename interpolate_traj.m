function [ exp_obj_traj ] = interpolate_traj( input_traj )
    step_size = 10;
    step_cnt = size(input_traj,1);
    
    for i=1:step_cnt
        for j=1:step_size
            if i==1
                exp_obj_traj((i-1)*step_size+j,:) = input_traj(i,:);
            else
                step_diff = input_traj(i,:)-input_traj(i-1,:);
                exp_obj_traj((i-1)*step_size+j,:) = input_traj(i-1,:) + step_diff*j/step_size;
            end
        end
    end
    

end

