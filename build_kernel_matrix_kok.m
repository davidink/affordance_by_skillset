function build_kernel_matrix_kok(dataFolder)
    
    bool_end = false;
    %dataFolder = 'data/reactive/';
    i=0;
    while ~bool_end 
        i = i+1;
        if(size(dir([dataFolder 'kernel' num2str(i) '*']),1)==0)
            bool_end = true;
            scene_cnt = i -1;
        end
    end

    kernel_global = [];
    kernel_local = [];
    kernel_conti = [];
    for i=1:scene_cnt
        frame_num = size(dir([dataFolder 'kernel' num2str(i) '*']),1);
        info = dir([dataFolder 'kernel' num2str(i) '*']);
        for j=1:frame_num
            load([dataFolder info(j).name]);
            kernel_global = [kernel_global;global_shape_features]; 
            kernel_local = [kernel_local;local_contact_shape_features]; 
            kernel_conti = [kernel_conti;relative_pose_features action_pose_features];
        end    
    end

     save([dataFolder 'feat_kernel.mat'], 'kernel_global', 'kernel_local','kernel_conti');
    %save('feat_kernel_react.mat', 'kernel_global', 'kernel_local');
end