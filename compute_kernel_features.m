function [ kernel_feat ] = compute_kernel_features(feat_kernel, grid_feat)
    
    load(feat_kernel);
%     load('feat_kernel2.mat');
%     kernel_global = kernel_global(1:15,:);
%     kernel_local = kernel_local(1:15,:);
    
    size_kernel_global = size(kernel_global,2);
    size_kernel_local = size(kernel_local,2);
    global_grid = grid_feat(:,1:size_kernel_global);    
    local_shape_feat = grid_feat(:,size_kernel_global+1:end);
    
    % Compute kernel features for global 
    kernel_ee = kernel_global;
    kernel_ee(kernel_ee==2) = 0;
    %kernel_ee(kernel_ee==1) = 1;    
    
    kernel_obj = kernel_global;
    kernel_obj(kernel_obj==1) = 0;
    kernel_obj(kernel_obj==2) = 1;
    
    input_grid_ee = global_grid;
    input_grid_ee(input_grid_ee==2) = 0;
    %input_grid_ee(input_grid_ee==) = 1;
    
    input_grid_obj = global_grid;
    input_grid_obj(input_grid_obj==1) = 0;
    input_grid_obj(input_grid_obj==2) = 1;
    
    %self kernel value for normalization
%     for i=1:size(kernel_ee,1)
%         for j=1:size(kernel_ee,1)
%             self_kernel_ee(i,j) = sum(~xor(kernel_ee(i,:),kernel_ee(j,:)))/size(kernel_ee,2);
%         end
%     end
    
    
    %input_kernel_ee = (kernel_ee * input_grid_ee')/(sum(input_grid_ee));
    %input_kernel_obj = (kernel_obj * input_grid_obj')/(sum(input_grid_obj));
    for i=1:size(kernel_ee,1)
        input_kernel_ee(i) = sum(~xor(kernel_ee(i,:),input_grid_ee))/size(input_grid_ee,2);
        input_kernel_obj(i) = sum(~xor(kernel_obj(i,:),input_grid_obj))/size(input_grid_obj,2);
    end
    
    
    kernel_feat = [input_kernel_ee input_kernel_obj];
    
    
    
    % compute kernel for local shape
    input_kernel_local = [];
    
    %for ee local shape
    x1 = [local_shape_feat(1:3) local_shape_feat(7:9)]';
    prec1 = diag([local_shape_feat(4:6) local_shape_feat(10:12)]');
    for i=1:size(kernel_local,1)
        x2 = [kernel_local(i,1:3) kernel_local(i,7:9)]';
        prec2 = diag([kernel_local(i,4:6) kernel_local(i,10:12)]');
        input_kernel_local_ee(1,i) = GaussBhatKern(x1,prec1,x2, prec2);
    end
    
    %for obj local shape
    offset = 12;
    x1 = [local_shape_feat(offset+1:offset+3) local_shape_feat(offset+7:offset+9)]';
    prec1 = diag([local_shape_feat(offset+4:offset+6) local_shape_feat(offset+10:offset+12)]');
    for i=1:size(kernel_local,1)
        x2 = [kernel_local(i,offset+1:offset+3) kernel_local(i,offset+7:offset+9)]';
        prec2 = diag([kernel_local(i,offset+4:offset+6) kernel_local(i,offset+10:offset+12)]');
        input_kernel_local_obj(1,i) = GaussBhatKern(x1,prec1,x2, prec2);
    end
    
    input_kernel_local = [input_kernel_local_ee input_kernel_local_obj];
    
    kernel_feat = [kernel_feat input_kernel_local];    
    
    
    
end

