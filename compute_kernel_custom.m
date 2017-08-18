function kernel_val = compute_kernel_custom(X,XP,theta)
    %%
    X_obj_shape = X(1,1:125);
    XP_obj_shape = XP(1,1:125);
    
    % only consider grid label as "obj=2"
    X_obj_shape(X_obj_shape==1) = 0;
    XP_obj_shape(XP_obj_shape==1) = 0;
    
    % compute kernel using hamming distance
    k_val_obj_shape = sum(~xor(X_obj_shape, XP_obj_shape))/size(X_obj_shape,2);
    
    %%
    X_cont_shape = X(1,126:149);
    XP_cont_shape = XP(1,126:149);
    
    % compute kernel using Bhattacharya kernel
    %for ee local shape
    x1 = [X_cont_shape(1:3) X_cont_shape(7:9)]';
    prec1 = diag([X_cont_shape(4:6) X_cont_shape(10:12)]');
    x2 = [XP_cont_shape(1:3) XP_cont_shape(7:9)]';
    prec2 = diag([XP_cont_shape(4:6) XP_cont_shape(10:12)]');
    k_val_cont_shape_obj1 = GaussBhatKern(x1,prec1,x2, prec2);
    
    %for obj local shape
    offset = 12;
    x1 = [X_cont_shape(offset+1:offset+3) X_cont_shape(offset+7:offset+9)]';
    prec1 = diag([X_cont_shape(offset+4:offset+6) X_cont_shape(offset+10:offset+12)]');
    x2 = [XP_cont_shape(offset+1:offset+3) XP_cont_shape(offset+7:offset+9)]';
    prec2 = diag([XP_cont_shape(offset+4:offset+6) XP_cont_shape(offset+10:offset+12)]');
    k_val_cont_shape_obj2 = GaussBhatKern(x1,prec1,x2, prec2);
    
    %%
    X_relative_pose = X(1,150:155);
    XP_relative_pose = XP(1,150:155);
    k_val_relative_pose = (exp(theta(4))^2)*exp(-(pdist2(X_relative_pose,XP_relative_pose).^2)/(2*exp(theta(3))^2));
    
    %%
    X_action_pose = X(1,156:end);
    XP_action_pose = XP(1,156:end);
    k_val_action_pose = (exp(theta(4))^2)*exp(-(pdist2(XP_action_pose,X_action_pose).^2)/(2*exp(theta(3))^2));
    
    %%
    kernel_val = theta(1)*k_val_obj_shape + theta(2)*k_val_cont_shape_obj1 + theta(2)*k_val_cont_shape_obj2 + k_val_relative_pose + k_val_action_pose;
end

