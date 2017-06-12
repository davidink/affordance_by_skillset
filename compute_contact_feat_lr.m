function [feats] = compute_contact_feat_lr(obj_point, obj_point_norm, target_obj_cloud, target_obj_cloud_norm)
    den_est = 0;
    surf_den_est = 0;
    sigma = 0.1;
    for i=1:size(target_obj_cloud,1)
        den_est = den_est + exp(-sqrt(sum((obj_point-target_obj_cloud(i,:)).^2))/sigma^2);
        surf_den_est = surf_den_est + (obj_point_norm*target_obj_cloud_norm(i,:)')*exp(-sqrt(sum((obj_point-target_obj_cloud(i,:)).^2))/sigma^2);
    end
    feats = [den_est surf_den_est 1];
end