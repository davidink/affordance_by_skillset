function [pcd pcd_norm] = convert_pcd_frame(input_pcd, input_pcd_norm, ref_frame)
    ref_frame;
    conv_pcd = ref_frame*[input_pcd'; ones(1,size(input_pcd,1))];
    conv_pcd_norms = [ref_frame(1:3,1:3) [0 0 0]';0 0 0 1]*[input_pcd_norm'; ones(1,size(input_pcd_norm,1))];
    pcd=conv_pcd(1:3,:)';
    pcd_norm = conv_pcd_norms(1:3,:)';
end

