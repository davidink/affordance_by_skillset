function draw_contextual_features(features,fig_hd)
    num_sample = size(features,1);
    context_features = features(:,57:end);
    
    figure(fig_hd);
    clf(fig_hd);
    axis equal;
    for i=1:num_sample
        context_feat = context_features(i,:);
        feat_array1 = context_feat(1:4:end);
        feat_array2 = context_feat(2:4:end);
        feat_array3 = context_feat(3:4:end);
        feat_array4 = context_feat(4:4:end);
        y_min = -0.3;
        z_min = -0.1;
        for j=1:sqrt(size(feat_array1,2))
            for k=1:sqrt(size(feat_array1,2))
                r= rectangle('Position',[-0.3+0.05*(j-1) -0.1+0.05*(k-1) 0.05 0.05]);
                if feat_array1((j-1)*11+k) == 1 %obj
                   r.FaceColor = [0 1 0];
                else if feat_array2((j-1)*11+k) == 1 %obs
                        r.FaceColor = [1 0 0 ];
                    else if feat_array3((j-1)*11+k) == 1 %table
                            r.FaceColor = [0.7 0.7 0.7];
                        else if feat_array4((j-1)*11+k) == 1 %free space
                                r.FaceColor = [0 0 0];
                            end
                        end
                    end
                end
            end
        end
    end
end