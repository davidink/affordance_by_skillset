function [feat_grp rslt_grp] = generate_cv_dataset(features,results,num_crossfold)
    
    num_div = floor(size(features,1)/num_crossfold);
    num_rm = size(features,1) - num_div * num_crossfold;
    
    rnd_idx = randperm(size(features,1));
    feat_rnd = features(rnd_idx,:);
    rslt_rnd = results(rnd_idx,:);
    
    for i=1:num_crossfold
        if num_rm > 0
            num_add = 1;
            num_rm = num_rm -1;
        else
            num_add =0;
        end
        feat_grp{i} = feat_rnd((i-1)*num_div+1:i*num_div+num_add,:);
        rslt_grp{i} = rslt_rnd((i-1)*num_div+1:i*num_div+num_add,:);
    end
   
end