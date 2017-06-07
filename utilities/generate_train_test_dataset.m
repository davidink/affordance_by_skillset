function [feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(features,results,num_train)
    
    train_idx = randperm(size(features,1),num_train);
    idx = zeros(size(features,1),1);
    for i=1:num_train
        idx(train_idx(1,i),1) = 1;
    end
    train_cnt=0;
    test_cnt=0;    
    for i=1:size(features,1)
        if idx(i,1) == 1 % train
            train_cnt = train_cnt+1;
            feat_train(train_cnt,:) = features(i,:);
            rslt_train(train_cnt,:) = results(i,:);
        else if idx(i,1) ==0 % test
                test_cnt = test_cnt +1;
                feat_test(test_cnt,:) = features(i,:);
                rslt_test(test_cnt,:) = results(i,:);
            end
        end
    end
end