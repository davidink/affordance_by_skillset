clear all
close all
load 'features_n_results_interaction1';
results_all = [results results_interaction];

for k=1:5

    [feat_train rslt_train feat_test rslt_test] = generate_train_test_dataset(features,results_all,55);

    figure;
    %plot(results(:,2),results(:,3),'.r');
    plot3(results(:,2),results(:,3),results(:,4),'.r');
    hold on;
    result_nointeraction = results(results_interaction==0,:);
    plot3(result_nointeraction(:,2),result_nointeraction(:,3),result_nointeraction(:,4),'ob');
    xlabel('y');
    ylabel('z');
    zlabel('rot');
    %axis equal;
    axis([-0.1 0.1 0 0.3]);

    %feat_train = feat_train(:,25:64);
    feat_train = feat_train(:,25:64);

    for i=1:size(feat_train,1)
        pred_label(i,rslt_train(i,5)+1) = 1;
    end

    mdl_WLR = LogisticRegression(size(feat_train,2),2,'linearFeatures');
    mdl_WLR.Train(feat_train,pred_label,1e-6);
    interaction_pred = mdl_WLR.computeLikelihood(feat_test(:,25:64),[]);

    for i=1:size(feat_test)
        for j=1:2
            if interaction_pred(i,j) == max(interaction_pred(i,:))
                int_pred_label(i,1) = j;
            end
        end
    end
    
    cnt_correct = 0;
    cnt_incrt = 0;
    for i=1:size(rslt_test)
        if rslt_test(i,5)+1 == int_pred_label(i)
            cnt_correct = cnt_correct +1;
        else
            cnt_incrt = cnt_incrt+1;
        end
    end
    accu(k) = cnt_correct/(cnt_correct+cnt_incrt);
end


