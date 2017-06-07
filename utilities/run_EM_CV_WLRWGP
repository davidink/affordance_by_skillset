clear all
close all
load 'features_n_results_interaction2';

results_all = [results results_interaction];

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
plot3(results(:,6),results(:,7),results(:,8),'.g');


num_submodel = 3;

%run EM

%initialize with GMM
options = statset('Display','final');
%gmmodel = fitgmdist(rslt_train(:,2:4),num_submodel,'Options',options);
%submodel_label = cluster(gmmodel, rslt_train(:,2:4));
gmmodel = fitgmdist([rslt_train(:,2:4) rslt_train(:,6:8)],num_submodel,'Options',options);
submodel_label = cluster(gmmodel, [rslt_train(:,2:4) rslt_train(:,6:8)]);

figure;
hold on;
for i=1:num_submodel
    idx = find(submodel_label==i);
    plot3(rslt_train(idx,2),rslt_train(idx,3),rslt_train(idx,4),'.');
end

for i=1:6
    if i <= 3
        run_idx = i+1;
    else
        run_idx = i+2;
    end
    [mdl_WLR{i} mdl_WGP{i}] = computeEM(num_submodel, feat_train, rslt_train(:,run_idx), submodel_label);
    WLR_test_pred{i} = mdl_WLR{i}.computeLikelihood(feat_test,[]);
    for j=1:num_submodel
        [WGP_test_pred{i,j}.pred_mu WGP_test_pred{i,j}.pred_ss] = wgp_pred(mdl_WGP{1,i}{1,j}.hyp, mdl_WGP{1,i}{1,j}.weight, mdl_WGP{1,i}{1,j}.input, mdl_WGP{1,i}{1,j}.target, feat_test);
    end
    
end

% find likelihood of submodel grouping
pred_test_likhd = ones(size(feat_test,1),num_submodel);
for i=1:size(feat_test,1)
    for j=1:num_submodel       
        for k=1:6
            pred_test_likhd(i,j) = pred_test_likhd(i,j) * WLR_test_pred{1,k}(i,j);
        end
    end
end

%find the optimal submodel labels
for i=1:size(feat_test,1)
    for j=1:num_submodel
        if pred_test_likhd(i,j) == max(pred_test_likhd(i,:))
            pred_test_label(i,1) = j;
        end
    end
end

%find the prediction value from selected labels
for i=1:size(feat_test,1)
    for k=1:6
            pred_test_mu(i,k) = WGP_test_pred{k,pred_test_label(i)}.pred_mu(i);
            pred_test_ss(i,k) = WGP_test_pred{k,pred_test_label(i)}.pred_ss(i);
    end    
end

%visualize final predicted submodel labels
figure;
hold on;
for i=1:num_submodel
    idx = find(pred_test_label==i);
    plot3(pred_test_mu(idx,1),pred_test_mu(idx,2),pred_test_mu(idx,3),'.');
    text(pred_test_mu(idx,1),pred_test_mu(idx,2),pred_test_mu(idx,3),[' ' num2str(i)]);
    plot3(rslt_test(idx,2),rslt_test(idx,3),rslt_test(idx,4),'bx');
end
for i=1:size(rslt_test,1)
    plot3([pred_test_mu(i,1);rslt_test(i,2)], [pred_test_mu(i,2);rslt_test(i,3)],[pred_test_mu(i,3);rslt_test(i,4)],'c');
end
xlabel('y');
ylabel('z');
zlabel('rot');
title('Predicted labels for first object');

figure;
hold on;
for i=1:num_submodel
    idx = find(pred_test_label==i);
    plot3(pred_test_mu(idx,4),pred_test_mu(idx,5),pred_test_mu(idx,6),'.');
    text(pred_test_mu(idx,4),pred_test_mu(idx,5),pred_test_mu(idx,6),[' ' num2str(i)]);
    plot3(rslt_test(idx,6),rslt_test(idx,7),rslt_test(idx,8),'bx');
end
for i=1:size(rslt_test,1)
    plot3([pred_test_mu(i,4);rslt_test(i,6)], [pred_test_mu(i,5);rslt_test(i,7)],[pred_test_mu(i,6);rslt_test(i,8)],'c');
end
xlabel('y');
ylabel('z');
zlabel('rot');
title('Predicted labels for second object');

%compute error of prediction for each 6 axes
rslt_test_axis = [rslt_test(:,2:4) rslt_test(:,6:8)];
pred_test_err = rms(rslt_test_axis - pred_test_mu);
for i=1:size(rslt_test_axis,1)
    for k=1:6
        pred_test_likhd(i,k) = log(pdf('Normal',rslt_test_axis(i,k),pred_test_mu(i,k),sqrt(pred_test_ss(i,k))));
        %1/sqrt(2*pi*pred_test_ss(1,1))*exp((rslt_test_axis(1,1)-pred_test_mu(1,1))^2*-1/(2*pred_test_ss(1,1)))
    end
end



% [mdl_WLR mdl_WGP_y] = computeEM(num_submodel, feat_train, rslt_train(:,2), submodel_label); % y
% WLR_y_test_pred = mdl_WLR_y.computeLikelihood(feat_test,[]);
% WLR_y_train_pred = mdl_WLR_y.computeLikelihood(feat_train,[]);
% for i=1:num_submodel
%     [WGP_y{i}.pred_mu WGP_y{i}.pred_ss] = wgp_pred(mdl_WGP_y{i}.hyp, mdl_WGP_y{i}.weight, mdl_WGP_y{i}.input, mdl_WGP_y{i}.target, feat_test);
% end
% 
% for i=1:size(feat_test)
%     for j=1:num_submodel
%         if WLR_y_test_pred(i,j) == max(WLR_y_test_pred(i,:))
%             y_pred_label(i,1) = j;
%         end
%     end
%     y_pred_mu(i,1) = WGP_y{y_pred_label(i)}.pred_mu(i);    
%     
% %     if WLR_y_pred(i,1) >= WLR_y_pred(i,2) 
% %     end
% end
% 
% [mdl_WLR_z mdl_WGP_z] = computeEM(num_submodel, feat_train, rslt_train(:,3), submodel_label); % z
% WLR_z_test_pred = mdl_WLR_z.computeLikelihood(feat_test,[]);
% %WLR_y_train_pred = mdl_WLR_y.computeLikelihood(feat_train,[]);
% for i=1:num_submodel
%     [WGP_z{i}.pred_mu WGP_z{i}.pred_ss] = wgp_pred(mdl_WGP_z{i}.hyp, mdl_WGP_z{i}.weight, mdl_WGP_z{i}.input, mdl_WGP_z{i}.target, feat_test);
% end
% 
% for i=1:size(feat_test)
%     for j=1:num_submodel
%         if WLR_z_test_pred(i,j) == max(WLR_z_test_pred(i,:))
%             z_pred_label(i,1) = j;
%         end
%     end
%     z_pred_mu(i,1) = WGP_z{z_pred_label(i)}.pred_mu(i);    
%     
% %     if WLR_y_pred(i,1) >= WLR_y_pred(i,2) 
% %     end
% end
% 
% [mdl_WLR_r mdl_WGP_r] = computeEM(num_submodel, feat_train, rslt_train(:,4), submodel_label); % rot
% WLR_r_test_pred = mdl_WLR_r.computeLikelihood(feat_test,[]);
% %WLR_y_train_pred = mdl_WLR_y.computeLikelihood(feat_train,[]);
% for i=1:num_submodel
%     [WGP_r{i}.pred_mu WGP_r{i}.pred_ss] = wgp_pred(mdl_WGP_r{i}.hyp, mdl_WGP_r{i}.weight, mdl_WGP_r{i}.input, mdl_WGP_r{i}.target, feat_test);
% end
% 
% for i=1:size(feat_test)
%     for j=1:num_submodel
%         if WLR_r_test_pred(i,j) == max(WLR_r_test_pred(i,:))
%             r_pred_label(i,1) = j;
%         end
%     end
%     r_pred_mu(i,1) = WGP_r{r_pred_label(i)}.pred_mu(i);    
%     
% %     if WLR_y_pred(i,1) >= WLR_y_pred(i,2) 
% %     end
% end

% figure;
% hold on;
% plot3(y_pred_mu,z_pred_mu,r_pred_mu,'.r');
% plot3(rslt_test(:,2),rslt_test(:,3),rslt_test(:,4),'ob');
% for i=1:size(rslt_test,1)
%     plot3([y_pred_mu(i);rslt_test(i,2)], [z_pred_mu(i);rslt_test(i,3)],[r_pred_mu(i);rslt_test(i,4)],'c');
% end
% xlabel('y');
% ylabel('z');
% zlabel('rot');
% %axis equal;
% axis([-0.1 0.1 0 0.3]);
% pred_err = [rms(y_pred_mu - rslt_test(:,2)) rms(z_pred_mu - rslt_test(:,3)) rms(r_pred_mu - rslt_test(:,4))];

