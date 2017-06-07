load 'rslts/sub_1_gmm_init_1.mat'
pred_err_all = [pred_err];
load 'rslts/sub_2_gmm_init_1.mat'
pred_err_all = [pred_err_all;pred_err];
load 'rslts/sub_3_gmm_init_1.mat'
pred_err_all = [pred_err_all;pred_err];
load 'rslts/sub_4_gmm_init_1.mat'
pred_err_all = [pred_err_all;pred_err];
load 'rslts/sub_5_gmm_init_1.mat'
pred_err_all = [pred_err_all;pred_err];
load 'rslts/sub_6_gmm_init_1.mat'
pred_err_all = [pred_err_all;pred_err];

load 'rslts/LR_interaction_pred_feat_all.mat'
accu_all(1,:) = accu;
load 'rslts/LR_interaction_pred_feat_cont.mat'
accu_all(2,:) = accu;