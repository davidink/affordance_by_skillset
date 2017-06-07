clear all
close all

num_submodel = 4;

colorOrder = [0.8 0 0; 0 0.8 0; 0 0 0.8; 0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0.8 0.8 0.8];

for i=1:num_submodel % submodel
    for j=1:8 % feat set
        pred_err = [];
        pred_err_ll = [];
        data_loc = ['rslts2/sub_' num2str(i) '_featset' num2str(j) '.mat'];                
        load(data_loc);
        pred_err = [pred_err;pred_test_err];            
        pred_err_ll = [pred_err_ll;pred_test_ll_mean];            
        
        pred_err_sub{i,j} = pred_err;
        pred_err_ll_sub{i,j} = pred_err_ll;
    end
end

for i=1:num_submodel
    for j=1:8
        avg_pred_err{i}(j,:) = mean(pred_err_sub{i,j});
        std_pred_err{i}(j,:) = std(pred_err_sub{i,j});
        avg_pred_err_ll{i}(j,:) = mean(pred_err_ll_sub{i,j});
        std_pred_err_ll{i}(j,:) = std(pred_err_ll_sub{i,j});
    end
end

for i=1:6 % num_axis
    figure
    hold on
    avg_data = [];
    std_data = [];
    for j=1:num_submodel %num_submodel
        avg_data = [avg_data avg_pred_err{j}(:,i)]; 
        std_data = [std_data std_pred_err{j}(:,i)];
    end
    hb = bar(1:num_submodel,avg_data'); %num_submodel
    % For each set of bars, find the centers of the bars, and write error bars
    pause(0.1); %pause allows the figure to be created
    for ib = 1:numel(hb)
        %XData property is the tick labels/group centers; XOffset is the offset
        %of each distinct group
        xData = hb(ib).XData+hb(ib).XOffset;
        errorbar(xData,avg_data(ib,:),std_data(ib,:),'k.')
    end

%     hb(1).FaceColor = [0.8 0 0];
%     hb(2).FaceColor = [0.0 0.8 0];
%     hb(3).FaceColor = [0.0 0 0.8];
%     hb(4).FaceColor = [0.5 0.5 0];
%     hb(5).FaceColor = [0.5 0 0.5];
%     hb(6).FaceColor = [0.0 0.5 0.5];
%     hb(7).FaceColor = [0.0 0.3 0.3];
%     hb(8).FaceColor = [0.7 0.7 0.1];

    Labels = {'sub_1', 'sub_2','sub_3','sub_4', 'sub_5','sub_6'};
    set(gca, 'XTick', 1:6, 'XTickLabel', Labels);
    
    legend('Object(A) / A','A+Contact(B) / A+B','Context(C) / C','A+B+C / A+B+C','A+B / C','A+B / A+B+C','C / A+B','C / A+B+C','Location','southwest');
    if i==1 
        title('RMSE of Prediction Models: \Delta y_{p1}');
    end
    if i==2 
        title('RMSE of Prediction Models: \Delta z_{p1}');
    end
    if i==3 
        title('RMSE of Prediction Models: \Delta \theta_{p1}');
    end
    if i==4 
        title('RMSE of Prediction Models: \Delta y_{p2}');
    end
    if i==5 
        title('RMSE of Prediction Models: \Delta z_{p2}');
    end
    if i==6 
        title('RMSE of Prediction Models: \Delta \theta_{p2}');
    end
    
    print(['figs/rms_axis_' num2str(i)],'-depsc');

%     y_lim = get(gca,'ylim');
%     ylabel('meter');
%     box off
%     % Create second Y axes on the right.
%     a2 = axes('YAxisLocation', 'Right');
%     % Hide second plot.
%     set(a2, 'color', 'none');
%     set(a2, 'XTick', []);
%     % Set scala for second Y.
%     set(a2, 'YLim', [y_lim(1)*theta_scale y_lim(2)*theta_scale]);
%     ylabel('radian');
end

% for i=1:6
%     figure
%     hold on
%     avg_data = [];
%     std_data = [];
%     for j=1:6
%         avg_data = [avg_data avg_pred_err_ll{j}(:,i)]; 
%         std_data = [std_data std_pred_err_ll{j}(:,i)];
%     end
%     hb = bar(1:6,avg_data');
%     % For each set of bars, find the centers of the bars, and write error bars
%     pause(0.1); %pause allows the figure to be created
%     for ib = 1:numel(hb)
%         %XData property is the tick labels/group centers; XOffset is the offset
%         %of each distinct group
%         xData = hb(ib).XData+hb(ib).XOffset;
%         errorbar(xData,avg_data(ib,:),std_data(ib,:),'k.')
%     end
% 
%     % hb(1).FaceColor = [1 0 0];
%     % hb(2).FaceColor = [0.0 1 0];
%     % hb(3).FaceColor = [0.0 0 1];
%     % hb(4).FaceColor = [0.5 0.5 0];
%     % hb(5).FaceColor = [0.5 0 0.5];
%     % hb(6).FaceColor = [0.0 0.5 0.5];
% 
%     Labels = {'sub_1', 'sub_2','sub_3','sub_4', 'sub_5','sub_6'};
%     set(gca, 'XTick', 1:6, 'XTickLabel', Labels);
% 
%     legend('Object/Object','Obj+Cnct/Obj+Cnct','Obj+Cnct+Cntxt/Obj+Cnct+Cntxt','Cntxt/Obj+Cnct','Cntxt/Obj+Cnct+Cntxt','Cntxt(discrete)/Obj+Cnct','Obj+Cnct/Cntxt','Obj+Cnct/Obj+Cnct+Cntxt','Location','northwest');
%     title('RMSE of Forward Affordance Prediction Models');
% 
%     y_lim = get(gca,'ylim');
%     ylabel('meter');
%     box off
%     % Create second Y axes on the right.
%     a2 = axes('YAxisLocation', 'Right');
%     % Hide second plot.
%     set(a2, 'color', 'none');
%     set(a2, 'XTick', []);
%     % Set scala for second Y.
%     set(a2, 'YLim', [y_lim(1)*theta_scale y_lim(2)*theta_scale]);
%     ylabel('radian');
% end

% for i=1:3
%     figure
%     hold on
%     hb = bar(1:6,avg_pred_err{i}');
%     % For each set of bars, find the centers of the bars, and write error bars
%     pause(0.1); %pause allows the figure to be created
%     for ib = 1:numel(hb)
%         %XData property is the tick labels/group centers; XOffset is the offset
%         %of each distinct group
%         xData = hb(ib).XData+hb(ib).XOffset;
%         errorbar(xData,avg_pred_err{i}(ib,:),std_pred_err{i}(ib,:),'k.')
%     end
% 
%     % hb(1).FaceColor = [1 0 0];
%     % hb(2).FaceColor = [0.0 1 0];
%     % hb(3).FaceColor = [0.0 0 1];
%     % hb(4).FaceColor = [0.5 0.5 0];
%     % hb(5).FaceColor = [0.5 0 0.5];
%     % hb(6).FaceColor = [0.0 0.5 0.5];
% 
%     Labels = {'y_p', 'z_p','theta','y_p', 'z_p','theta'};
%     set(gca, 'XTick', 1:6, 'XTickLabel', Labels);
% 
%     legend('Object','Obj+Contact','Obj+Contact+Context','Context/Obj+Cont','Obj+Cont/Context','Context(discrete)/Obj+Cont','Context/Ojb+Cont','Obj+Cont/Obj+Cont+Context','Location','northwest');
%     title('RMSE of Forward Affordance Prediction Models');
% 
%     y_lim = get(gca,'ylim');
%     ylabel('meter');
%     box off
%     % Create second Y axes on the right.
%     a2 = axes('YAxisLocation', 'Right');
%     % Hide second plot.
%     set(a2, 'color', 'none');
%     set(a2, 'XTick', []);
%     % Set scala for second Y.
%     set(a2, 'YLim', [y_lim(1)*theta_scale y_lim(2)*theta_scale]);
%     ylabel('radian');
% end
