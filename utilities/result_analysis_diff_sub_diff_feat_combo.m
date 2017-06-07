clear all
close all

for i=1:6 % submodel
    for j=1:8 % feat set
        pred_err = [];
        pred_err_ll = [];
        for k=1:3 % num of cv
            data_loc = ['rslts/sub_' num2str(i) '_featset' num2str(j) '_cv' num2str(k) '.mat'];                
            load(data_loc);
            pred_err = [pred_err;pred_test_err];            
            pred_err_ll = [pred_err_ll;pred_test_ll_mean];            
        end
        pred_err_sub{i,j} = pred_err;
        pred_err_ll_sub{i,j} = pred_err_ll;
    end
end

%Swap set 7 and set 5
for i=1:6
    temp = pred_err_sub{i,5};
    pred_err_sub{i,5} = pred_err_sub{i,7};
    pred_err_sub{i,7} = temp;
    temp = pred_err_ll_sub{i,5};
    pred_err_ll_sub{i,5} = pred_err_ll_sub{i,7};
    pred_err_ll_sub{i,7} = temp;
    
end

for i=1:6
    for j=1:8
        avg_pred_err{i}(j,:) = mean(pred_err_sub{i,j});
        std_pred_err{i}(j,:) = std(pred_err_sub{i,j});
        avg_pred_err_ll{i}(j,:) = mean(pred_err_ll_sub{i,j});
        std_pred_err_ll{i}(j,:) = std(pred_err_ll_sub{i,j});
    end
end

% visualize
theta_scale = 10;
for i=1:6
    avg_pred_err{i}(:,3) = avg_pred_err{i}(:,3)/theta_scale;
    std_pred_err{i}(:,3) = std_pred_err{i}(:,3)/theta_scale; 
    avg_pred_err{i}(:,6) = avg_pred_err{i}(:,6)/theta_scale;
    std_pred_err{i}(:,6) = std_pred_err{i}(:,6)/theta_scale; 
end

for i=1:6
    figure
    hold on
    avg_data = [];
    std_data = [];
    for j=1:6
        avg_data = [avg_data avg_pred_err{j}(:,i)]; 
        std_data = [std_data std_pred_err{j}(:,i)];
    end
    hb = bar(1:6,avg_data');
    % For each set of bars, find the centers of the bars, and write error bars
    pause(0.1); %pause allows the figure to be created
    for ib = 1:numel(hb)
        %XData property is the tick labels/group centers; XOffset is the offset
        %of each distinct group
        xData = hb(ib).XData+hb(ib).XOffset;
        errorbar(xData,avg_data(ib,:),std_data(ib,:),'k.')
    end

    % hb(1).FaceColor = [1 0 0];
    % hb(2).FaceColor = [0.0 1 0];
    % hb(3).FaceColor = [0.0 0 1];
    % hb(4).FaceColor = [0.5 0.5 0];
    % hb(5).FaceColor = [0.5 0 0.5];
    % hb(6).FaceColor = [0.0 0.5 0.5];

    Labels = {'sub_1', 'sub_2','sub_3','sub_4', 'sub_5','sub_6'};
    set(gca, 'XTick', 1:6, 'XTickLabel', Labels);

    legend('Object/Object','Obj+Cnct/Obj+Cnct','Obj+Cnct+Cntxt/Obj+Cnct+Cntxt','Cntxt/Obj+Cnct','Cntxt/Obj+Cnct+Cntxt','Cntxt(discrete)/Obj+Cnct','Obj+Cnct/Cntxt','Obj+Cnct/Obj+Cnct+Cntxt','Location','northwest');
    title('RMSE of Forward Affordance Prediction Models');

    y_lim = get(gca,'ylim');
    ylabel('meter');
    box off
    % Create second Y axes on the right.
    a2 = axes('YAxisLocation', 'Right');
    % Hide second plot.
    set(a2, 'color', 'none');
    set(a2, 'XTick', []);
    % Set scala for second Y.
    set(a2, 'YLim', [y_lim(1)*theta_scale y_lim(2)*theta_scale]);
    ylabel('radian');
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
