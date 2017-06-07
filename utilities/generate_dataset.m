clear all
close all

x1 = [-1:0.02:1];
x2 = [-1:0.02:1];

[xx1, xx2] = meshgrid(x1,x2);

for i=1:size(x1,2)
    for j=1:size(x2,2)
        if x1(i) >= x2(j)
            y(i,j) = sin(2*pi*x1(i));
        else
            y(i,j) = cos(2*pi*x2(j));
        end
    end
end

y = y + 0.1*randn(size(y,1));
figure;
surf(xx1,xx2,y);

X_all = [];
Y_all = [];
for i=1:size(xx1,1)
    for j=1:size(xx2,1)
        X_all = [X_all;xx1(i,j) xx2(i,j)];
        Y_all = [Y_all;y(i,j)];
    end
end

figure
plot3(X_all(:,1),X_all(:,2),Y_all,'.r');
xlabel('x1');
ylabel('x2');

sample_cnt = 1000;

idx = randperm(size(Y_all,1),sample_cnt);
X = X_all(idx',:);
Y = Y_all(idx');
% 
% X = [x1(idx1)' x2(idx2)'];
% 
% for i=1:sample_cnt
%     Y(i,1) = y(idx1(i),idx2(i));
% end

figure

figure
plot3(X(:,1), X(:,2), Y, '.r'); 

save('dataset1.mat','X','Y','X_all','Y_all');