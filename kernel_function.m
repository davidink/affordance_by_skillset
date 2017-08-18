function kernel = kernel_function(XN,XM,theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%      val = (exp(theta(2))^2)*exp(-(pdist2(XN,XM).^2)/(2*exp(theta(1))^2));
%     disp(['function called: ' num2str(theta(1)) ' ' num2str(theta(2))]);
    for i=1:size(XN,1)
        for j=i:size(XM,1)
            kernel(i,j) = compute_kernel_custom(XN(i,:),XM(j,:),theta);
            kernel(j,i) = kernel(i,j);
        end
    end
%     val = kernel;
end

