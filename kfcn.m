function val = kfcn()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
     val = (exp(theta(2))^2)*exp(-(pdist2(XN,XM).^2)/(2*exp(theta(1))^2));
end

