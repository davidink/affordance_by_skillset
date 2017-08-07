function result=GaussBhatKern(x1,prec1,x2, prec2)
    rho=0.5;
    D=length(x1);
    eps = 1e-16;
    prec1 = prec1 + eps*eye(size(prec1,1));
    prec2 = prec2 + eps*eye(size(prec1,1));
    sigmaC=inv(prec1+prec2);
    precC=(prec1+prec2);
    xC=prec1*x1+prec2*x2;
    result=(sqrt((2*pi)^((1-2*rho)*D/2)))*(rho^(-0.5*D))*(sqrt(det(sigmaC)))*(sqrt(det(prec1)^rho))*(sqrt(det(prec2)^rho))*exp(-rho*0.5*(x1'*prec1*x1+x2'*prec2*x2-xC'*sigmaC*xC));
end