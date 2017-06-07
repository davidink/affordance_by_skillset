function [euler,R]=qGetE(inR)
% R=qGetR(quat);
R=inR;
Z=R(:,3);
Z(2)=0;

if norm(Z)==0
    euler(1)=0;
else
    Z=Z/norm(Z);
    euler(1)=atan2(Z(1),Z(3));
end

e1R=[[cos(euler(1)) 0 sin(euler(1))];[0 1 0];[-sin(euler(1)) 0 cos(euler(1))]];
R=e1R'*R;


Z=R(:,3);
Z(1)=0;

if norm(Z)==0
    euler(1)=0;
else
    Z=Z/norm(Z);
    euler(2)=atan2(-Z(2),Z(3));
end

e2R=[[1 0 0];[0 cos(euler(2)) -sin(euler(2))];[0 sin(euler(2)) cos(euler(2))]];
R=e2R'*R;

euler(3)=atan2(R(2,1),R(1,1));

e3R=[[cos(euler(3)) -sin(euler(3)) 0];[sin(euler(3)) cos(euler(3)) 0];[0 0 1]];
R=e3R'*R;

end