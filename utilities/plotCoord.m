function plotCoord(pos,R,shift)
hold on
plot3(pos(1),pos(2),pos(3),'.k');
plot3([pos(1) pos(1)+shift*R(1,1)],[pos(2) pos(2)+shift*R(2,1)],[pos(3) pos(3)+shift*R(3,1)],'r');
plot3([pos(1) pos(1)+shift*R(1,2)],[pos(2) pos(2)+shift*R(2,2)],[pos(3) pos(3)+shift*R(3,2)],'g');
plot3([pos(1) pos(1)+shift*R(1,3)],[pos(2) pos(2)+shift*R(2,3)],[pos(3) pos(3)+shift*R(3,3)],'b');
end