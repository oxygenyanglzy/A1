function dy=odefun(t,y)
    dy=zeros(3,1);
   % dy(1)=16*(y(2)-y(1));
   dy(1)=26*(y(2)-y(1));
%     dy(2)=y(1)*(46-y(3))-y(2);
   dy(2)=y(1)*(47-y(3))-y(2);
   dy(3)=y(1)*y(2)-4*y(3);
end
