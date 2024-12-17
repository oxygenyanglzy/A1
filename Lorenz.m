
%LORENZ 此处显示有关此函数的摘要
%   此处显示详细说明

function [component_one,component_two,component_three]=Lorenz(initial_value)%1E-17
    tspan=[20,100];%时间截断范围
    [~,Y0]=ode45(@odefun,tspan,initial_value);
    component_one=Y0(:,1)';
    component_two=Y0(:,2)';
    component_three=Y0(:,3)';
%     figure(1);
%     plot3(Y0(:,1),Y0(:,2),Y0(:,3));
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     title('Lorenz混沌系统相空间图');
%     grid on;
end

