function A_encry= gongg(  )
%GONGG �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
initial_value = [0.001,0.005,0.002];
[~,~,A] = Lorenz(initial_value);
% ������ʼʱ�̵����ݶ��Ƚ�С������޳���ʼʱ�̵�����
A = A(100:end);
% ��һ����������
A_normal = (A-min(A))./(max(A)-min(A));
% ���Ʒ�ΧΪ[0,255]
A_uint8 = uint8(A_normal*255);
A_encry = reshape(A_uint8(1:32*32),32,32);
end

