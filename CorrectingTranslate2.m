function correctedwatermarkedimg=CorrectingTranslate2(attackimg,cornerpoint)
%�������Ķ��ƶ�
[m,n,~]=size(attackimg);
up=0; left=0; 
if cornerpoint(1,1,1)>1  %������
    left=-(cornerpoint(1,1,1)-1);
elseif cornerpoint(1,1,3)<n %������
    left=n-cornerpoint(1,1,3);
end
if cornerpoint(1,2,1)>1 %������
    up=-(cornerpoint(1,2,1)-1);
elseif cornerpoint(1,2,3)<m %������
    up=m-cornerpoint(1,2,3);
end

correctedwatermarkedimg=TranslatingAttack(attackimg,up,left);
% imshow(correctedwatermarkedimg);
end 