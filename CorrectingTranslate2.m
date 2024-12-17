function correctedwatermarkedimg=CorrectingTranslate2(attackimg,cornerpoint)
%¿´¿´ÍùÄÄ¶ùÒÆ¶¯
[m,n,~]=size(attackimg);
up=0; left=0; 
if cornerpoint(1,1,1)>1  %Íù×óÒÆ
    left=-(cornerpoint(1,1,1)-1);
elseif cornerpoint(1,1,3)<n %ÍùÓÒÒÆ
    left=n-cornerpoint(1,1,3);
end
if cornerpoint(1,2,1)>1 %ÍùÉÏÒÆ
    up=-(cornerpoint(1,2,1)-1);
elseif cornerpoint(1,2,3)<m %ÍùÏÂÒÆ
    up=m-cornerpoint(1,2,3);
end

correctedwatermarkedimg=TranslatingAttack(attackimg,up,left);
% imshow(correctedwatermarkedimg);
end 