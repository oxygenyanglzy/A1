function [attackimg,aname]=GeometricAttack(watermarkedimg)
para=20;attackimg=RotatingAttack(watermarkedimg,para); name='rotate';%Ğı×ª
% para=1;tar1=[50,1;512,60;450,512;1,400]; attackimg=SlantAttack(watermarkedimg,tar1); name='slant'; %ÇãĞ±
% para=2;tar2=[1,120;350,1;512,420;100,512]; attackimg=SlantAttack(watermarkedimg,tar2); name='slant';%ÇãĞ±
% para=0.09;attackimg=AffineAttack(watermarkedimg,1,para,para); name='affine';% ·ÂÉä
% para=0.1;attackimg=AffineAttack(watermarkedimg,2,para,0); name='y-affine';% Ë®Æ½·ÂÉä ±ßÔµ¼ì²âºóÕÒ²»µ½µã
% para=0.1;attackimg=AffineAttack(watermarkedimg,3,0,para); name='x-affine';% ´¹Ö±·ÂÉä
% para=-1;attackimg=TranslatingAttack(watermarkedimg,-35,-25); name='translate';% Æ½ÒÆ ±ßÔµ¼ì²âºóÕÒ²»µ½µã
% imshow(attackimg);
aname=[name,num2str(para)]; %rotate30
