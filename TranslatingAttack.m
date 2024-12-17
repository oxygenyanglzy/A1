function img=TranslatingAttack(watermarkedimg,x,y)
%----------平移变换---------------2019.12.06完成
I=watermarkedimg;
se = translate(strel(1), [x y]);
img = imdilate(I,se);
% imwrite(uint8(img),['Lena-Translation(',num2str(x),',',num2str(y),').jpg']);
% % imwrite(uint8(I1),'ImageDatebases\Lena-Translation(-20,0).jpg');
% figure(2),imshow(uint8(img),'InitialMagnification','fit');%图像过大时，使用该自适应函数显示
end