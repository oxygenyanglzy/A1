function [watermarked] = Imagecorrection(attackimg)
%IMAGECORRECTION 此处显示有关此函数的摘要
%   此处显示详细说明
[cornerpoint,attacktype]=LocatingVertices4(attackimg);
if attacktype %平移
    watermarked=CorrectingTranslate2(attackimg,cornerpoint);
else %其他情况
    src1=[cornerpoint(:,:,1);cornerpoint(:,:,2);cornerpoint(:,:,3);cornerpoint(:,:,4)];%改，根据不同角度确定对应的是cornerpoint的哪个点
    pt_M=512;
    tar1=[1,1;pt_M,1;pt_M,pt_M;1,pt_M];
    TForm1 = fitgeotrans(src1,tar1,'Projective'); %Projective
    attackimg_pt = imwarp(attackimg, TForm1);
    %     imwrite(attackimg_pt,'attackimg_pt.bmp');
    %     imshow(attackimg_pt);
    % % 去黑边及消除锯齿效应
    watermarked=RemovingBlackEdges4(attackimg_pt,512,512); % 去黑边
    %     imwrite(uint8(watermarked),'removeblackedge_F16.bmp');
    %figure(3),subplot(122),imshow(watermarked),title('remove black matte image');
    watermarked=sawtoothProcess(double(watermarked)); %消除锯齿效应
    %     imwrite(uint8(watermarked),'reducesawtootheffect.bmp');
    %figure(4),subplot(121),imshow(uint8(watermarked)),title('remove sawtooth image');
end

