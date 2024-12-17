function value = performance_value(inputArg1,inputArg2)
for level=1:3
    ws=double(watermarked(:,:,level));
    coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
    wa=coeffs{1,1};%%su3
    for  counter=1:lenw
        i=position(1,counter);
        j=position(2,counter);
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u2,s2]=schur(block);
        Rext=fix(s2(1,1)/Qfac);
        ext=mod(Rext,2);
        if ext==0
            watermark='0';
        else
            watermark='1';
        end
        ExWater(1,counter)=watermark;
    end
    wateronedim(level,:)=ExWater;
    
end
%水印恢复
%extractlevelwatermark=watermarkrestore(wateronedim);
for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
end
%toc;
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
end

