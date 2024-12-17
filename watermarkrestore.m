function extractlevelwatermark=watermarkrestore(wateronedim)
for level=1:3
     B=wateronedim(level,:);
     lB=size(B,2);
      for i=1:lB/8
          extract(1,i)=bin2dec(B(1,(i-1)*8+1:(i-1)*8+8));
      end
      extracts=reshape(extract,32,32);
      extractlevelwatermark(:,:,level)=extracts';
 end