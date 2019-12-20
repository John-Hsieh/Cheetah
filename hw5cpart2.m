
%Start processing the graph data
image = imread('cheetah.bmp');
original1 = image;
image = double(image)/255;
blocks = {};
imagerow = size(image,1);
imagecol = size(image,2);
xint =size(image,2)-8;
yint =size(image,1)-8;
ooooo=1;
for i=1:xint
    for j = 1:yint
        block = image(j:j+7,i:i+7);%(8*(i-1)+1:8*i,8*(j-1)+1:8*j);
        blocks{ooooo} = block;
        ooooo = ooooo+1;
    end    
end
totalblock = ooooo-1;
for i = 1:totalblock
    blocks{1,i} = dct2(blocks{1,i}); 
end

%----------------------------------------------------------------------------------------------------------------------

zigblocks = {};
for i=1:totalblock%对每个进行zigzag
    counter = 1;
    zigsingle = [1,(size(blocks{i},1)*size(blocks{i},2))];
    %inner loop, deal with first part
    for j = 1:size(blocks{i},2)
       if rem(j,2)==0
            for k = 0:j-1
                zigsingle(counter)=blocks{i}(1+k,j-k);
                counter = counter + 1;
            end
       else
           for k = 0:j-1
                zigsingle(counter)=blocks{i}(j-k,1+k);
                counter = counter + 1;
            end
        end
    end
    for j = 2:size(blocks{i},1)%j=start of y
        loopstop = j+size(blocks{i},2);
        if loopstop>size(blocks{i},1)
            loopstop = size(blocks{i},1);
        end
        if rem(j,2)==0
            for k = j:loopstop%K=start of x
                zigsingle(counter)=blocks{i}(size(blocks{i},2)-k+j,k);
                counter = counter+1;
            end
        else
            for k = j:loopstop
                zigsingle(counter)=blocks{i}(k,size(blocks{i},2)-k+j);
                counter = counter+1;
            end
        end
    end
    zigsingle = reshape(zigsingle,[1,(size(zigsingle,1)*size(zigsingle,2))]);
    zigblocks{i}=zigsingle;
end

%----------------------------------------------------------------------------------------------------------------------
imagefeature = ones(size(zigblocks,2),64);%ones(size(zigblocks,2),8);
for k = 1:size(zigblocks,2)
    imagefeature(k,:) = zigblocks{k};%(originbig); %- max(blocks{k})*ones(size(blocks{k}));%are we trying to minus from the whole graph?
end
%----------------------------------------------------------------------------------------------------------------------
%Image features extracted.

graphaftertotal = {}; 
counter = 1;
graphafter = ones(1,totalblock); 

errorstorer = zeros(5,feanum);
counter = 1;
for fea = 1:feanum
        component = feaDim(fea);
        graphafter = ones(1,totalblock); 
        for j = 1:totalblock
            probability1 = 0;
            probability2 = 0;
            for k = 1:component
                probability1 = probability1+mvnpdf(imagefeature(j,:),FGmeanhold{fea}(k,:),FGstdhold{fea}(:,:,k))*FGweighthold{fea}(1,k);
                probability2 = probability2+mvnpdf(imagefeature(j,:),BGmeanhold{fea}(k,:),BGstdhold{fea}(:,:,k))*BGweighthold{fea}(1,k);
            end
            if(sum(probability1)>=sum(probability2))
                graphafter(1,j) = 0;
            end
        end
        hohoho = mask(1:247,1:262);
        errorstorer(1,fea) = sum(abs(graphafter-reshape(hohoho,[1,247*262])))/(247*262);
        counter = counter+1;
end
hahaha1 = ones(1,25) - errorstorer;%the error rate
graphafter = reshape(graphafter,[yint,xint]);
%graphafter = rot90(graphafter,1);
pcolor(graphafter)