load TrainingSamplesDCT_subsets_8
load TrainingSamplesDCT_8_new
load prior_1
load prior_2
load Alpha

grass = TrainsampleDCT_BG;
cheetah = TrainsampleDCT_FG;
meangrass = mean(grass);
meancheetah = mean(cheetah);
hahaha1 = cov(grass);
hahaha2 = cov(cheetah);
u = meangrass';

image = imread('cheetah.bmp');
original1 = image;
image = double(image)/255;
blocks = {};
imagerow = size(image,1);
imagecol = size(image,2);
xint =size(image,2)-8;
yint =size(image,1)-8;
ooooo=1
for i=1:xint
    for j = 1:yint
        block = image(j:j+7,i:i+7);%(8*(i-1)+1:8*i,8*(j-1)+1:8*j);
        blocks{ooooo} = block;
        ooooo = ooooo+1;
    end    
end
totalblock = ooooo-1;
for i = 1:totalblock
    blocks{1,i} = dct2(blocks{1,i}); %dct2(reshape(( double(blocks{i})),[1,(size(blocks{i},1)*(size(blocks{i},2)))]));
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
graphafter = zeros(1,totalblock); 
for i = 1:totalblock
    testfeature = imagefeature(i,:);
    probability1 = meangrass*hahaha1*imagefeature(i,:)'-0.5*meangrass*hahaha1*meangrass';
    probability2 = meancheetah*hahaha2*imagefeature(i,:)'-0.5*meancheetah*hahaha2*meancheetah';
    if(sum(probability1)>=sum(probability2))
        graphafter(1,i) = 1;
    end
end

graphafter = reshape(graphafter,[yint,xint]);
%graphafter = rot90(graphafter,1);
pcolor(mask)
    