load TrainingSamplesDCT_subsets_8
load prior_1
load Alpha

dataset1 = D3_BG;
dataset2 = D3_FG;
D1mean_BG = zeros(1,64);
D1mean_FG = zeros(1,64);
BGn = size(dataset1,1);
FGn = size(dataset2,1);

mask = imread('cheetah_mask.bmp');
mask = double(mask)/255;
%ML solution:
%realmean,not used--------------------------------------------

for i =1:BGn;
    D1mean_BG = D1mean_BG+ dataset1(i,:);
end
for i =1:FGn;
    D1mean_FG = D1mean_FG+ dataset2(i,:);
end
D1mean_BG = D1mean_BG/BGn;
D1mean_FG = D1mean_FG/FGn;
D1variance_BG = zeros(64,64);
D1variance_FG = zeros(64,64);

for i =1:BGn
    D1variance_BG = D1variance_BG+(dataset1(i,:)-D1mean_BG)*(dataset1(i,:)-D1mean_BG)';
end

for i =1:FGn
     D1variance_FG = D1variance_FG+(dataset2(i,:)-D1mean_FG)*(dataset2(i,:)-D1mean_FG)';
end
%realcov,not used--------------------------------------------

D1variance_FG = D1variance_FG/FGn;
D1variance_BG = D1variance_BG/BGn;
Listcovariance_BG = {};
Listmean_BG = {};
Listcovariance_FG = {};
Listmean_FG = {};
covariance = zeros(64,64);
for i = 1:size(alpha,2)
    for j = 1:64
        covariance(j,j) = alpha(i)*W0(j);
    end
    Listcovariance_BG{i} =  D1variance_BG;
    Listcovariance_FG{i} =  D1variance_FG;
    
    Listmean_BG{i} = D1mean_BG'; 
    Listmean_FG{i} = D1mean_FG'; 
end



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
    imagefeature(k,:) = zigblocks{k};
end
%----------------------------------------------------------------------------------------------------------------------
graphaftertotal = {}; 
counter = 1;
for i = 1:size(alpha,2)
    graphafter = ones(1,totalblock); 
    for j = 1:totalblock
        testfeature = imagefeature(i,:);
        probability1 = Listmean_BG{1,i}'*Listcovariance_BG{1,i}*imagefeature(j,:)'-0.5* Listmean_BG{1,i}'*Listcovariance_BG{1,i}* Listmean_BG{1,i};
        probability2 = Listmean_FG{1,i}'*Listcovariance_FG{1,i}*imagefeature(j,:)'-0.5* Listmean_FG{1,i}'*Listcovariance_FG{1,i}* Listmean_FG{1,i};
        if(sum(probability1)>=sum(probability2))
            graphafter(1,j) = 0;
        end
    end
    graphaftertotal{counter} = graphafter;
    counter= counter+1;
end
hohoho = mask(1:247,1:262);
difference = zeros(1,size(graphaftertotal,2));
for i = 1:size(graphaftertotal,2)
    difference(1,i) = sum(abs(graphaftertotal{1,i}-reshape(hohoho,[1,247*262])))/(247*262);
end
min = 999999999;
minposition = 1;
for i = 1:size(graphaftertotal,2)
    if (difference(1,i) <= min)
        min = difference(1,i);
        minposition=i;
    end
end
    
graphafter = reshape(graphaftertotal{1,minposition},[yint,xint]);
%graphafter = rot90(graphafter,1);
pcolor(graphafter)
    