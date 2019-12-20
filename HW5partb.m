load TrainingSamplesDCT_8_new
mask = imread('cheetah_mask.bmp');
mask = double(mask)/255;
cheetah_image = imread('cheetah.bmp');

%prior:
[FG_row,~] = size(TrainsampleDCT_FG);
[BG_row,~] = size(TrainsampleDCT_BG);
FG_prior = (FG_row/(FG_row+BG_row));
BG_prior = (BG_row/(FG_row+BG_row));
%general variables
testDim = [1,2,4,8,16,24,32,40,48,56,64];
feaDim = [1,2,4,8,16,24,32,40,48,56,64];
feanum = size(feaDim,2);
component = 8;
mixture = 5;    %5groups of values
gaussiannum = 8;
featuredim = 64;
itercount = 100;
%each should have 8 dimension = component, 


FGtotalweight = ones(1,8,5);
BGtotalweight = ones(1,8,5);

FGtotalmean = ones(component,featuredim,5);
BGtotalmean = ones(component,featuredim,5);

FGtotalstd = zeros(64,64,gaussiannum,5);
BGtotalstd = zeros(64,64,gaussiannum,5);
for classes = 1:5
    FGmean = TrainsampleDCT_FG(randi(FG_row,1,8),:);%5mixtures, each = 64*8
    BGmean = TrainsampleDCT_BG(randi(BG_row,1,8),:);

    FGstd = zeros(64,64,gaussiannum); %need to be square covariance matrix
    BGstd = zeros(64,64,gaussiannum);

    for i = 1:8
        FGstd(:,:,i) = diag(rand(1,gaussiannum*gaussiannum));
        BGstd(:,:,i) = diag(rand(1,gaussiannum*gaussiannum));
    end
    FGweight = ones(1,8);
    FGweight = FGweight/sum(FGweight);
    BGweight = ones(1,8);
    BGweight = BGweight/sum(BGweight);

    FGhmatrix = rand(FG_row,gaussiannum);%zeros(FG_row,component);
    BGhmatrix = rand(BG_row,gaussiannum);

    FGposs = zeros(size(FGhmatrix,1),component);
    BGposs = zeros(size(BGhmatrix,1),component);
    %TrainsampleDCT_FG = TrainsampleDCT_FG(:,component);
    for iter = 1:itercount
        for origin = 1:component
            FGposs(:,origin) = mvnpdf(TrainsampleDCT_FG,FGmean(origin,:),FGstd(:,:,origin));
        end

        for i = 1:gaussiannum
            FGhmatrix(:,i) = mvnpdf(TrainsampleDCT_FG,FGmean(i,:),FGstd(:,:,i))*FGweight(i);
            FGhmatrix(:,i) = FGhmatrix(:,i)./(FGposs*FGweight');

            FGtotal = sum(FGhmatrix(:,i));
            FGweight(i) = FGtotal/size(FGhmatrix,1);
            FGmean(i,:) = sum(repmat(FGhmatrix(:,i),1,featuredim).*TrainsampleDCT_FG)'/FGtotal;

            FGstdUP = zeros(featuredim,featuredim);

            for k = 1:FG_row
               difference = TrainsampleDCT_FG(k,:)-FGmean(i,:);
               FGstdUP = FGstdUP + FGhmatrix(k,i)*diag(diag(difference'*difference));
            end
            FGstd(:,:,i) = FGstdUP/FGtotal;
        end
    end
    
    
    for iter = 1:itercount
        for origin = 1:component
            BGposs(:,origin) = mvnpdf(TrainsampleDCT_BG,BGmean(origin,:),BGstd(:,:,origin));
        end

        for i = 1:gaussiannum
            BGhmatrix(:,i) = mvnpdf(TrainsampleDCT_BG,BGmean(i,:),BGstd(:,:,i))*BGweight(i);
            BGhmatrix(:,i) = BGhmatrix(:,i)./(BGposs*BGweight');

            BGtotal = sum(BGhmatrix(:,i));
            BGweight(i) = BGtotal/size(BGhmatrix,1);
            BGmean(i,:) = sum(repmat(BGhmatrix(:,i),1,featuredim).*TrainsampleDCT_BG)'/BGtotal;

            BGstdUP = zeros(featuredim,featuredim);

            for k = 1:BG_row
               difference = TrainsampleDCT_BG(k,:)-BGmean(i,:);
               BGstdUP = BGstdUP + BGhmatrix(k,i)*diag(diag(difference'*difference));
            end
            BGstd(:,:,i) = BGstdUP/BGtotal;
        end
    end
    FGtotalweight(:,:,classes) = FGweight(:,:);
    BGtotalweight(:,:,classes) = BGweight(:,:);
    
    FGtotalmean(:,:,classes) = FGmean(:,:);
    BGtotalmean(:,:,classes) = BGmean(:,:);
    
    FGtotalstd(:,:,:,classes) = FGstd(:,:,:);
    BGtotalstd(:,:,:,classes) = BGstd(:,:,:);
end
%Same thing for BG

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
for i=1:totalblock%zigzag
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

%-----------------------------------------------------------------------
imagefeature = ones(size(zigblocks,2),64);%ones(size(zigblocks,2),8);
for k = 1:size(zigblocks,2)
    imagefeature(k,:) = zigblocks{k};%(originbig); %- max(blocks{k})*ones(size(blocks{k}));
end
%-----------------------------------------------------------------------
%Image features extracted.

graphaftertotal = {}; 
counter = 1;
graphafter = ones(1,totalblock); 

errorstorer = zeros(mixture,feanum,mixture);
counter = 1;
for fea = 1:feanum
    for classes1 = 1:5%5mixers
        for classes2 = 1:5
            graphafter = ones(1,totalblock); 
            for j = 1:totalblock
               probability1 = 0;
               probability2 = 0;
              for k = 1:component
                  probability1 = probability1+mvnpdf(imagefeature(j,1:fea),FGtotalmean(k,1:fea,classes1),FGtotalstd(1:fea,1:fea,k,classes1))*FGtotalweight(1,k,classes1);
                  probability2 = probability2+mvnpdf(imagefeature(j,1:fea),BGtotalmean(k,1:fea,classes2),BGtotalstd(1:fea,1:fea,k,classes1))*BGtotalweight(1,k,classes2);
              end
              if(sum(probability1)>=sum(probability2))
                 graphafter(1,j) = 0;
              end
              
            end
            hohoho = mask(1:247,1:262);
            errorstorer(classes1,fea,classes2) = sum(abs(graphafter-reshape(hohoho,[1,247*262])))/(247*262);
        end
    end
end
hahaha1 = ones(5,11) - errorstorer;%the error rate
graphafter = reshape(graphafter,[yint,xint]);
%graphafter = rot90(graphafter,1);
pcolor(graphafter)