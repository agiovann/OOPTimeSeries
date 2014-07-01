classdef MovieTSImaging < MovieTS
    
    properties (SetAccess = public)
        
    end
    
    methods (Static)
        
    end
    
    methods (Access = protected)
        
    end
    
    methods (Access = public)
        %%
        function obj = MovieTSImaging(movMat,timeVect,movNames)
            obj=obj@MovieTS(movMat,timeVect,movNames);
        end
        %%
        function [trmat,timeVectMissing]=extractHRROIs(obj,dendmask,varargin)
            p = inputParser;
            defaultCutOff = .2;
            defaultInterpMethod = 'linear';
            defaultNumErasedLinesBeginning=0;
            defaultNumErasedLinesEnd=0;
            
            addRequired(p,'obj');
            addRequired(p,'dendmask',@isnumeric);
            addOptional(p,'method',defaultInterpMethod);
            addOptional(p,'NumErasedLinesBeginning',defaultNumErasedLinesBeginning);
            addOptional(p,'NumErasedLinesEnd',defaultNumErasedLinesEnd);

            addOptional(p,'CutOffLinesFewPixels',defaultCutOff);
            
            parse(p,obj,dendmask,varargin{:});
            interpMethod= p.Results.method;
            cutoffPixels=p.Results.CutOffLinesFewPixels;
            numErasedLines=p.Results.NumErasedLinesBeginning;    
            numErasedEnd=p.Results.NumErasedLinesEnd;                       

            
            numDends=max((size(dendmask,3)-1),1);
            frameRateHR=obj.frameRate/obj.linesPerFrame;
            timeVectHR=(1:numel(obj.Data(:))/obj.pixelsPerLine) * frameRateHR;
            trmat=mytscollection(timeVectHR);
            dendmask(end,:,:)=NaN;
            for kk=1:numDends                
                mask=dendmask(:,:,kk);
                % fix mask lines with few pixels
                numPixHorz=(nansum(mask,2));
                medianNumPixHorz=median(numPixHorz);
                badLines=find(numPixHorz>0 & numPixHorz<(cutoffPixels*medianNumPixHorz)); % have too few pixels to provide meaningful information
                mask(badLines,:)=0;
%                 mask(mask==0)=NaN;
                objFiltered=filterROI(obj,mask);                
                dat=objFiltered.Data;
                dat=nanmedian(dat,2);
                sigHR=dat(:);
                interval=numel(sigHR)/obj.numMovies;
                trBadGuys=agAddRectangleAroundIndexes(1:interval:numel(sigHR),numErasedEnd,numErasedLines,numel(sigHR));                
                sigHR(trBadGuys>0)=NaN;
                
                timeVectMissing=timeVectHR(isnan(sigHR));
                sigHRNotNaN=sigHR(~isnan(sigHR));
                timeVectHRNotNAN=timeVectHR(~isnan(sigHR));
                if ~isempty(timeVectHRNotNAN)
                    sigHRInterp=interp1(timeVectHRNotNAN,sigHRNotNaN,timeVectMissing,interpMethod);
                    sigHRInterp(isnan(sigHRInterp))=nanmean(sigHRInterp);
                    sigHR(isnan(sigHR))=sigHRInterp;
                end
                trmatHR=mytimeseries(sigHR,timeVectHR);
                trmat=trmat.addts(trmatHR,['D_' num2str(kk)]);
            end
        end
        %%
        function obj=normalizeMovie(obj,varargin)
            if(nargin>1)
                winSpanSec=varargin{1};
            else
                winSpanSec=10;
            end
            numFramesWindow=round(winSpanSec/obj.frameRate);
            [obj.Data]=agNormalizeMovie(obj.Data,numFramesWindow);
        end
        %%
        function obj=normalizeMovieQuantile(obj,quantileMin,quantileWinSec)
            if(nargin<2)
                quantileMin=.5;
            end
            if nargin<3
                quantileWin=round(12/obj.frameRate);
            else
                quantileWin=round(quantileWinSec/obj.frameRate);
            end
            convert=0;
            mov=obj.Data;
            if isa(mov,'uint8')
                convert=1;
                mov=single(mov);
            end
            if(quantileWinSec==0)
                movBL=quantile(mov,quantileMin,3);
                movBL=repmat(movBL,[1 1 size(mov,3)]);
                movNorm=(mov-movBL);
            else
                movBL=agRunningQuantile(mov,quantileMin,quantileWin);
                movNorm=[(mov-movBL)];
            end     
            if convert
               movNorm=mat2gray(movNorm,[0 double(max(movNorm(:)))]);
               movNorm=gray2ind(movNorm,255);
            end
            obj.Data=movNorm;
            
        end
        %%
        function newObj=getLineScanSignal(obj)
            data=ipermute(obj.Data,[2 1 3]);
            sig=reshape(data,obj.linesPerFrame,obj.pixelsPerLine*obj.numFramesPerMovie*obj.numMovies);
            sig=squeeze(nanmean(sig))';
            tim=(1:numel(sig))*obj.frameRate/obj.linesPerFrame;
            newObj=mytimeseries(sig,tim);
        end

        %%
        function newObj=getLineScanSignalOGB(obj)
            data=ipermute(obj.Data,[2 1 3]);            
            sig=reshape(data(:),obj.pixelsPerLine,obj.linesPerFrame*obj.numFramesPerMovie*obj.numMovies);
            sig=squeeze(nanmean(sig))';
            tim=(1:numel(sig))*obj.frameRate/obj.linesPerFrame;
            newObj=mytimeseries(sig,tim);
        end
       
        %%
        function trmat=extractROIs(obj,dendmask)
            numDends=max((size(dendmask,3)-1),1);
            trmat=mytscollection(obj.Time);
            for kk=1:numDends
                mask=dendmask(:,:,kk);
                maskmov=filterROI(obj,mask);
                %                 tr=nanmean(maskmov);
                tr=squeeze(nansum(nansum(maskmov.Data,1),2))/nansum(nansum(mask,1),2);
                tr=mytimeseries(tr,obj.Time);
                trmat=trmat.addts(tr,['D_' num2str(kk)]);
            end
        end
        %%
        function [movCorrected,numsamples,xshifts,yshifts,newTrigs]=motionCorrectTrialWise(movToCorrect,correctionChannel,trigsSecs,timeBefore,timeAfter,refframe,maxshiftobj)            
            duration=timeBefore+timeAfter;
            movCorrected=[];
            numsamples=[];
            xshifts=[];
            yshifts=[];            
            data=[];
            for kk=1:numel(trigsSecs)
                disp(kk)
                mm=movToCorrect.getFramesAlignedToTrigger(trigsSecs(kk),timeBefore,timeAfter-movToCorrect.frameRate);
                mmCh=correctionChannel.getFramesAlignedToTrigger(trigsSecs(kk),timeBefore,timeAfter-movToCorrect.frameRate);
                %     mm=mm-mm.zProject(.02);
                [mm,ns,xs,ys]=mm.motionCorrect(mmCh,refframe,maxshiftobj);
                 numsamples=[numsamples ns];
                 xshifts=[xshifts xs];
                 yshifts=[yshifts ys];
                data=cat(3,data,mm.Data);                
                close all
            end
            newTrigs=timeBefore:(duration):(duration*numel(trigsSecs));
            timeVect=0:movToCorrect.frameRate:(movToCorrect.frameRate*(size(data,3)-1));
            idxTriggersMaintained=ceil(trigsSecs/movToCorrect.frameRate/movToCorrect.numFramesPerMovie);
            movCorrected=MovieTSImaging(data,timeVect,movToCorrect.movieNames(idxTriggersMaintained));
%             movCorrected.dirFiles=movToCorrect.dirFiles;
        end

%%
        function [movToCorrect,numsamples,xshifts,yshifts]=motionCorrect(movToCorrect,correctionChannel,refframe,maxshiftobj)
            [movMotCorrTmp,xshifts,yshifts,numsamples,maxxshift,maxyshift,minxshift,minyshift]=agmotcorshift01({movToCorrect.Data,correctionChannel.Data},1, refframe, maxshiftobj,1);
            movToCorrect.Data=normrnd(mean(movToCorrect.Data(:)),std(movToCorrect.Data(:))/3,size(movToCorrect.Data));
            % set the same dimensional matrix
            ny=(1+maxyshift:movToCorrect.linesPerFrame+minyshift);
            nx=(1+maxxshift:movToCorrect.pixelsPerLine+minxshift);
            movToCorrect.Data(ny,nx,:)=movMotCorrTmp;
        end
        %%
        function [obj,xshifts,yshifts,numsamples]=motionCorrectGauss3D(obj,refframe,maxshiftobj,x,y,z)
            MAXSHIFT=6;
            if(nargin<3)
                maxshiftobj=MAXSHIFT;
                x=2;
                y=2;
                z=round(3/obj.frameRate);
            end
            warning('EmployingGauss 3D filtered Movie for Motion Correction')
            %                 movRunMean=gaussianSmooth3D(obj.Data,x,y,z);
            movRunMean=imfilter(obj.Data,fspecial3('gaussian',[x y z]));
            [movMotCorrTmp,xshifts,yshifts,numsamples,maxxshift,maxyshift,minxshift,minyshift]=agmotcorshift01({obj.Data,movRunMean},1, refframe, maxshiftobj,1);
            
            obj.Data=normrnd(mean(obj.Data(:)),std(obj.Data(:))/3,size(obj.Data));
            % set the same dimensional matrix
            ny=(1+maxyshift:obj.linesPerFrame+minyshift);
            nx=(1+maxxshift:obj.pixelsPerLine+minxshift);
            obj.Data(ny,nx,:)=movMotCorrTmp;
        end
        %%
        function obj=smooth3(obj,type,dims,sd)
            % example (obj.smooth3('gaussian',[1 1 3],1)
            x=dims(1);
            y=dims(2);
            z=dims(3);
            if mod(z,2)==0
                z=z+1;
            end
            obj.Data=smooth3(obj.Data,type,[x y z],sd);
        end
        %%
        function dendmaskOrig=extractGCMaskNoPCAICA_2(obj,trigs,minPixMask,threshSTD)
            newobj=getFramesAlignedToTrigger(obj,trigs,2,2);
            signalAvg=(squeeze(median(median(newobj.Data,1))));
            
            bl=quantile(signalAvg,.08);
            max=quantile(signalAvg,.92);
            sigSize=(max-bl);
            thr=bl+sigSize*.25;
            idxGood=find(signalAvg>thr);
            
%             imagesc(mean(newobj.smooth3('gaussian',[3 3 3],.65).Data(:,:,idxGood),3))
            movSTD=(std(newobj.smooth3('gaussian',[3 3 1],.65).Data(:,:,idxGood),[],3));
            movSTDThresh=(movSTD>threshSTD);                      
            maskGC=[];
            I=double(movSTDThresh);
            BW = im2bw(I, graythresh(I));
            cc = bwconncomp(BW, 4);
            pixels=cc.PixelIdxList;
            for tt=1:cc.NumObjects
                px=pixels{tt};
                if numel(px)>minPixMask
                    m=zeros(size(BW));
                    m(px)=1;
                    maskGC=cat(3,maskGC,m);
                end
            end
            dendmaskOrig=cat(3,maskGC,sum(maskGC,3));    
            
            colormap(gray)
            subplot(2,2,1)
            imagesc(movSTD)
            axis image
            subplot(2,2,2)
            imagesc(movSTDThresh)
            axis image
            subplot(2,2,3)
            imagesc(sum(maskGC,3))
            %
            axis image
            
        end
        %%
        %%
        function dendmaskOrig=extractGCMaskNoPCAICA_3(obj,trigs,factResize,thresh)
            newobj=getFramesAlignedToTrigger(obj,trigs,2,2);
%             signalAvg=(squeeze(median(median(newobj.Data,1))));
            
            
            factResize=1; % for 20X 3X
            sizeRectout=4*factResize;
            sizeRectIn=2*factResize;
            sizeKernel=6*factResize;
            
            s=strel('rectangle',[sizeRectout sizeRectout]);                        
            kern=s.getnhood; 
            largeMat=blkdiag(zeros((sizeKernel-sizeRectout)/2),kern,zeros((sizeKernel-sizeRectout)/2));
            s=strel('rectangle',[sizeRectIn sizeRectIn]);                        
            kern2=s.getnhood;
            largeMat2=blkdiag(zeros((sizeKernel-sizeRectIn)/2),kern2,zeros((sizeKernel-sizeRectIn)/2));
            kernMat=largeMat-largeMat2;
            imshow(kernMat)
            
            typeAvg='mean';
            filtSize=1;
            isAdapthist=0;
            isimadjust=0;
            resizeFact=4;
            I=getProcessedImages(newobj,typeAvg,filtSize,isAdapthist,isimadjust,resizeFact);
            I2=I;   
            I2=adapthisteq(I2,'NumTiles',[newobj.linesPerFrame/8 newobj.pixelsPerLine/8],'clipLimit',0.007,'Distribution','uniform');

            se = strel('disk',3);
            I2=imsubtract(imadd(I2,imtophat(I2,se)), imbothat(I2,se));


            imshow(imresize(I2/2,3))
            thrDeconv=3;
            imshow(imresize(imfilter(imregionalmax(deconvblind(I2,double(kernMat))>3,26),double(sse)),3))
            deconvolvedImage=deconvblind(I2,double(kernMat));
            deconvPeaks=imregionalmax(deconvolvedImage,sizeKernel);
%             deconvPeaksThr=deconvPeaks>thrDeconv;
            deconvPeaks=deconvolvedImage.*(deconvPeaks>0);
            reconvImage=imfilter(deconvPeaks>2.5,kernMat,'same');
            imshow(reconvImage)
            figure
            imshow(I2/2)
            
            bl=quantile(I2(:),.08);
            max=quantile(I2(:),.92);
            sigSize=(max-bl);
            thr=bl+sigSize*.25;
            idxGood=find(I2>thr);
            
%             imagesc(mean(newobj.smooth3('gaussian',[3 3 3],.65).Data(:,:,idxGood),3))
%             movSTD=(std(newobj.smooth3('gaussian',[3 3 1],.65).Data(:,:,idxGood),[],3));
%             movSTDThresh=(movSTD>threshSTD);                      
            maskGC=[];
            I=double(I2>thresh);
            BW = im2bw(I, graythresh(I));
            cc = bwconncomp(BW, 4);
            pixels=cc.PixelIdxList;
            for tt=1:cc.NumObjects
                px=pixels{tt};
                if numel(px)>minPixMask
                    m=zeros(size(BW));
                    m(px)=1;
                    maskGC=cat(3,maskGC,m);
                end
            end
            dendmaskOrig=cat(3,maskGC,sum(maskGC,3));    
            
            colormap(gray)
            subplot(2,2,1)
            imshow(imresize(I2,3))
            axis image
            subplot(2,2,2)
            imshow(imresize(BW,3))
            axis image
            subplot(2,2,3)
            imagesc(dendmaskOrig(:,:,end))
            %
            axis image
            
        end
        %%
        function dendmaskOrig=extractGCMaskNoPCAICA(obj,threshPixels,minPixMask,threshQuant,gaussFiltWidth,gaussFiltSTD)
            %             threshPixels=.45;
            %             minPixMask=30;
            %             gaussFiltWidth=[5 5];
            %             gaussFiltSTD=.65;
            az=obj.zProject(threshQuant);
            importantPixels=(az>threshPixels).*mean(az,3);
            importantPixelsSM=filter2(fspecial('Gaussian',gaussFiltWidth,gaussFiltSTD),importantPixels);
            colormap(gray)
            subplot(2,2,1)
            imagesc(importantPixels)
            axis image
            subplot(2,2,2)
            imagesc(importantPixelsSM)
            axis image
            subplot(2,2,3)
            imagesc(importantPixelsSM>threshPixels)
            %
            axis image
            maskGC=[];
            I=double(az>threshPixels);
            BW = im2bw(I, graythresh(I));
            cc = bwconncomp(BW, 4);
            pixels=cc.PixelIdxList;
            for tt=1:cc.NumObjects
                px=pixels{tt};
                if numel(px)>minPixMask
                    m=zeros(size(BW));
                    m(px)=1;
                    maskGC=cat(3,maskGC,m);
                end
            end
            dendmaskOrig=cat(3,maskGC,sum(maskGC,3));
        end
        %%
        function dendmaskOrig=extractGCMaskHighSkewness(obj,threshSkewness,minPixMask,filterSize,gaussFiltSTD)
            %            threshSkewness=5;
            %             filterSize=[3 3 3];
            %             gaussFiltSTD=.65;
            movGauss=obj.smooth3('gaussian',filterSize,gaussFiltSTD);
            
            skewImage=k(movGauss.Data,[],3);
            colormap(gray)
            subplot(2,2,1)
            imagesc(movGauss.zProject)
            axis image
            subplot(2,2,2)
            imagesc(skewImage)
            axis image
            subplot(2,2,3)
            imagesc(skewImage>threshSkewness)
            %
            axis image
            maskGC=[];
            I=double(skewImage>threshSkewness);
            BW = im2bw(I, graythresh(I));
            cc = bwconncomp(BW, 4);
            pixels=cc.PixelIdxList;
            for tt=1:cc.NumObjects
                px=pixels{tt};
                if numel(px)>minPixMask
                    m=zeros(size(BW));
                    m(px)=1;
                    maskGC=cat(3,maskGC,m);
                end
            end
            dendmaskOrig=cat(3,maskGC,sum(maskGC,3));
        end
        %%
        function [dendmaskOrig,icamask]=extractPCDendritesPCAICA(obj,numICAComp,numSVDComp,minPixMask,minFreqAccepted)                        
            frameRate=obj.frameRate;        
            mov=obj.Data;
            [icasig, A, W, E, D,spcomps,x]=do_ica(mov,'tanh',numICAComp,numSVDComp);
            icamask=[];
            icamask.spcomps=spcomps;
            icamask.tpcomps=W;
            [dendmaskOrig,trmatOrig,coeffOfDifferenceNormTrMatOrig]=agautomDendExtraction02(mov,...
                frameRate,icamask.spcomps,...
                'minRatioDim',2,'minPixMask',minPixMask,'maxNlogLAccepted',10000,'minFreqAccepted',minFreqAccepted,'minSNRFarz',0,...
                'useWeightedMask',0,'isInteractive',1,'maxCorrAllowed',0.0);
        end
        %%
        function [dendmaskOrig,dendmaskF0]=extractPCDendritesManually(obj,I)
            if isempty(I)
               I=(obj.getdatasamples(1)); 
            end
            figure
            imshow(imresize(I,4))
            axis image
            dd=[];
            
            % extract manually ROIs
            cont=true;
            while cont
                a=roipoly;
                if isempty(a) || (sum(a(:))<3)
                    cont=false;
                else
                    dd=cat(3,dd,imresize(a,1/4));
                end
            end
            dendmaskOrig=cat(3,dd,sum(dd,3));
            
            % extract baseline dendritic mask
            figure
            dd=[];
            imshow(imresize(I*1,4))
            a=roipoly;
            % a=impoly;
            dd=cat(3,dd,imresize(a,1/4));
            dendmaskF0=cat(3,dd,sum(dd,3));                        
        end
        %%
        function I=getProcessedImages(obj,typeAvg,filtSize,isAdapthist,isimadjust,resizeFact)
            mov=obj.Data;            
            figure
            colormap(gray)
            switch typeAvg
                case 'mean'
                    img=(mean(mov,3));
                case 'median'
                    img=(median(mov,3));      
                case 'max'
                    img=(max(mov,[],3));
                case 'iqr'              
                    img=iqr(mov,3);
                case 'zProject'
                    img=obj.zProject;
            end
            I=mat2gray(img,[quantile(img(:),.0001) quantile(img(:),.9999)]);
            I = medfilt2(I,[filtSize filtSize]);
            % SE = strel('disk',[2]);
            % I = imerode(I,SE);
            % I = imadjust(I,[],[],1.5);
            % I=adapthisteq(I,'Range','original','Distribution','exponential','numtiles',[6 10],'ClipLimit',.0000001);
            if isAdapthist
                I=adapthisteq(I);
            end
            if isimadjust
                  I=imadjust(I,stretchlim(I),[0 1]);
            end                        
            imshow(imresize(I,resizeFact))
            axis image
        end
    end
    
end

