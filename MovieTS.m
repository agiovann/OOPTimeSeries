classdef MovieTS < timeseries
    % This class implements the general movie class extending the
    % timeseries object. It is used to manipulate movies created with
    % scanimage on different channels. Scanimage channels can also be used
    % to collect and synchronize a wide variety of data, like triggers,
    % electrophysiology, behavior, etc...
    % Movies can also be collected in chunks and then joined into a single
    % movie.
    properties (SetAccess = public)
        frameRate
        % if movies are collected in chunks and then joined this variable contain the name of each
        movieNames
        % number of joined movies
        numMovies
        % number of frames per movie if joined
        numFramesPerMovie
        % lines in each frame
        linesPerFrame
        % pixels in each line
        pixelsPerLine
        % if some frames are removed, this keeps track of them
        removedFrames
        % directory where the movie pbject will be stored.e
        dirFiles
    end
    
    methods (Static)
        function [movieMat,timeVector]=loadMovieChain(movies,channelId,numFramesSkipBeginning)
            % Join multiple scanimage  movies colected from the same region of interest. 
            %[movieMat,timeVector]=loadMovieChain(movies,channelId,numFramesSkipBeginning)
            % This functions join many movies colected from the same region
            % of interest.
            % movies: cells containing the names of the scanimage tiff movies to be
            % joined.
            % channelId: index of the channel to be extracted from the
            % movies. If 0 it will return empty vectors.
            % numFramesSkipBeginning: sometumes the first frames collected
            % are dark. This allows to exclude some frames from the
            % beginining of each movie.
            % movieMat: 3D array containing the movie
            % timeVector: time vector of the movie exctracted from the
            % metainformation contained in the tiff files
            parameters=loadMovie(movies,channelId,numFramesSkipBeginning);
            if(channelId~=0)
                timeVector=(0:(parameters.numFrames-1))*parameters.frameRate;
                movieMat=parameters.mov;
            else
                timeVector=[];
                movieMat=[];
            end
        end
        
        
        function [movieMat,timeVector,trialType,trigsSecs]=loadMoviesBehavior(movies,timeSkipBeginning,delayfromTriggerMaster8)
            % DEPRECATED. load movies stored in matlab format from the camera filming the
            % behavior and extracts information about the eyelid closure and locomotion. It removes corrupted frames
            % and interpolate missing frames
            % movies: cells containing the matlab movies, if empty a gui
            % will open for selecting the movies manually
            % timeSkipBeginning: time to skip at the beginning of each
            % movie if necessary
            % delayfromTriggerMaster8: represents the delay implemented in
            % master 8 between the trigger starting imaging acquisition and
            % the delivery of the CS onset.
            % movieMat: concatenatd and preprocessed  movies
            % timeVector: time stam p for movieMat
            %
            if isempty(delayfromTriggerMaster8)
                delayfromTriggerMaster8=1.936;
            end
            
            %% manually load movies
            if isempty(movies)
                filenamesP=uipickfiles;
                % open a gui to select the files
                movies={};
                for ll=1:numel(filenamesP)
                    [a,b,c]=fileparts(filenamesP{ll})
                    movies{ll}=['./' b c];
                end
            end
            %%
            load(movies{1},'mov','timeVect','lastTrial')
            % we are resampling so that all the movies are comparable
            ntrials=numel(movies);
            frate=round(200*mode(diff(timeVect)))/200;
            durationTrial=ceil(100*timeVect(end))/100;
            nFrPerMovie=ceil(durationTrial/frate);
            numFramesSkipBeginning=round(timeSkipBeginning/frate);
            movieMat=uint8(zeros(size(mov,1),size(mov,2),ntrials*(nFrPerMovie-numFramesSkipBeginning)));
            mov=double(mov); % necessary because of the ROI extraction process
            %% For removing bad frames
            meanMovLev=mean(mean(mean(mov,1),2));
            stdMovLev=std(mean(mean(mov,1),2));
            threshMov=meanMovLev-6*stdMovLev;
            newFramesPerMovie=nFrPerMovie-numFramesSkipBeginning; %
            newDurationTrial=durationTrial-timeSkipBeginning;
            trialType=[];
            for mmvv=1:ntrials
                if mod(mmvv,5) == 0
                    disp(mmvv) %counter
                end
                load(movies{mmvv},'mov','timeVect','lastTrial')
                newMov=zeros(size(mov,1),size(mov,2),nFrPerMovie)*NaN;
                timeVect(timeVect==0)=[];
                timeVect(find(mean(mean(mov,1),2)<threshMov))=[];
                newMov(:,:,round(timeVect/frate))=mov(:,:,round(timeVect/frate));
                %% remove corrupted frames
                badGuys=find(isnan(mean(mean(newMov,1),2)) | mean(mean(newMov,1),2)<threshMov);
                for bg=badGuys'
                    newMov(:,:,bg)=nanmax(newMov(:,:,min(max(1,[bg-1:bg+1]),nFrPerMovie)),[],3);
                    if isnan(mean(mean(newMov(:,:,bg),1),2)) || mean(mean(newMov(:,:,bg),1),2)<threshMov
                        newMov(:,:,bg)=nanmax(newMov(:,:,min(max(1,[bg-2:bg+2]),nFrPerMovie)),[],3);
                    end
                    if isnan(mean(mean(newMov(:,:,bg),1),2)) || mean(mean(newMov(:,:,bg),1),2)<threshMov
                        newMov(:,:,bg)=nanmax(newMov(:,:,min(max(1,[bg-3:bg+3]),nFrPerMovie)),[],3);
                    end
                    if isnan(mean(mean(newMov(:,:,bg),1),2)) || mean(mean(newMov(:,:,bg),1),2)<threshMov
                        error('nan frame!')
                    end
                end
                %%
                newMov(:,:,1:numFramesSkipBeginning)=[];
                trialType(mmvv)=lastTrial;
                movieMat(:,:,((mmvv-1)*newFramesPerMovie+1):mmvv*newFramesPerMovie)=newMov;
                plot(squeeze(mean(mean(movieMat,1),2)))
            end
            timeVector=frate:frate:(newDurationTrial*ntrials);
            trigsSecs=delayfromTriggerMaster8:newDurationTrial:(newDurationTrial*ntrials);
        end
        
        
        function [movieMat,timeVector]=loadMovieOnePhoton(movies,channelId,numFramesSkipBeginning)
            %% WIP probably should not be here. loads movies collected with one photon imaging
            parameters=loadMovie(movies,channelId,numFramesSkipBeginning);
            if(channelId~=0)
                timeVector=(0:(parameters.numFrames-1))*parameters.frameRate;
                movieMat=parameters.mov;
            else
                timeVector=[];
                movieMat=[];
            end
        end
        
        %%
        %         %%
        %         function obj=loadMovieObj(name)
        %             load(name,'obj');
        %         end
        %% 
        function infoImg=getImageInformation(movieName)
            % get metainformation stored by scanimage in eachmovie
            infoImg=agReadImageInformation(movieName);
        end
        
        
        
    end
    
    methods (Access = protected)
    end % methods
    methods (Access = public)
        %% constructor
        function obj = MovieTS(movMat,timeVect,movNames)
            % Constructor getting a matrix and the time vector. If multiple
            % movies are joined the names are stored
            %  obj = MovieTS(movMat,timeVect,movNames)
            % movMat: 3d matrix containing the movie
            % timeVect: time stamp for each frame
            % movNames: cell containing the names of the movies originally
            % joined to obtain movMat
            
            obj=obj@timeseries(movMat,timeVect);
            obj.movieNames=movNames;
            obj.numMovies=length(obj.movieNames);
            obj.numFramesPerMovie=obj.Length/obj.numMovies; %assumes same number of frames per movie
            obj.numMovies=length(obj.movieNames);
            obj.numFramesPerMovie=obj.Length/obj.numMovies;
            ss=obj.getdatasamplesize;
            obj.linesPerFrame=ss(1);
            obj.pixelsPerLine=ss(2);
            obj.frameRate=median(diff(timeVect)); % assumes equally spaced timestamp
            %% names the directory where to store the object based on the movNames
            if ~isempty(movNames)
                if numel(movNames)==1
                    obj.dirFiles=movNames{1}(1:end-4);
                else
                    obj.dirFiles=[movNames{1}(1:end-4) '_X' num2str(length(movNames))];
                end
                mkdir(obj.dirFiles)
            else
                obj.dirFiles='';
            end
        end
        %% remove frames around triggers
        function obj=movieWithoutPeriTrigger(obj,timeBefore,timeAfter,triggersSec)
            % obj=movieWithoutPeriTrigger(obj,timeBefore,timeAfter,triggersSec)
            
            for jj=1:length(triggersSec)
                trigTime=triggersSec(jj);
                obj=obj.delsample('Index',find(obj.Time>=(trigTime-timeBefore) & obj.Time<=(trigTime+timeAfter)));
            end
        end
        %%
        function play(obj,gain,frameRate)
            if(nargin<3)
                frameRate=obj.frameRate;
            end
            %play movie play(obj,gain)
            mv=obj.Data;
            if isa(obj.Data,'uint8')
                mv=round(mv*gain);
                implay(mv,1/frameRate);
            else
                implay(mv./max(mv(:))*gain,1/frameRate);
            end
            
        end
        
        function hh=plot(obj,varargin)
            hh=plot(mean(obj),varargin{:});
        end
        
        
        function scanSig=getScanSignal(obj)
            scanTmp=ipermute(obj.Data,[2 1 3]);
            scanSig=scanTmp(:);
            timeVectHR=(1:length(scanSig))*(obj.frameRate/obj.linesPerFrame/obj.pixelsPerLine);
            scanSig=mytimeseries(scanSig,timeVectHR);
        end
        
        function obj=filter(obj,b,a,dimension)
            obj.Data=filter(b,a,obj.Data,NaN(1,max(length(a),length(b))-1),dimension);
        end
        
        function obj=power(obj,expon)
            obj.Data=obj.Data.^expon;
        end
        
        function maxTS=max(obj)
            dat=obj.Data;
            dat=squeeze(max(max(dat,[],1),[],2));
            maxTS=mytimeseries(dat,obj.Time);
        end
        
        function obj=clipMovie(obj,maxValue)
            dat=obj.Data;
            dat=single(dat);
            %              maxVal=quantile(dat,maxValuePercentile,3);
            dat(dat>maxValue)=maxValue;
            obj.Data=dat;
        end
        
        %         function meanTS=nanmean(obj)
        %             dat=obj.Data;
        %             dat=squeeze(nanmean(nanmean(dat,1),2));
        %             meanTS=timeseries(dat,obj.Time);
        %         end
        function [newObj,newTrigs]=concatenateMoviesAroundTrigger(obj,trigsSecs,timeBefore,timeAfter)
            timeBefore=round(timeBefore/obj.frameRate)*obj.frameRate
            timeAfter=round(timeAfter/obj.frameRate)*obj.frameRate
            duration=timeBefore+timeAfter;
            movCorrected=[];
            numsamples=[];
            xshifts=[];
            yshifts=[];
            data=[];
            for kk=1:numel(trigsSecs)
                disp(kk)
                mm=obj.getFramesAlignedToTrigger(trigsSecs(kk),timeBefore,timeAfter-obj.frameRate);
                data=cat(3,data,mm.Data);
            end
            newTrigs=(timeBefore):(duration):(duration*numel(trigsSecs));
            timeVect=0:obj.frameRate:(obj.frameRate*(size(data,3)-1));
            idxTriggersMaintained=ceil(trigsSecs/obj.frameRate/obj.numFramesPerMovie);
            strToCall='(data,timeVect,obj.movieNames(idxTriggersMaintained));';
            strToCall=['newObj=' class(obj) strToCall];
            eval(strToCall)
        end
        
        function obj=diff(obj)
            dat=obj.Data;
            diffdat=diff(cat(3,dat(:,:,1),dat),1,3);
            obj.Data=diffdat;
        end
        
        
        function meanTS=mean(obj)
            dat=obj.Data;
            dat=squeeze(mean(mean(dat,1),2));
            meanTS=mytimeseries(dat,obj.Time);
        end
        
        function newts=f(obj,funct)
            % apply function to each element in the timeserie
            cellMat=mat2cell(obj.Data,obj.linesPerFrame,obj.pixelsPerLine,repmat(1,[1 obj.Length]));
            newvalues=squeeze(cellfun(funct,cellMat));
            newts=mytimeseries(newvalues,obj.Time);
        end
        
        function obj=filterROI(obj,mask)
            mask=double(mask);
            mask(mask==0)=NaN;
            obj.Data=double(obj.Data).*repmat(mask,[1 1 obj.Length]);
        end
        
        function showFrame(obj,frameIdx)
            imagesc(obj.getdatasamples(frameIdx))
        end
        
        function outp=zProject(obj,varargin)
            if(numel(varargin)<1)
                quantil=.5;
            else
                quantil=varargin{1};
            end
            if(numel(varargin)<2)
                if isa(obj.Data,'uint8')
                    maxLum=mean(obj.Data(:))+4*std(single(obj.Data(:)));
                else
                    maxLum=mean(obj.Data(:))+4*std(obj.Data(:));
                end
            else
                maxLum=varargin{2};
            end
            
            if isa(obj.Data,'uint8')
                img=quantile(obj.Data,quantil,3);
            else
                img=quantile(single(obj.Data),quantil,3);
            end
            if(nargout==0)
                imagesc(img,[0 maxLum]);
                colormap(gray)
            else
                outp=img;
            end
        end
        
        function saveFig(obj,figId,name,format)
            saveas(figId,['./' fullfile(obj.dirFiles,name)],format);
        end
        
        function saveToDir(obj,name)
            save(['./' fullfile(obj.dirFiles,name)],'obj');
        end
        
        function newobj=getFramesAlignedToTrigger(obj,trigsSec,secBefore,secAfter)
            if ~isempty(trigsSec)
                samplesBefore=round(secBefore/obj.frameRate);
                samplesAfter=round(secAfter/obj.frameRate);
                triggersSamples=round(trigsSec/obj.frameRate)+1;
                sampVectTrial=(-samplesBefore:samplesAfter);
                timeVectTrial=sampVectTrial*obj.frameRate;
                data=obj.Data;
                %                 colormap(gray)
                allimageAvg=zeros(size(data,1),size(data,2),length(sampVectTrial));
                for sm=1:length(sampVectTrial)
                    frIdx=min(max(1,triggersSamples+sampVectTrial(sm)),obj.Length);
                    %                     disp(frIdx')
                    imageAvg=nanmedian(data(:,:,frIdx),3);
                    allimageAvg(:,:,sm)=imageAvg;
                    %                     imagesc(imageAvg,[0 3])
                end
                clstr=eval(['@' class(obj)]);
                %                 newobj=clstr(allimageAvg,timeVectTrial,{[obj.dirFiles '_WVF.tif']});
                newobj=clstr(allimageAvg,timeVectTrial,{''});
                
            else
                newobj=[];
            end
            
        end
        
        function saveToAvi(obj,filename,frate,multFact)
            writerObj = VideoWriter(filename);
            writerObj.FrameRate = frate;
            open(writerObj);
            cm=colormap(gray(256));
            mov=obj.Data;
            mov=mov./max(mov(:));
            mov=uint8(mov*256-1).*multFact;
            for jj=1:obj.Length
                M(jj)=im2frame(mov(:,:,jj),cm);
            end
            writeVideo(writerObj,M);
            close(writerObj);
        end
        
        %          function saveToTiff(obj,filename,frate)
        %             writerObj = VideoWriter(filename);
        %             writerObj.FrameRate = frate;
        %             open(writerObj);
        %             cm=colormap(gray(256));
        %             mov=uint8(obj.Data*256-1);
        %             for jj=1:obj.Length
        %                 M(jj)=im2frame(mov(:,:,jj),cm);
        %             end
        %             writeVideo(writerObj,M);
        %             close(writerObj);
        %         end
        %
        
    end
end
%%

