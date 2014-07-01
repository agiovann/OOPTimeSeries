classdef MovieTSImagingMultiPlane 
    
    properties (SetAccess = public)
        moviesColl;
        frameRate;
        numMovies;
        numFramesPerMovie;
        linesPerFrame;
        pixelsPerLine;        
        dirFiles;
    end
    
    methods (Static)
        
    end
    
    methods (Access = protected)
        
    end
    
    methods (Access = public)
        %%
        function obj = MovieTSImagingMultiPlane(objMovie,numFramesPerMovie,movName)
            movTot=objMovie.Data;
            moviesColl={};
            counter=0;
            for kk=1:numFramesPerMovie:objMovie.Length
                counter=counter+1;
                timeVect=(1:numFramesPerMovie)*objMovie.frameRate;
                mov=objMovie.Data(:,:,kk:(kk+numFramesPerMovie-1));
                mov=mov./(quantile(mov(:),.8));
                moviesColl{counter}=MovieTSImaging(mov,timeVect,[]);
            end
            obj.moviesColl=moviesColl;
            obj.frameRate=objMovie.frameRate;
            obj.numFramesPerMovie=numFramesPerMovie;
            obj.numMovies=length(moviesColl);
            obj.pixelsPerLine=objMovie.pixelsPerLine;
            obj.linesPerFrame=objMovie.linesPerFrame;
            obj.dirFiles=movName;
        end
        
        function [obj]=motionCorrectWSmooth(obj,refframe,maxshiftobj,xyz,sigmaSm)                                     
            moviesColl={};
            for mm=1:obj.numMovies
                movClOrig=obj.moviesColl{mm};
                movSmooth=movClOrig.smooth3('gaussian',[xyz(1) xyz(2) xyz(3)],sigmaSm);
                [moviesColl{mm},numsamples]=movClOrig.motionCorrect(movSmooth,refframe,maxshiftobj);            
            end
           obj.moviesColl=moviesColl;
        end
        
        function newObj=getWholeMovie(obj)
            alldata=[];
            for mm=1:obj.numMovies
                alldata=cat(3,alldata,obj.moviesColl{mm}.Data);
            end
            timev=(1:size(alldata,3))*obj.frameRate;
            newObj=MovieTSImaging(alldata,timev,[]);
            newObj.dirFiles=obj.dirFiles;
        end
        
        function [obj,trSec]=getFramesAlignedToTrigger(obj,trigsSec,secBefore,secAfter)
            moviesColl={}; 
            triggersTrace=zeros(1,obj.numFramesPerMovie*obj.numMovies);
            triggersTrace(max(1,ceil(trigsSec/obj.frameRate)))=1;
            triggersTrace=reshape(triggersTrace, obj.numFramesPerMovie,obj.numMovies);
            trSec={};
            for mm=1:obj.numMovies   
                trSec{mm}=find(triggersTrace(:,mm))'*obj.frameRate;
                moviesColl{mm}=getFramesAlignedToTrigger(obj.moviesColl{mm},trSec{mm},secBefore,secAfter);
%                 plot(obj.moviesColl{mm})
%                 hold on
%                 plot(obj.moviesColl{mm}.Time,triggersTrace(:,mm),'r')
%                 uiwait
            end
            obj.moviesColl=moviesColl;
        end
        
        function newObj=zProject(obj,varargin)
            mat=[] 
            for mm=1:obj.numMovies                
                mat=cat(3,mat,zProject(obj.moviesColl{mm},varargin{:}));
            end
            newObj=MovieTSImaging(mat,(1:obj.numMovies)*obj.numFramesPerMovie*obj.frameRate,'');
        end
        function play(obj,gain,frate)
            play(getWholeMovie(obj),gain,frate)
        end
        
        function obj=smooth3(obj,type,dims,sd)
            moviesColl={}; 
            for mm=1:obj.numMovies                
                moviesColl{mm}=smooth3(obj.moviesColl{mm},type,dims,sd);
            end
            obj.moviesColl=moviesColl;
        end
        
         function saveToDir(obj,name)
            mkdir(['./' obj.dirFiles])
            save(['./' fullfile(obj.dirFiles,name)],'obj');
        end
    end
    
end

