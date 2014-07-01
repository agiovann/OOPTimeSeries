classdef MovieTSMagnet < MovieTS    
    properties (SetAccess = public)
        
    end
    
    methods (Static)
         
    end
    
    methods (Access = protected)
       
    end
    
    methods (Access = public)
        function obj = MovieTSMagnet(movMat,timeVect,movNames)
            obj=obj@MovieTS(movMat,timeVect,movNames);
        end
        
%         function [trigsSecs,trigSamp]=extractTriggers(obj,minDistTrigsSeconds)
%              minDistTrigs=ceil(minDistTrigsSeconds/obj.frameRate);
%              maxTS=max(obj);
%             [allTriggersAlignedUS,~,~,~]=agExtractTriggers(maxTS.Data,minDistTrigs);
%             trigsSecs=(allTriggersAlignedUS-1)*obj.frameRate;
%             trigSamp=allTriggersAlignedUS;
%         end
        
        function traceDS=getTrace(obj,decFact)
            traceDS= decimate(obj.getScanSignal,decFact);
        end
        
        function [wvfEyelid,sigEyelid]=getEyelidMatrix(obj,trigsSecs,tbef,taft,removeBL)
            if ~exist('removeBL')
                removeBL=1;
            end
            sigEyelid=obj.getTrace(128*4);               
            wvfEyelid=sigEyelid.extractWaveforms(trigsSecs,tbef,taft); 
%             wvfEyelid=mytscollection(wvfEyelid);
            if removeBL
                  wvfEyelid=wvfEyelid.removeBLWaveforms(1.25,-.25);
            end
        end
        
        function [amplitudes,wvfEyelid,sigEyelid]=getNormalizedMatrixCRs(obj,trigsSecs,triggers,trigsUS,tbef,taft)
            sigEyelid=obj.getTrace(128*4);               
            wvfEyelidUS=sigEyelid.extractWaveforms(trigsSecs(vertcat(triggers{trigsUS})),1,1);                        
            wvfEyelidUS=wvfEyelidUS.removeBLWaveforms(1);
            maxEyelid=median(wvfEyelidUS.max);
            sigEyelid=sigEyelid.detrend(1)/maxEyelid;
            
            
            for jj=1:numel(triggers)
                if ~isempty(triggers{jj})
                    wvf=sigEyelid.extractWaveforms(trigsSecs(triggers{jj}),1,1);
                    wvf=wvf.removeBLWaveforms(1);                    
                    amplitudes{jj}=max(wvf.getsampleusingtime(-tbef,taft));
                end    
            end
            wvfEyelid=sigEyelid.extractWaveforms(trigsSecs,1,1);            
        end
        
        function plot(obj)                  
            plot(getTrace(obj,128*10))
        end
    end
   
end
        