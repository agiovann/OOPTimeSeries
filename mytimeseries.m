classdef mytimeseries < timeseries
    % The following properties can be set only by class methods
    properties (SetAccess = public)
        frameRate
        qualityInfo;
    end
    methods (Static)
    end
    
    methods (Access = protected)
        
    end % methods
    methods (Access = public)
        function obj = mytimeseries(varargin)
            % Constructor getting a cell ararray of movie names, join them
            % togethrt and set the variables inside the class
            obj=obj@timeseries(varargin{:});
            obj.frameRate=median(diff(obj.Time));
            obj.qualityInfo=obj.QualityInfo;
            
        end
        
        function obj=filter(obj,b,a)
            initVal=nanmedian(obj.Data(:));
            obj.Data=filter(b,a,obj.Data,initVal*ones(1,max(length(a),length(b))-1));
        end
        
        function newobj=decimate(obj,decimationFactor)
            intrpData=decimate(obj.Data,decimationFactor);
            newobj=mytimeseries(intrpData,downsample(obj.Time,decimationFactor));
        end
        
        function obj=nanmeanAG(obj,dimension)
            dat=obj.Data;
            dat=nanmean(dat,dimension);
            obj.Data=dat;
        end
        function obj=nanmedianAG(obj,dimension)
            dat=obj.Data;
            dat=nanmedian(dat,dimension);
            obj.Data=dat;
        end
        function obj=removeSignalWRToTriggers(obj, timeBef,timeAft,trigsSecs)
            dat=obj.Data;
            trigsSamp=round(trigsSecs/obj.frameRate);
            frBef=round(timeBef/obj.frameRate);
            frAft=round(timeAft/obj.frameRate);
            windowUS1=agAddRectangleAroundIndexes(trigsSamp,frBef,frAft,obj.Length);
            dat=dat.*windowUS1';
            obj.Data=dat;
        end
        
        function newobj=getSpikeShape(obj,eventsC,timeBef,timeAft)
%            trigsSecs=find(eventsC.Data)*obj.frameRate;
           trigsSecs=eventsC.Time(find(eventsC.Data>0));
           newobj=obj.extractWaveforms(trigsSecs,timeBef,timeAft); 
        end
        
        
        function obj=getDFFQuantile(obj,quantilemin,quantilewinsec,minValClip)
            if nargin<2
                quantilemin=.08;
            end
            if nargin < 3
                quantilewin=round(12/obj.frameRate);
            else
                quantilewin=round(quantilewinsec/obj.frameRate);
            end
            sig=obj.Data;
            if(numel(find(sig<minValClip))>0)
               warning(['***** Changed ' num2str(numel(find(sig<minValClip))) ' values to more then minimum value of 10']) 
               sig(sig<minValClip)=minValClip;
            end
            if quantilewin==0
                bl=quantile(sig,quantilemin);
            else                
                if(numel(quantilewin)>1)
                    quantilewin=quantilewin(1);
                    quantileWinDouble=quantilewin*2;
                    sigShifted=circshift(sig,quantileWinDouble);
                    sigShifted(1:quantileWinDouble)=quantile(sig(1:quantileWinDouble),quantilemin)';
                    bl=agRunningQuantile(sigShifted,quantilemin,quantileWinDouble)';
                    bl=circshift(bl,-quantileWinDouble);
                else
                    error('Only One Argument Old version!')
                    bl=agRunningQuantile(sig,quantilemin,quantilewin)';
                end
%                 
            end
            
            sig=(sig-bl);
            %             sig=sig-nanmedian(sig);
%             sig=sig./abs(bl);
            sig=sig./smooth(bl,round(5./obj.frameRate));

            
            obj.Data=sig;
        end
        
        function obj=removeBLQuantile(obj,quantilemin,quantilewinsec)
            if nargin<2
                quantilemin=.08;
            end
            if nargin < 3
                quantilewin=round(1/obj.frameRate);
            else
                quantilewin=round(quantilewinsec/obj.frameRate);
            end
            sig=obj.Data;
            bl=agRunningQuantile(sig,quantilemin,quantilewin)';
            sig=(sig-bl);
            sig=sig-nanmedian(sig);
            obj.Data=sig;
        end
        
        function [wvfTS]=extractWaveforms(obj,trigsSecs,secBefore,secAfter)
            if ~isempty(trigsSecs)
                samplesBefore=round(secBefore/obj.frameRate);
                samplesAfter=round(secAfter/obj.frameRate);
                triggersSamples=floor(trigsSecs/obj.frameRate)+1;
                timeVectTrial=(-samplesBefore:samplesAfter)*obj.frameRate;
                wvf=extractWaveforms(squeeze(obj.Data)',triggersSamples,samplesBefore,samplesAfter,1);
                wvfTS=mytimeseries(timeseries(squeeze(wvf'),timeVectTrial));
            else
                wvfTS=mytimeseries([],[]);
            end
        end
        function saveToDir(obj,name)
            save(['./' name],'obj');
        end
        
        function obj=decimateNoFilter(obj,decFactor)
                tim=obj.Time(1:decFactor:end);
                dat=obj.Data;
                dat=dat(1:decFactor:end,1);
                obj=mytimeseries(dat,tim);
        end
        
        function newObj=decimateSpikes(obj,decFactor)
            spikesId=find(obj.Data);
            newSpikesId=round(spikesId/decFactor);
            newTime=obj.Time(1:decFactor:end);
            newSpikes=zeros(numel(newTime),1);
            newSpikes(max(1,newSpikesId))=1;
            newObj=mytimeseries(newSpikes,newTime);
        end
        
        
        function newobj=downsample(obj,decimationFactor)
            newdata=downsample(obj.Data,decimationFactor);
            newobj=mytimeseries(newdata,downsample(obj.Time,decimationFactor));
        end
        
        function newobj=removeBLWaveforms(obj,timeWinBefore,timeInitTrigger)
            if ~exist('timeInitTrigger')
                timeInitTrigger=0;
            end            
            sig=obj.Data;
            tim=obj.Time;
            newsig=bsxfun(@minus,sig',nanmedian(sig(tim<timeInitTrigger & tim>  (timeInitTrigger-timeWinBefore),:))')';
%             newsig=bsxfun(@minus,sig',quantile(sig(tim<timeInitTrigger & tim>-timeWinBefore,:),.08)')';
            newobj=mytimeseries(newsig,obj.Time);            
        end
        
        function newobj=getDFFQuantileWVF(obj,timeWinBefore,timeTrig)
            sig=obj.Data;
            tim=obj.Time;
            bl=nanmedian(sig(tim<timeTrig & tim>(timeTrig-timeWinBefore),:));
            newsig=bsxfun(@minus,sig',bl')';
            newsig=bsxfun(@times,newsig',1./bl')';
            newobj=mytimeseries(newsig,obj.Time);
        end
        
        function newobj=appfun(obj,fhandle,parameters)
            newdata=feval(fhandle,obj.Data,parameters{:});
            newobj=mytimeseries(newdata,resample(obj.Time,numel(obj,Data),numel(newData)));
        end
        
        function newobj=getSubCollection(obj,indexes)
            newdata=obj.Data;
            newdata=newdata(:,indexes);
            newobj=mytimeseries(newdata,obj.Time);
        end
        
        function newColl=getSubCollectionGroups(obj,allTriggersStr,triggersNames)
           
            if exist('triggersNames')
                triggers=allTriggersStr;
            else
                 triggers={};
                 triggersNames={};
                for kk=1:numel(allTriggersStr)
                    triggers{kk}=allTriggersStr(kk).indexTR;
                    triggersNames{kk}=allTriggersStr(kk).trialType;
                end
            end
            newColl=mytscollection(obj.Time);
            for pp=1:numel(triggers)
                if ~isempty(triggers{pp})
                    ts=obj.getSubCollection(triggers{pp});
                    newColl=newColl.addts(ts,triggersNames{pp});
                end
            end
        end

      
        
        
        function [kinetics]=computeRaiseAndDecayTime(obj,minPerc,maxPerc,minDFF,stimLength)
            % obj must be aligned to trigger so that time 0 is when the
            % triggers begins
            kinetics=struct;
            kinetics.riseTime=NaN;
            kinetics.decayHalfTime=NaN;
            kinetics.timeToBL=NaN;
            kinetics.timePeak=NaN;
            kinetics.timPeak90=NaN;
            kinetics.valPeak90=NaN;
            kinetics.valPeakMax=NaN;
            kinetics.startDecayTim=NaN;
            kinetics.timeStart=NaN;
            kinetics.startDecayVal=NaN;
            kinetics.Area=0;
            kinetics.TimeHalfArea=NaN;
            wvf=obj;
            tim=wvf.Time';
            avgTrace=wvf.Data';
            [val,peak]=max(avgTrace.*(tim>0));
            if ~(val<minDFF || isnan(val))
                peak90=find(avgTrace>=(maxPerc*val));
                peak90=peak90(1);
                %
                timPeak=tim(peak);
                timPeak90=tim(peak90);
                
                kinetics.timePeak=timPeak;
                kinetics.timPeak90=timPeak90;
                kinetics.valPeak90=avgTrace(peak90);
                kinetics.valPeakMax=val;
                % riseTime omputation
                start=peak90;
                while start>1 & (avgTrace(start)>(minPerc*val))
                    start=start-1;
                end
                timeStart=tim(start);
                kinetics.timeStart=timeStart;
                kinetics.riseTime=timPeak90-timeStart;
                
                % decay val computation
                [startDecayVal,startDecayIdx]=max(avgTrace.*(tim>stimLength));
%                 plot((avgTrace.*(tim>stimLength)))
                
                startDecayTim=tim(startDecayIdx(1));
                kinetics.startDecayVal=startDecayVal;
                idxHalfDecay=find(avgTrace<(.5*startDecayVal) & tim>startDecayTim);
                if(isempty(idxHalfDecay))
                    idxHalfDecay=NaN;
                else
                    idxHalfDecay=idxHalfDecay(1);
                end
                if ~(startDecayVal<minDFF || isnan(idxHalfDecay))
                    kinetics.startDecayTim=startDecayTim;
                    
                    timeHalfDecay=tim(idxHalfDecay);
                    kinetics.decayHalfTime=timeHalfDecay-startDecayTim;
                    % back to baseline
                    idxBL=find(avgTrace<(.1*val) & tim>timeHalfDecay);
                    if ~isempty(idxBL)
                        idxBL=idxBL(1);
                    else
                        idxBL=numel(avgTrace);
                    end 
                    if ~isnan(idxBL)
                        timBL=tim(idxBL);
                        kinetics.timeToBL=timBL-startDecayTim;
                        traceToIntegrate=avgTrace(tim<timBL & tim>timeStart);
                        kinetics.Area=sum(traceToIntegrate);
                        halfarea=kinetics.Area/2;
                        idxtmp=find(cumsum(traceToIntegrate)>halfarea);
                        idxHalfArea=idxtmp(1);
                        frate=median(diff(tim));
                        kinetics.TimeHalfArea=timeStart+idxHalfArea*frate;

                    end  
                end
            end
        end
        
        function obj=nonlinearityRemove(obj,nonlinpoly,nonlininputfact)
            newdata=obj.Data;    
             if isempty(nonlinpoly)
                load Nonlinearpolinomial;
            end
            if isempty(nonlininputfact)
                nonlininputfact=nonlininputfact;
            end
                                    
            if isvector(newdata)
                newdata=ppval(nonlinpoly,newdata);                                
            else
                for tr=1:size(newdata,2)
                   newdata(:,tr)=ppval(nonlinpoly,newdata(:,tr));
                end
            end
            newdata=newdata./nonlininputfact;
            obj.Data=newdata;   
        end
        
        function obj=nonlinearityRemoveSigmoid(obj,indicator,basalFluorescence)
            % basalFluorescence = fluoresce of the cell at basal calcium
            % level
            % indicator = type of indicator used, for now only GCaMP6f
            dat=obj.Data;    
            switch indicator
                case 'GCaMP6f'
                    kd=290e-9;
                    Rf=29.2;
                    nh=6;
                otherwise
                    error('Specify the indicator parameters')
            end
            %plot ca curve and compare the result of inveting the curve
            %itself
            plotCaCurve=0;
            if(plotCaCurve)
                CalciumConc=10.^[-8:.01:-4.5];
                plot(log10(CalciumConc),agsigmoid([Rf log10(kd) nh],log10(CalciumConc)))
                CalciumConc-10.^(aginvertSigmoid([Rf log10(kd) nh],agsigmoid([Rf log10(kd) nh],log10(CalciumConc))))
            end
            blCaLevelFluo=basalFluorescence;                        
            dat=dat-mode(round(50*(dat./max(dat)))/50)*max(dat);
            dat=dat+blCaLevelFluo;
            dat(dat<0)=NaN;                       
            newtrace=(10.^aginvertSigmoid([Rf log10(kd) nh],dat));
            newtrace=newtrace-mode(round(50*(newtrace./max(newtrace)))/50)*max(newtrace);
            newtrace(newtrace<0)=NaN;
            newtrace(isnan(newtrace))=normrnd(0,iqr(newtrace)/30,size(find(isnan(newtrace))));
            obj.Data=newtrace;   
        end
        
        function newobj=deconv(obj,tauFilt,durFilt)
            newdata=obj.Data;            
           
%             tauFilt=.3;
%             durFilt=1;
            tfilt=0:obj.frameRate:durFilt;
            flt=exp(-(tfilt/tauFilt));
            diffTerms=length(flt)-1;            
            if isvector(newdata)
                newdata=deconv(newdata,flt); 
                newdata=[newdata; repmat(median(newdata),[diffTerms 1])];
            else
                for tr=1:size(newdata,2)                   
                    newdata(:,tr)=deconv(newdata(:,tr),flt);
                    newdata(:,tr)=[newdata(:,tr); repmat(median(newdata(:,tr)),[diffTerms 1])];
                end
            end            
            newobj=mytimeseries(newdata,obj.Time);
        end
        
        
        function obj=getTraceWithoutPeriTrigger(obj,timeBefore,timeAfter,triggersSec)
            for jj=1:length(triggersSec)
                trigTime=triggersSec(jj);
                obj=obj.delsample('Index',find(obj.Time>=(trigTime-timeBefore) & obj.Time<=(trigTime+timeAfter)));
            end
         end
        
        
        function [newObj,timePeak,valPeak]=getNormalizedEyelidMatrix(wvfEyelidNorm,triggersIdx)      
            if ~exist('triggersIdx') || isempty(triggersIdx)
                triggersIdx=[1:size(wvfEyelidNorm.Data,2)];
            end
            wvfEyelidData=wvfEyelidNorm.Data;
            timeVect=wvfEyelidNorm.Time;
            wvfEyelidDataUSCSUS=wvfEyelidData(:,triggersIdx);
            threshCR=max(median(wvfEyelidDataUSCSUS'))*.1;
            idxNOCR=find(max(wvfEyelidDataUSCSUS(abs(timeVect)<0.02,:)<threshCR));
            avgCSUS=median(wvfEyelidDataUSCSUS(:,idxNOCR),2);
            timeSubIdx=find(timeVect>0.04 & timeVect<.15);
            [valPeak,idxPeak]=max(avgCSUS(timeSubIdx,:));
            timePeak=timeVect(timeSubIdx(idxPeak));
            newObj=wvfEyelidNorm/valPeak;
        end
        
      function [amplitudesAll,eyelidMatrixNorm,idxMinAll]=getCRAmplitudes(eyelidMatrixNorm,timePeaks,lengthFilterSec,triggersIdx)
          if ~exist('triggersIdx') 
              triggersIdx=[1:size(eyelidMatrixNorm.Data,2)];
          end
          amplitudesAll=[];
          idxMinAll=[];
          for jj=1:numel(timePeaks)
              timePeak=timePeaks(jj);
              lengthFilter=round(lengthFilterSec/eyelidMatrixNorm.frameRate);
              eyelidMatrixNorm=eyelidMatrixNorm.filter(ones(1,lengthFilter)/lengthFilter,1);
              timeVect=eyelidMatrixNorm.Time;
              wvfEyelidData=eyelidMatrixNorm.Data;
              [valMin,idxMin]=min(abs(timeVect-timePeak));
              amplitudes=wvfEyelidData(idxMin,triggersIdx);
              amplitudesAll=[amplitudesAll; amplitudes];
              idxMinAll=[idxMinAll; idxMin];
          end                   
      end

      function [triggersOut,namesTriggersOut,amplitudesOut,filteredEyelid,amplitudesBefore,amplitudesAtUS,idxToRemove]=getCRsAnalysis(obj,idxUS,idxCSUS,idxCS,numIQRCRs,ISI,idxIsRunning,doExcludeNegativeCRs,timeSupposedMax)
          if ~exist('timeSupposedMax')
              timeSupposedMax=.03;              
          end
          filterLengthSec=.15;
          timeBefore=2*ISI;
          timePeakUS=.14;                    

          wvfEyelidNorm=obj;
          
          [amplitudesAtUS,filteredEyelid]=getCRAmplitudes(wvfEyelidNorm,timeSupposedMax,filterLengthSec);
          [amplitudesAtCR]=getCRAmplitudes(wvfEyelidNorm,timePeakUS,filterLengthSec);

          amplitudesBefore=getCRAmplitudes(wvfEyelidNorm,-timeBefore-ISI-timeSupposedMax,filterLengthSec);
          
          amplCRsCSUS=getCRAmplitudes(wvfEyelidNorm,timeSupposedMax,filterLengthSec,idxCSUS);
          amplCRsCS=getCRAmplitudes(wvfEyelidNorm,timePeakUS,filterLengthSec,idxCS);
          amplUS=getCRAmplitudes(wvfEyelidNorm,timePeakUS,filterLengthSec,idxUS);
          
          if numIQRCRs>1
              thresholdCRs=median(amplitudesBefore)+numIQRCRs*iqr(amplitudesBefore);
          else
              thresholdCRs=numIQRCRs;
          end
          idxCSUSwCR=idxCSUS(find(amplCRsCSUS>thresholdCRs));
          idxCSUSwoCR=idxCSUS(find(amplCRsCSUS<=thresholdCRs));
          idxCSwCR=idxCS(find(amplCRsCS>thresholdCRs));
          idxCSwoCR=idxCS(find(amplCRsCS<=thresholdCRs));
          
          if (doExcludeNegativeCRs)
              idxCRNegative=find(amplitudesAtUS<=-thresholdCRs);
          else
              idxCRNegative=[];
          end
                    
          idxToRemove=union(idxCRNegative,idxIsRunning);
                    
          idxCSUSwCR=setdiff(idxCSUSwCR,idxToRemove);
          idxCSUSwoCR=setdiff(idxCSUSwoCR,idxToRemove);
          idxCSwCR=setdiff(idxCSwCR,idxToRemove);
          idxCSwoCR=setdiff(idxCSwoCR,idxToRemove);
          idxUS=setdiff(idxUS,idxToRemove);
          
          amplCSUSwCR=amplitudesAtUS(idxCSUSwCR);
          amplCSUSwoCR=amplitudesAtUS(idxCSUSwoCR);
          amplCSwCR=amplitudesAtCR(idxCSwCR);
          amplCSwoCR=amplitudesAtCR(idxCSwoCR);
          amplUS=amplitudesAtCR(idxUS);
          
          
%           idxCSUS=[idxCSUSwCR idxCSUSwoCR];
%           idxCS=[idxCSwCR idxCSwoCR];
          
          triggersOut{1}=idxUS;
          triggersOut{2}=idxCSUSwoCR;
          triggersOut{3}=idxCSUSwCR;
          triggersOut{4}=idxCSwCR;
          triggersOut{5}=idxCSwoCR;
          

          amplitudesOut{1}=amplUS;
          amplitudesOut{2}=amplCSUSwoCR;
          amplitudesOut{3}=amplCSUSwCR;
          amplitudesOut{4}=amplCSwCR;
          amplitudesOut{5}=amplCSwoCR;
            
          namesTriggersOut={'US','CSUSNoCR','CSUSwCR','CSwCR','CSNoCR'};

            
      end
      
     
    end
end