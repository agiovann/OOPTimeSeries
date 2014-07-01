classdef MovieTSTrigger < MovieTS
    properties (SetAccess = public)
        
    end
    
    methods (Static)
        
    end
    
    methods (Access = protected)
        
    end
    
    methods (Access = public)
        function obj = MovieTSTrigger(movMat,timeVect,movNames)
            obj=obj@MovieTS(movMat,timeVect,movNames);
        end
        
        %         function [trigsSecs,trigSamp]=extractTriggers(obj,minDistTrigsSeconds)
        %              minDistTrigs=ceil(minDistTrigsSeconds/obj.frameRate);
        %              maxTS=max(obj);
        %             [allTriggersAlignedUS,~,~,~]=agExtractTriggers(maxTS.Data,minDistTrigs);
        %             trigsSecs=(allTriggersAlignedUS-1)*obj.frameRate;
        %             trigSamp=allTriggersAlignedUS;
        %         end
        
        function [triggers,namesTriggers]=extractTriggersType(obj,trigsSec,varargin)
            maxTypes=5;
            maxTS=(obj.getScanSignal);
            maxTS=mytimeseries(maxTS.resample(0:.001:obj.Time(end),'zoh'));
            wvf=maxTS.extractWaveforms(trigsSec,0,1);
            trigMat=wvf.Data';
            trigMat(isnan(trigMat))=0;
            trigMat=trigMat./max(trigMat(:));
            [coeff,score,latent,tsquared,explained,mu]=pca(trigMat);
            clusters=[];
            if nargin==3
                trigTypes=varargin{1};
                [T,C,sumd]=kmeans(score(:,1),trigTypes,'distance','cityblock','replicates',50);
                clusters=T;
            else
                previous=Inf;
                for tt=1:maxTypes
                    try
                        [T,C,sumd]=kmeans(score(:,1),tt,'distance','cityblock','replicates',50);
                        if (previous-mean(sumd'))>1
                            clusters=T;
                            previous=mean(sumd');
                        else
                            break
                        end
                        previous=mean(sumd');
                    catch
                        break
                    end
                end
            end
            classes=unique(clusters);
            triggers={};
            
            for kk=1:length(classes)
                triggers{kk}=find(clusters==classes(kk));
                namesTriggers{kk}=num2str(classes(kk));
            end
        end
        
        function [trigsSecs,trigSamp]=extractHRTriggers(obj,minDistTrigsSeconds)
            downSampleFactor=10;
            minDistTrigs=ceil(minDistTrigsSeconds/obj.frameRate*obj.linesPerFrame*obj.pixelsPerLine/downSampleFactor);
            maxTS=getScanSignal(obj);
%             timeOrig=maxTS.Time;
            maxTS=resample(maxTS,maxTS.Time(1:downSampleFactor:end));
            [allTriggersAlignedUS,~,~,~]=agExtractTriggers(maxTS.Data,minDistTrigs);
            trigsSecs=(allTriggersAlignedUS-1)*obj.frameRate/obj.linesPerFrame/obj.pixelsPerLine*downSampleFactor;
            trigSamp=allTriggersAlignedUS;
        end
        %
        %          trigSignal=movCLTrig.getScanSignal;
        %     trigSignal=mytimeseries(trigSignal.Data,trigSignal.Time);
        %     trigSignal=trigSignal/max(trigSignal)*.1;
        %     [~,trigsSecs]=movCLTrig.extractHRTriggers(3);
        function [trigsSecs,triggers]=extractTriggersPulsesBen(obj)
            % extract triggers with Ben's algorithm
            slopeThreshold = 0.3; %when undifferentiated signal was normalized
            pulseSeparationThreshold = 200; %samples
            triggerSeparationThreshold = 10000; %samples
            
            movTrig=obj;
            
            frameRate = movTrig.frameRate;
            scanSignal = movTrig.getScanSignal;
            timePoints = scanSignal.Time;
            
            %extract rises and falls:
            signal = scanSignal.Data;
            signal = signal./max(signal);
            signalDiff = diff(signal);
            risefallIdxs = find(abs(signalDiff) > slopeThreshold);
            riseIdxs = risefallIdxs(signalDiff(risefallIdxs) > 0)';
            fallIdxs = risefallIdxs(signalDiff(risefallIdxs) < 0)';
            
            %find rise starts:
            riseIdxs = fliplr(riseIdxs);
            riseIdxs = [riseIdxs(find(abs(diff(riseIdxs)) > pulseSeparationThreshold)), riseIdxs(end)];
            riseIdxs = fliplr(riseIdxs);
            %find fall ends:
            fallIdxs = [fallIdxs(find(diff(fallIdxs) > pulseSeparationThreshold)), fallIdxs(end)];
            
            %separate trigger events:
            triggerBreakIdxs = [1, find(diff(riseIdxs) > triggerSeparationThreshold)+1];
            triggerEvents = {};
            for i=1:numel(triggerBreakIdxs)
                if (i==numel(triggerBreakIdxs))
                    triggerEvents{i} = [riseIdxs(triggerBreakIdxs(i):end); fallIdxs(triggerBreakIdxs(i):end)]';
                else
                    triggerEvents{i} = [riseIdxs(triggerBreakIdxs(i):triggerBreakIdxs(i+1)-1); fallIdxs(triggerBreakIdxs(i):triggerBreakIdxs(i+1)-1)]';
                end
            end
            
            %collect trigger data
            triggers = {};
            trigsSecs = [];
            for i=1:numel(triggerEvents)
                event = triggerEvents{i};
                triggers(i).firstSample = event(1,1);
                triggers(i).lastSample = event(end,end);
                triggers(i).startTime = timePoints ( triggers(i).firstSample );
                triggers(i).endTime = timePoints ( triggers(i).lastSample );
                triggers(i).duration = triggers(i).endTime - triggers(i).startTime;
                triggers(i).numPulses = numel(event(:,1));
                triggers(i).pulseBordersInSamples = event;
                triggers(i).pulseBorders = timePoints ( event );
                triggers(i).pulseDurationsInSamples = event(:,2) - event(:,1);
                triggers(i).pulseDurations = timePoints ( event(:,2) - event(:,1) );
                triggers(i).ipisInSamples = event(2:end,1) - event(1:end-1,2);
                triggers(i).ipis = timePoints ( event(2:end,1) - event(1:end-1,2) );
                
                trigsSecs = [trigsSecs triggers(i).startTime];
            end
        end
        function plot(obj)
            newData=max(obj);
            newData=newData/max(newData);
            plot(newData)
        end
        
        
        function [trigsSecs,triggers]=extractTriggersLargeMovies(obj,minDist, minValuePerc)            
            mm=mean(obj);
            dat=mm.Data;
            dat=dat./max(dat(:));
            diffdat=diff(dat);
            [vals,triggers]=findpeaks(diffdat,'MINPEAKHEIGHT',minValuePerc,'MINPEAKDISTANCE',round(minDist/obj.frameRate))
            trigsSecs=obj.Time(triggers);
            
        end
    end
    
end
