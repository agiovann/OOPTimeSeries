classdef mytscollection < tscollection
    % The following properties can be set only by class methods
    properties (SetAccess = public)
        frameRate
    end
    methods (Static)
        function obj=createCollectionFromTSCollection(ts)
            obj=mytimeseries;
            obj=obj.fillWithTsCollection(coll,ts);
        end
        function fixPlotForAI(filename)
            hold on
            xl=get(gca,'xlim');
            xl=xl(1)+.5;
            yl=get(gca,'ylim');
            yl=yl(1)+.5;
            plot([xl xl+2],[yl yl],'k','LineWidth',2)
            %            text(xl+1,yl+1,'2s')
            plot([xl xl],[yl yl+1],'k','LineWidth',2)
            %            text(xl+2,yl+2,'100% \DeltaF/F')
            axis off
%             axis auto
            saveas(gcf,filename,'epsc')
        end
    end
    
    methods (Access = protected)
        
    end % methods
    methods (Access = public)
        function obj = mytscollection(varargin)
            % Constructor getting a cell ararray of movie names, join them
            % togethrt and set the variables inside the class
            if(isempty(varargin))
                error('Provide Time Vector At Least')
            end
            argConstr=varargin;
            frameRate=[];
            if(~isempty(varargin) && isfloat(varargin{1}))
                if(nargin == 2)
                    matrix=varargin{1};
                    timeVect=varargin{2};
                    ARGc={};
                    for numTr=1:size(matrix,1)
                        ARGc{numTr}=mytimeseries(matrix(numTr,:)',timeVect,'Name',['D_' num2str(numTr)]);
                    end
                    argConstr={};
                    argConstr{1}=ARGc;
                else
                    argConstr=varargin(1);
                end
            end
            if(~isempty(varargin) && isequal( class(varargin{1}),'tscollection'))
                ts=varargin{1};
                counter=0;
                ARGc={};
                for names=ts.gettimeseriesnames
                    myts= ts.get(names);
                    counter=counter+1;
                    ARGc{counter}= mytimeseries(myts.Data,myts.Time,'Name',['D_' num2str(counter)]);
                end
                argConstr= {};
                argConstr{1}=	ARGc;
            end
            obj=obj@tscollection(argConstr{:});
            obj.frameRate=nanmedian(diff(obj.Time));
        end
        
        
        function obj = fillWithTsCollection(obj,ts)
            for names=ts.gettimeseriesnames
                myts= ts.get(names);
                myts= mytimeseries(myts.Data,myts.Time);
                obj=obj.addts(myts,names);
            end
            obj.frameRate=nanmedian(diff(obj.Time));
        end
        
        function h = plotAllLines( obj, varargin)
            names  = gettimeseriesnames(obj);
            colors=colormap(lines);
            colors =cat(1,[1 1 1],colors);
            hold on
            hs=[];
            for n=1:length(names)
                ts = obj.get(names{n});
                if ~isempty(ts.Data)
                    h=plot(ts,'Color',colors(n+1,:),varargin{:})
                    hs(n)=h(1);
                else
                    hs(n)=plot(0,0,'w');
                end
                
            end
            legend(hs,names)
        end
        
        function newobj=resampleWithSpikes(obj,timeVect)
              names  = gettimeseriesnames(obj);
              newobj=mytscollection(timeVect);
              frate=median(diff(timeVect));
              frateRatio=obj.frameRate/frate;
              
              for n=1:length(names)
                      ts = obj.get(names{n});        
                      spikes=ts.Data;
                      ispikes=find(spikes>0);
                      newspikes=zeros(size(timeVect));
                      newspikes(round(ispikes*frateRatio)-1)=1;                      
                      newobj=newobj.addts(mytimeseries(newspikes,timeVect),names(n));
              end
        end
       

        function h = plot( obj, varargin)
            % TSC_PLOT - plot TSCollection
            % example:
            %  tsc_plot(tsc, 'linewidth',2)
            %             h = figure;
            names  = gettimeseriesnames(obj);
            colors=colormap(lines);
            colors =cat(1,[1 1 1],colors);
            recentMax=0;
            liness=[];
            dataMat=obj.Data;
            minDist=nanmean(dataMat(:));
            for n=1:length(names)
                ts = obj.get(names{n});               
                if ~isnan(max(range(ts.Data)))
                     recentMax=recentMax+max(minDist,max(range(ts.Data)));
                     allLines=plot(obj.Time, (recentMax+ts.Data),  'Color',colors(n+1,:), varargin{:});
                     liness(n)=allLines(1);
                    
                else
                    recentMax=recentMax+max(minDist,recentMax/n);
                    allLines=plot(obj.Time,repmat(recentMax,size(obj.Time)),'Color',colors(n+1,:), varargin{:});
                    liness(n)=allLines(1);
                end
                hold on
            end
            hold off
            xlabel('Time (s)')
            axis tight
            %             set(gca,'YTick',[])
            legend(liness, names,'interpreter','none','Fontsize',8,'Location','East');
            %             title(obj.name);
        end
        
        function h = plotTracesAndEvents( obj,events, varargin)
            % TSC_PLOT - plot TSCollection
            % example:
            %  tsc_plot(tsc, 'linewidth',2)
            %       
            if isequal(class(events),'mytscollection')
                events=events.resample(obj.Time);
                events=(events.Data')>0;                
            end
            h = figure;
            if nargin==3
                trigsSecs=varargin{1};
            else
                trigsSecs=[];
            end
            names  = gettimeseriesnames(obj);
            %             colors = hsv(length(names));
            colors = lines(length(names));
            
            colors(2,:)=[0 0 0];
            recentMax=0;
            linesH=[];
            numiqrs=4;
            for n=length(names):-1:1
                ts = obj.get(names{n});
                hold on
                evmod=events(n,events(n,:)>0).*(recentMax+numiqrs*iqr(ts.Data));
                timemod=obj.Time(events(n,:)>0);
                factSpike=quantile(ts.Data,.95)/5;
                plot([timemod timemod]',[evmod' evmod'-factSpike]','k','linewidth',1.5)
                linesH(n)=plot(obj.Time, (recentMax+ts.Data), 'Color',colors(n,:), varargin{2:end});
%                 recentMax=recentMax+max(range(ts.Data));
                recentMax=recentMax+quantile((ts.Data),.999);
            end
            hold off
            xlabel('Time (s)')
            axis tight
            yL = get(gca,'YLim');
            if ~isempty(trigsSecs)
                line([trigsSecs trigsSecs],yL,'Color',[.8 .8 .8],'LineWidth',2);
                %             set(gca,'YTick',[])
            end
            legend(linesH, names,'interpreter','none','Fontsize',8,'Location','NorthEastOutside');
            %             title(obj.name);
        end
        
        
        function mts = plotMEDIAN( obj, varargin)
            % TSC_PLOT - plot TSCollection
            % example:
            %  tsc_plot(tsc, 'linewidth',2)
            %             h = figure;
            names  = gettimeseriesnames(obj);
            %             colors = hsv(length(names));
            colors=colormap(lines);
            colors =cat(1,[1 1 1],colors);
            mts=mytscollection(obj.Time);
            for n=1:length(names)
                ts = obj.get(names{n});
                plot(obj.Time, nanmedian(ts.Data,2),'Color', colors(n+1,:), varargin{:});
                mts=mts.addts(nanmedian(ts.Data,2),['D_' num2str(n)]);
                hold on
            end
            hold off
            xlabel('Time (s)')
            legend(gca, names,'interpreter','none','Fontsize',8);
            %             title(obj.name);
        end
        
        
        
        function h = plotMEAN( obj, varargin)
            % TSC_PLOT - plot TSCollection
            % example:
            %  tsc_plot(tsc, 'linewidth',2)
            %             h = figure;
            names  = gettimeseriesnames(obj);
            colors=colormap(lines);
            colors =cat(1,[1 1 1],colors);
            h=[];
            for n=1:length(names)
                ts = obj.get(names{n});
                if ~isempty(ts.Data)
                    h(n)=plot(obj.Time, nanmean(ts.Data,2),'Color', colors(n+1,:), varargin{:});
                else
                    h(n)=plot(0,0,'w');
                end
                hold on
            end
            hold off
            xlabel('Time (s)')
            legend(h, names,'interpreter','none','Fontsize',8);
            %             title(obj.name);
        end
        
        
        function newObj=filter(obj,b,a)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts=filter(ts,b,a);
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=decimate(obj,decimationFactor)
            newObj=mytscollection(downsample(obj.Time,decimationFactor));
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts=decimate(ts,decimationFactor);
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=diff(obj)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts=mytimeseries(diff([0; ts.Data]),ts.Time);
                newObj=newObj.addts(ts,names);
            end            
        end
        
        function newObj=getSpikeShape(obj,eventsC,timeBef,timeAft)
            names  = gettimeseriesnames(obj);
            ts=getSpikeShape(obj.get(names{1}),eventsC.get(names{1}),timeBef,timeAft);
            newObj=mytscollection(ts.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);                
                ts=getSpikeShape(ts,eventsC.get(names),timeBef,timeAft);
                newObj=newObj.addts(ts,names);
            end            
        end
        
        
        function newObj=nanmeanAG(obj,dimension)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts=nanmeanAG(ts,dimension);
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=nanmedianAG(obj,dimension)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts=nanmedianAG(ts,dimension);
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=minus(obj,obj2)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                switch class(obj2)
                    case 'mytscollection'
                        ts2= obj2.get(names);
                        ts=ts-ts2;
                    otherwise
                        ts=ts-obj2;
                end
                newObj=newObj.addts(ts,names);
            end
        end
        
        
        function newObj=plus(obj,obj2)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                switch class(obj2)
                    case 'mytscollection'
                        ts2= obj2.get(names);
                        ts=ts+ts2;
                    otherwise
                        ts=ts+obj2;
                end
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=rdivide(obj,obj2)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                
                switch class(obj2)
                    case 'mytscollection'
                        ts2= obj2.get(names);
                        ts=ts./ts2;
                    otherwise
                        ts=ts./obj2;
                end
                
                
                
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=horzcat(obj,obj2)
            if  isempty(obj)
                newObj=obj2;
            else
                if ~isempty(obj2)
%                     obj2.Time=obj2.Time+(obj.Time(end)-obj.Time(1));
                      obj2.Time=obj2.Time+(obj.Time(end)-obj.Time(1)+obj.frameRate);

                    newObj=[];
                    for names  = gettimeseriesnames(obj);
                        ts = obj.get(names);
                        ts2= obj2.get(names);
                        if isempty(newObj)
                            newObj=mytscollection(ts.append(ts2).Time);
                        end
                        newObj=newObj.addts(ts.append(ts2),names);
                    end
                else
                    newObj=obj;
                end
            end
        end
        
        function newObj=vertcat(obj,obj2)
            if ~isempty(obj2)
                newObj=mytscollection(obj.Time);
                for names  = gettimeseriesnames(obj);
                    ts = obj.get(names);
                    ts2= obj2.get(names);
                    nts=mytimeseries([ts.Data ts2.Data],ts.Time);
                    newObj=newObj.addts(nts,names);
                    
                end
            else
                newObj=obj;
            end
        end
        
        
        function [newObj]=extractWaveforms(obj,trigsSecs,secBefore,secAfter)
            newObj=tscollection;
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts = extractWaveforms(ts,trigsSecs,secBefore,secAfter);
                newObj=newObj.addts(ts,names);
            end
            newObj=mytscollection(newObj);
        end
        function newObj=getDFFQuantile(obj,varargin)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= getDFFQuantile(ts,varargin{:});
                newObj=newObj.addts(ts,names);
            end
            
            
        end
        
        function newObj=removeBLQuantile(obj,varargin)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= removeBLQuantile(ts,varargin{:});
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=normalizeIQR(obj)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= mytimeseries(ts.Data./iqr(ts.Data),obj.Time);
                newObj=newObj.addts(ts,names);
            end
        end
        
        
         function obj=getTraceWithoutPeriTrigger(obj,timeBefore,timeAfter,triggersSec)
            for jj=1:length(triggersSec)
                trigTime=triggersSec(jj);
                obj=obj.delsamplefromcollection('Index',find(obj.Time>=(trigTime-timeBefore) & obj.Time<=(trigTime+timeAfter)));
            end
         end
        
         
         function [amplitudesTot,matNorm]=getCRAmplitudes(obj,timePeaks,lengthFilterSec,triggersIdx)         
             matNorm=mytscollection(obj.Time);
             if ~iscell(triggersIdx)
                 amplitudesTot=struct;
                 for names  = gettimeseriesnames(obj);
                     [amps,eyelidMatrixNorm,~]=getCRAmplitudes(obj.get(names),timePeaks,lengthFilterSec,triggersIdx);
                     amplitudesTot.(names{1})=amps;
                      matNorm= matNorm.addts(eyelidMatrixNorm,names);
                 end
             else
                 amplitudesTot={};
                 for gg=1:numel(triggersIdx)
                     amplitudesAll=struct;
                     for names  = gettimeseriesnames(obj);
                         [amps,~,~]=getCRAmplitudes(obj.get(names),timePeaks,lengthFilterSec,triggersIdx{gg});
                         amplitudesAll.(names{1})=amps;
                        
                     end
                     amplitudesTot{gg}=amplitudesAll;
                 end
             end
             
         end
         
         function newObj=clip(obj,minVal,maxVal,multFact)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);                
                ts = obj.get(names);
                dat=ts.Data;
%                 dat(dat<minVal)=minVal;
%                 dat(dat>maxVal)=maxVal;
                dat=agClip(dat,minVal,maxVal,multFact);
                ts= mytimeseries(dat,obj.Time);
                newObj=newObj.addts(ts,names);
            end             
         end
             
         
         function obj=getTraceOnlyPeriTrigger(obj,timeBefore,timeAfter,triggersSec)            
            indexes=round(triggersSec/obj.frameRate);
            samplesBefore=round(timeBefore/obj.frameRate);
            samplesAfter=round(timeAfter/obj.frameRate);
            vectorLength=length(obj);
            indexesTrigger=agAddRectangleAroundIndexes(indexes,samplesBefore,samplesAfter,vectorLength);           
            obj=obj.delsamplefromcollection('Index',find(1-indexesTrigger));           
        end
        
        function mat=Data(obj)
            mat=[];
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                mat=[mat ts.Data];
            end
        end
        function [eventsC,eventsMQuality,templatesAll]=extractCalciumEventsMRWithTriplets(obj,minProbabAmpl,minProbabDec,isSinglet,varargin)
            [eventsM,eventsMQuality,templatesAll]=extractCalciumEventsMRWithTriplets(obj.Data',obj.frameRate,minProbabAmpl,minProbabDec,isSinglet);
            eventsC=mytscollection(obj.Time);
            tnames=obj.gettimeseriesnames;
            for kk=1:numel(tnames)
                eventsC=eventsC.addts(mytimeseries(eventsM(kk,:)',obj.Time),tnames{kk});
            end
        end
        
        function  eventsC=extractCalciumEventsFilterOEM(obj,iqrMultiplier)
            largerThanSurroundFraction=0;
            frate=median(diff(obj.Time));
            eventsM=agExtractCacliumEventsFilterOEM(obj.Data',frate,iqrMultiplier,largerThanSurroundFraction);
            eventsC=mytscollection(obj.Time);
            tnames=obj.gettimeseriesnames;
            for kk=1:numel(tnames)
                eventsC=eventsC.addts(mytimeseries(eventsM(kk,:)',obj.Time),tnames{kk});
            end
        end
        
        function  [eventsC,myNlogL]=extractCalciumEventsDeconvOrig(obj,probSpikes)
            [eventsM,~,~,myNlogL,~,~,~,~,~]=extractCalciumEvents01(obj.Data',obj.frameRate,probSpikes);            
             eventsC=mytscollection(obj.Time);
            tnames=obj.gettimeseriesnames;
            for kk=1:numel(tnames)
                eventsC=eventsC.addts(mytimeseries(eventsM(kk,:)',obj.Time),tnames{kk});
            end
        end
        
        
        function  [eventsC,decSig,CaSig]=extractCalciumEventsFastNonnegDec(obj,expFreq,numStdNoise,minPeakDec)
            tau=.17;            
            [decSig,P_best,Vv,CaSig]=fastNonegativeDeconvolution(obj,tau,expFreq,numStdNoise);           
            names=decSig.gettimeseriesnames();
            eventsC=mytscollection(obj.Time);
            for nm=names
               ds=decSig.get(nm).Data;
               if ~isempty(minPeakDec)
                [probSpikes,idx]=findpeaks(ds,'MINPEAKHEIGHT',minPeakDec);
               else
                [probSpikes,idx]=findpeaks(ds);
               end
               ev=zeros(size(ds));
               ev(idx)=1;               
               eventsC=eventsC.addts(mytimeseries(ev,eventsC.Time),nm);
            end
        end
        
        function gt(obj,val)
           obj.Data=obj.Data>val; 
        end
        function  [decSig,P_best,Vv,CaSig]=fastNonegativeDeconvolution(obj,tau,expFreq,numStdNoise,varargin)           
            V.dt=obj.frameRate;
            P=struct;
            if ~isempty(tau)
                P.gam=1-obj.frameRate./tau;
            end
            if ~isempty(expFreq)
                P.lam=expFreq;
            end
%             P.b=0;
%               P.a=5e-9; 
%               P.a=1e-1;
%             V.fast_poiss=0;
%             V.est_sig=1;
%             V.est_gam=1;
%             V.est_lam=1;
%             V.est_a=1;
%             V.est_b=1;
%             obj=obj.removeBLQuantile(.3,10);
            nms=obj.gettimeseriesnames;
            decSig=mytscollection(obj.Time);
            CaSig=mytscollection(obj.Time);
            P_best={};
            Vv={};
%             if ~isempty(numStdNoise)
%                 P.sig=quantile(iqr(obj.Data)/1.3*numStdNoise,.75);
%             end
            
            
            for ll=1:numel(nms)
                tr=obj.get(nms{ll});
%                 P.b=quantile(tr.Data,0.08);
%                 P.b=0;
                if ~isempty(numStdNoise)
                    P.sig=std(tr)*numStdNoise;                
                end
                [n_best,P_best{ll},Vv{ll},Cc]=ag_fast_oopsi(tr.Data,V,P);
                decSig=decSig.addts(mytimeseries(n_best,obj.Time),nms(ll));
                CaSig=CaSig.addts(mytimeseries(Cc,obj.Time),nms(ll));
            end
        end
        
        function newObj=mtimes(obj,val)
            newObj=mytscollection(obj.Time); 
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= ts*val;
                newObj=newObj.addts(ts,names);
            end

        end
        function newObj=nonlinearityRemove(obj,varargin)
            newObj=mytscollection(obj.Time);
           
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= ts.nonlinearityRemove(varargin{:});
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newobj=replaceSamplesFromCollection(obj,indexes,namefunc)
            if isempty(namefunc)
                namefunc='norm';                
            end
            
            switch namefunc
                case 'norm'
                    nms=obj.gettimeseriesnames;
                    newobj=mytscollection(obj.Time);
                    for jj=1:length(nms)
                        ts=obj.get(nms{jj});
                        datt=ts.Data;
                        stdd=iqr(datt)/20;
                        meand=median(datt);
                        datt(indexes)=normrnd(meand, stdd,size(indexes)); 
                        ts.Data=datt;
                        newobj=newobj.addts(ts,ts.Name);
                    end 
                case 'nan'
                    nms=obj.gettimeseriesnames;
                    newobj=mytscollection(obj.Time);
                    for jj=1:length(nms)
                        ts=obj.get(nms{jj});
                        datt=ts.Data;                        
                        datt(indexes)=NaN;
                        ts.Data=datt;
                        newobj=newobj.addts(ts,ts.Name);
                    end
                    
                otherwise
                    error('unknown funcion name')
            end
            
            
        end
        function newObj=nonlinearityRemoveSigmoid(obj,indicator,basalFluorescence)
            newObj=mytscollection(obj.Time);
           
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= ts.nonlinearityRemoveSigmoid(indicator,basalFluorescence);
                newObj=newObj.addts(ts,names);
            end
         end
        
        function saveToDir(obj,name)
            save(['./' name],'obj');
        end
        
        
        
        function newObj=decimateSpikes(obj,decFactor)
            newObj=mytscollection(obj.Time(1:decFactor:end));
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= ts.decimateSpikes(decFactor);
                newObj=newObj.addts(ts,names);
            end
        end
        
         
        
        function newObj=deconv(obj,tauFilt,durFilt)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts= ts.deconv(tauFilt,durFilt);
                newObj=newObj.addts(ts,names);
            end
        end
        
        function newObj=removeBLWaveforms(obj,timeWinBefore,timeInitTrigger)
            if ~exist('timeInitTrigger')
               timeInitTrigger=0; 
            end
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts = ts.removeBLWaveforms(timeWinBefore,timeInitTrigger);
                newObj=newObj.addts(ts,names);
            end            
        end
        
        
        
        function newObj=getDFFQuantileWVF(obj,timeWinBefore,timeTrigger)
            newObj=mytscollection(obj.Time);
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                ts = ts.getDFFQuantileWVF(timeWinBefore,timeTrigger);
                newObj=newObj.addts(ts,names);
            end
            
        end
        
        
        function newObj=getSubCollection(obj,idxs)
            newObj=mytscollection(obj.Time);

            if iscell(idxs)                
                for kk=1:numel(idxs)
                    nms=idxs(kk);
                    newObj=newObj.addts(obj.get(nms));
                end
            else
                names  = gettimeseriesnames(obj);
                for kk=1:numel(idxs)
                    newObj=newObj.addts(obj.get(names(idxs(kk))));
                end
            end
            
%             allIdxs=1:numel(names);
%             complIdx=setdiff(allIdxs,idxs);
%             newObj=obj.removets(names(complIdx));
        end
        
        function rastersTot=rasterPlots(obj)
            nms=obj.gettimeseriesnames;
            rowcols=ceil(sqrt(numel(nms)));
            rastersTot=[];
            for nn=1:numel(nms)
                subplot(rowcols,rowcols,nn)  
                dat=obj.get(nms{nn});
                imagesc([obj.Time(1) obj.Time(end)],[],dat.Data')
                rastersTot=cat(3,rastersTot,dat.Data');
            end
        end
        
        function newobj=meanAllTS(obj)
            newdata=[];
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                newdata=[newdata ts.Data];
            end
            newobj=mytimeseries(nanmean(newdata,2),obj.Time);
        end
        
        
       
       function grpdTraces=plotGroupedByTriggersWvf(obj,allTriggersStr,trigsSecs,timeBeforeForBLremoval,func)
           if ~exist('func')
               func=@(X) nanmean(X,2);
           end
           
           dat=obj.Data;
           dat=(feval(func,dat));
           if ~isempty(dat)
               joinedTrace=mytimeseries(dat,obj.Time);
               jtWVF=joinedTrace.extractWaveforms(trigsSecs,1.7,2);       
               jtWVF=jtWVF.removeBLWaveforms(.5,timeBeforeForBLremoval);
               grpdTraces=jtWVF.getSubCollectionGroups(allTriggersStr);
               plotMEAN(grpdTraces);
           else
               grpdTraces=[];
           end
               
        end
        
        %         function obj=addts(obj,varargin)
        %             obj=addts@tscollection(obj,varargin{:});
        % %             if isempty(obj.frameRate)
        % %                 ts=varargin{1};
        % %                 obj.frameRate=ts.frameRate;
        % %             end
        %         end
        function kinetics=computeRaiseAndDecayTime(obj,minPerc,maxPerc,minDFF,stimLength)
            kinetics=[]; 
            for names  = gettimeseriesnames(obj);
                ts = obj.get(names);
                kinetics= [kinetics ts.computeRaiseAndDecayTime(minPerc,maxPerc,minDFF,stimLength)];               
            end
        end
    end
end