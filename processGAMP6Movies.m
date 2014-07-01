function [movClOrig,movClNorm,movCLTrig,movCLMagn]=processGAMP6Movies(movFile,numFramesToSkipBeginning,isMotCorr,timeMotCorr)
% %%
% clear movClNorm trNorm trOrig trCOrig
% clear movClOrig trBL trDFF ans trC trtscNorm
% %% OR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movClNorm=MovieTSImaging(movNorm.mov,(0:(movImaging.numFrames-1))*movImaging.frameRate,'');
% movClOrig=MovieTSImaging(movImaging.mov,(0:(movImaging.numFrames-1))*movImaging.frameRate,'');
% movCLTrig=MovieTSTrigger(movTrig.mov,(0:(movTrig.numFrames-1))*movImaging.frameRate,'');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
    isMotCorr=1;
end
imgChannel=1;
trigChannel=2;
magnChannel=3;
    % clear all classes
    % movFile={'slight movement se -433 -386 95  face stim005.tif'}
    %mov={'-1015 -266 24 16 lines pu cell001.tif'};
    [m,t]=MovieTSImaging.loadMovieChain(movFile,imgChannel,numFramesToSkipBeginning);
    movClOrig=MovieTSImaging(m,t,movFile);
    if(isMotCorr)
        try
            [movClOrigMC,ns]=movClOrig.motionCorrect(movClOrig.smooth3('gaussian',[3,3, round(timeMotCorr/movClOrig.frameRate)],1),0,8);
            if numel(find(ns<10))>(movClOrigMC.Length/20)
                            [movClOrigMC,ns]=movClOrig.motionCorrect(movClOrig.smooth3('gaussian',[3,3, round(3/movClOrig.frameRate)],.85),0,8);

            end
        catch
            try
                [movClOrigMC,ns]=movClOrig.motionCorrect(movClOrig.smooth3('gaussian',[3 3 round(3.5/movClOrig.frameRate)],.65),0,8);
                if numel(find(ns<10))>(movClOrigMC.Length/20)
                    [movClOrigMC,ns]=movClOrig.motionCorrect(movClOrig.smooth3('gaussian',[3 3 round(7/movClOrig.frameRate)],.65),1,8);
                end
            catch
                movClOrigMC=movClOrig;
            end
        end
        movClOrig=movClOrigMC;
    end
%     movClNorm=movClOrig.normalizeMovieQuantile(.5,5);
    movClNorm=(movClOrig-movClOrig.zProject(.5));
    %     movClNorm=movClOrig.normalizeMovie(10);    
try
    [m,t]=MovieTSTrigger.loadMovieChain(movFile,trigChannel,numFramesToSkipBeginning);
    movCLTrig=MovieTSTrigger(m,t,movFile);
catch e
    warning('NO TRIGGERS IN FILE!')
    pause(2)
    movCLTrig=MovieTSTrigger(zeros(size(movClOrig.Data)),movClOrig.Time,movFile);   
end
try
    [m,t]=MovieTSMagnet.loadMovieChain(movFile,magnChannel,numFramesToSkipBeginning);
    movCLMagn=MovieTSMagnet(m,t,movFile);
catch e
    warning('NO MAGNET IN FILE!')
    pause(2)
    movCLMagn=MovieTSMagnet(zeros(size(movClOrig.Data)),movClOrig.Time,movFile);   
end
%%
% close all
% imagesc(movClOrig.zProject)
% colormap gray
% movClOrig.saveFig(gcf,'avgMovie','fig')
% close
%%
% plot(movClOrig)
% movClOrig.saveFig(gcf,'AvgFOVTrace','fig')
% close
%%
% movClOrig.zProject(.5);
% axis image;
% movClOrig.saveFig(gcf,'zProject','fig')
% close
%%
% plot(movClNorm) 
% movClOrig.saveFig(gcf,'AvgFOVTraceNorm','fig')
% close
%%
movClOrig.saveToDir('movOrig');
movClNorm.saveToDir('movNorm')
movCLTrig.dirFiles=movClNorm.dirFiles;
movCLTrig.saveToDir('movTrig');
movCLMagn.dirFiles=movClNorm.dirFiles;
movCLMagn.saveToDir('movMagn');
