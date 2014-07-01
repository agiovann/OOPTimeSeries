function [movClOrig]=processGAMP6MoviesFast(movFile,numFramesToSkipBeginning,numIQR,timeFromBeginning)
%% %%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgChannel=1;
trigChannel=2;
magnChannel=3;
[m,t]=MovieTSImaging.loadMovieChain(movFile,imgChannel,numFramesToSkipBeginning);
movClOrig=MovieTSImaging(m,t,movFile);
%%
clc
tic
thr=mean(movClOrig.Data(:))+numIQR*std(movClOrig.Data(:));
toc
tic
dendmaskOrig=movClOrig.extractGCMaskNoPCAICA(thr,30,.8,[5 5],.65);
toc
tic
close
traces=movClOrig.extractROIs(dendmaskOrig);
toc
tic
traces=traces.getDFFQuantile(.08,2);
toc
trialDuration=movClOrig.numFramesPerMovie*movClOrig.frameRate;
totalDuration=movClOrig.Length*movClOrig.frameRate;
trigsSecs=timeFromBeginning:trialDuration:totalDuration;
subplot(2,1,1)
plotMEDIAN(traces.extractWaveforms(trigsSecs,timeFromBeginning,2));
subplot(2,1,2)
agcolorfulManyLines(dendmaskOrig);
