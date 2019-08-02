%%%%%%%This code (includes Step1, Step2, Step3)show how cluster individual theta cycle of LFP recorded in hippocampal CA1 pyrmidal layer
%%%%%%%into four clusters:slow gamma; median gamma; early high gamma and late high gamma through k-means clustering
%%%%%%%author: Lu Zhang;math2437@hotmail.com; tested with Matlab 2017b;

%%%%toolbox needed:
%%%%Toolbox 1:buzsaki Lab code for processing hc-11 data; not necessary if use your own
%%%%data:https://github.com/buzsakilab/buzcode

%%%%Toolbox 2:Comunity Clustering;not necessary if only try k-means
%%%%clustering;http://netwiki.amath.unc.edu/GenLouvain/GenLouvain


%%%%%%%Step 1 includes loading hc-11 data;please download the data first if not use your own data;
%% https://crcns.org/data-sets/hc/hc-11/about-hc-11

%%%%%%%Step 2 prepare data samples for clustering; 
%%%%%%%LFP, filtered LFP in theta band,LFP theta phase is needed;

%%%%%%%Step 3-1:k-means clustering and 4 theta-gamma asignment;

%%%%%%%Step 3-2:community clustering;
%%%%%%%%%%%%%%%Noted that this process is time coumsing when data sample is large.
%%%%%%%%%%%%%%%In this demo, only the first 5000 samples were used; It
%%%%%%%%%%%%%%%cost ~1 hour calculation if use all samples in this demo;


%%%%Step 1; prepare the data from hc-11 data set; This is not necesary if
%%%%you have your own data instead.
%%%%samprate=1250;%%sampling rate

%%%%%%%%%Load LFP Chan 30 from example data in hc-11 data set
clear all
% PathFile='Y:\rozell-singer-exchange\CRCNSdata\hc-11\data\Cicero_09012014\Cicero_09012014.eeg';
PathFile=Path to \CRCNSdata\hc-11\data\Cicero_09012014\Cicero_09012014.eeg;

ChanID=68;  %%%%%%load Chan 68, the ChanID here is ChID+1 in Neuroscope;see hc-11 document
ChanNum=134;%%%%%%Total Chan Number in this recording, see hc-11 document
samprate=1250;%%%%%sampling rate of LFP in .eeg file
lfp = bz_LoadBinary(PathFile,'frequency',samprate,'start',0,'duration',Inf,'channels',5,'nChannels',ChanNum);
lfp=double(lfp);
%%%%%%%%%Load LFP Chan 30 from hc-11 data set

%%%%%%%%%filter LFP in theta band and calcuate theta phase
%%%%%%%%%pass band 4-12 Hz?
%%%Creat filter;
Fs = samprate;  % Sampling Frequency
Fstop1 = 3;                % First Stopband Frequency
Fpass1 = 4;                % First Passband Frequency
Fpass2 = 12;               % Second Passband Frequency
Fstop2 = 13;               % Second Stopband Frequency
rs=50;                     % Stopband Attenuation 50DB
rp=1;                      % Passband Ripple 1DB
Dstop1 = 10^(-rs/20);  % First Stopband Attenuation
Dpass  = (10^(rp/20)-1)/(10^(rp/20)+1);   % Passband Ripple
Dstop2 = 10^(-rs/20);  % Second Stopband Attenuation
dens   = 20;               % Density Factor

[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);
bk_thetafilter.tf.num=Hd.Numerator;
bk_thetafilter.tf.den=1;
%%%Creat filter;


    lfpTheta = FiltFiltM(bk_thetafilter.tf.num,bk_thetafilter.tf.den,lfp);
    hilb = hilbert(lfpTheta);
    amp = abs(hilb);
    Thetaphase = angle(hilb);   %%%%%%%%%instanous theta phase;


%%%%Step 2; creat frequency-phase spectrum for each individual theta cycle
WaveParam.ntw=11;      %11 points smoothing time-window 
WaveParam.nsw=3;        %3 points smoothing fre-window 
WaveParam.DownSample=2; %DownSampling parameter;
WaveParam.wname='morl'; %Morlet Wavelet
WaveParam.Samplingrate=samprate;
WaveParam.Freq=[20:2:180]; %Frequeny band interested in Hz
WaveParam.Zscore=1;  %Normalize waveletPower by zscore in time domain for a given freqency point;

%%%%Analysis period in seconds;First row for start;second row for end;
% BehaviorFile='Y:\rozell-singer-exchange\CRCNSdata\hc-11\data\Cicero_09012014\Cicero_09012014_sessInfo.mat';
BehaviorFile=Path to \CRCNSdata\hc-11\data\Cicero_09012014\Cicero_09012014_sessInfo.mat;

load(BehaviorFile);
timerange=getfield(sessInfo.Epochs,'Wake')';   %%%%%Using only the Linear Track Trial ;
%%%%Analysis period from 0 to 40000 seconds;Using large end (15000) could ensure using the whole data.

PhaseBinNum=20; %%%%%%Divid each theta cycle by 20 phase bin;
PhaseStep=2*pi/PhaseBinNum;
PhaseBin=-pi:PhaseStep:pi;
PhaseBin=[PhaseBin(1:end-1)+PhaseBin(2:end)]/2;

ParamOut = SpecThetaExtract(lfp,lfpTheta,Thetaphase,samprate,timerange,WaveParam,PhaseStep);
Div=ParamOut.Div;



%%%%Step 3; Cluster data using either community clustering or k means
%%%%method;

%%%%Step3-1: k-means clustering using correlation distance;
%%%%Adjcent matrix defined by sample correlations;
Cnum=4; %%%%%%%number of clusters 
[SSkmeans,TempCom]=kmeans(ParamOut.sample',Cnum,'distance','correlation','Maxiter',10000);

%%%%%%%average samples within one state;
clear Comkmeans
for iC=1:Cnum
    Comkmeans(:,:,iC)=reshape(nanmean(ParamOut.sample(:,SSkmeans==iC),2),Div(1),Div(2));
end
%%%%%%%average samples within one state;
 figure;
%%%%%visulize the clusters;
for iC=1:Cnum              
                
            subplot(1,Cnum,iC);
            imagesc(PhaseBin,ParamOut.Fplot,Comkmeans(:,:,iC)');axis xy;
            set(gca,'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'})
            if iC==1
               ylabel('Frequency Hz') 
            elseif iC==2
               xlabel('Theta phase rad'); 

            end
            
end




thRateRatio=0.95; %%%%%%%%%0.95PeakValue was set as threshold for gamma field
GammaField=GravitySpec(Comkmeans,thRateRatio,PhaseBin,ParamOut.Fplot);
%%%%%asign four clusters into Slow gamma;Median gamma;Early Fast gamma and
%%%%%Late Fast gamma;
SortI=PhaseFreSort(GammaField.GravityFre,GammaField.GravityPhase);
ComkmeansSorted=Comkmeans(:,:,SortI);
GravityFre=GammaField.GravityFre(SortI);
GravityFre=GammaField.GravityFre(SortI);

SSsorted=SSkmeans;
 for ii=1:Cnum
     SSsorted(SSkmeans==SortI(ii))=ii;
 end
ComColor= [10,103,155;176,34,242;78,211,34;191,126,0;178,51,51;51,204,178]/255;

figure;
[Contx,Conty]=meshgrid(PhaseBin,ParamOut.Fplot);
for iC=1:Cnum              
                
            subplot(1,Cnum,iC);
            imagesc(PhaseBin,ParamOut.Fplot,ComkmeansSorted(:,:,iC)');axis xy;
            Th=max(max(squeeze(ComkmeansSorted(:,:,iC))))*thRateRatio;
            hold on;
            contour(Contx,Conty,ComkmeansSorted(:,:,iC)',[0 Th],'-','color',ComColor(iC,:),'linewidth',2);
            plot(GammaField.GravityPhase(SortI(iC)),GammaField.GravityFre(SortI(iC)),'^','markersize',6,'color',ComColor(iC,:));

            set(gca,'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'})
            if iC==1
               ylabel('Frequency Hz') 
            elseif iC==2
               xlabel('Theta phase rad'); 

            end
            
 end

 
%%%%%SSsorted gives you theta-gamma cluster ID for individual theta trough timestamps represented by ParamOut.ThetaTS 
 

%%%%Step3-2: community clustering;http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%%%%Adjcent matrix defined by sample correlations;
sampleNeedIndex=1:5000; %%%%%%%cluster the first 5000 samples;
SampleCorr=corr(ParamOut.sample(:,sampleNeedIndex));
AdjMat=SampleCorr+1;  %%%%%%%%%%Non-negative weighted correlation.
clear SampleCorr;
for itt=1:size(AdjMat,1)
    AdjMat(itt,itt)=0;%%%%%%diagonal zeros for Adjcent matrix
end
k = full(sum(AdjMat));
twom = sum(k);
%         B = full(Adj - gamma*k'*k/twom);
        tic
        k = full(sum(AdjMat));
        twom = sum(k);
        gamma = 1;  %%Community Clustering Paramter%%
        limit =100000; %%memory consideration for community clustering 

        B = @(i) AdjMat(:,i) - gamma*k'*k(i)/twom;
        disp('Clustering ...iterated_genlouvain.m');
        tic        %%%%%%Clustering
        [SScommunity,QQ,n_it]=iterated_genlouvain(B,limit);
        toc
%         QQ = QQ/twom;
        clear B AdjMat;

        for iC=1:max(SScommunity)
            Com(:,:,iC)=reshape(nanmean(ParamOut.sample(:,sampleNeedIndex(SScommunity==iC)),2),Div(1),Div(2));
        end
        
figure;
for iC=1:max(SScommunity)             
                
            subplot(1,max(SScommunity),iC);
            imagesc(PhaseBin,ParamOut.Fplot, Com(:,:,iC)');axis xy;
            set(gca,'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'})
            if iC==1
               ylabel('Frequency Hz') 
            elseif iC==2
               xlabel('Theta phase rad'); 

            end
            
 end



