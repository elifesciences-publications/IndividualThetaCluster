function ParamOut = SpecThetaExtract(LFP,LFPtheta,thetaphase,samprate,timerange,WaveParam,PhaseStep)
%   SpecThetaWaveLet get wavelet Power in time&frequency domain representation
%   Extract Wavelet Spectrum for Each Single Theta Circle and Normalize
%   Spectrum by Theta Phase.
%   Coded by Lu Zhang; math2437@hotmail.com
%   Last updated: (19th March 2019)

%Output:
%ParamOut.sample is the vectorized spectrum normalized in theta phase
%ParamOut.ThetaTS(i) is timeStamps(theta trough) of sample(:,i)
%ParamOut.Div dimension parameter;
%ParamOut.Fplot Frequency band for the spectrum data;

%Input
%LFP: raw LFPs
%LFPtheta: filtered LFPs in theta Band
%thetaphase: phase of filtered LFPs, usually calculated by hilbert transform;
%samprate: sampling rate;
%timerange: analysis period.[1 40 120;30 88 160] indicate 1-30s,40-88s and 120-160s for three analysis period in total;

% % Default setting
%WaveParam.ntw=100;      smoothing time-window 100
%WaveParam.nsw=3;        smoothing fre-window 3
%WaveParam.smoothW=10;   smoothing Power
%WaveParam.DownSample=2; DownSampling
%WaveParam.wname='morl'; Wavelet Name
%WaveParam.samprateD=samprateD;
%WaveParam.Freq=[1:150]; Frequeny band interested
%PhaseStep: Phase bin size;


PhaseBin=-pi:PhaseStep:pi;

ThetaSpec={};
ThetaPhase={};
LFPOutput={};
ThetaTimeStamps=[];


if isempty(timerange)
    sample=[];
    ThetaTS=[];
    F=WaveParam.Freq;
    Fplot=F(end:-1:1);
    Div=[];
    return
end

tic
Duration=round(diff(timerange));
NumThetaPeriod=length(Duration);

Istart=round(timerange(1,:)*samprate)+1;
Iend=round(timerange(2,:)*samprate)+1;
timerangeI=round(timerange*samprate)+1;

%%%%%%%%ii loops      multiple theta period in one recording file
for ii=1:length(Duration)

    %%%%%%%%Inlcuding 1s before the theta start and 1s after theta ends
    %%%%%%%%to avoid edge effects in time of Wavelet
    Istart(ii)=max(Istart(ii)-samprate,1);
    Iend(ii)=min(Iend(ii)+samprate,length(LFPtheta(:,1)));
    tempI=Istart(ii):Iend(ii);
    if isempty(tempI)
        continue
    end
    
    TExclude=[timerangeI(1,ii)-Istart(ii) Iend(ii)-timerangeI(2,ii)]/samprate;
    TimeInclude(:,ii)=[Istart(ii);Iend(ii)]/samprate;
    TInclude=([Istart(ii) Iend(ii)]-Istart(ii))/samprate;
    TInclude=[TInclude(1)+TExclude(1) TInclude(2)-TExclude(2)];
    
    LFPtemp=LFP(tempI);
    LFPthetatemp=LFPtheta(tempI);
    LFPphasetemp=thetaphase(tempI);

    LFPtemp=zscore(decimate(LFPtemp,WaveParam.DownSample));
    LFPthetatemp=downsample(LFPthetatemp,WaveParam.DownSample);
    LFPphasetemp=downsample(LFPphasetemp,WaveParam.DownSample);

    samprateD=samprate/WaveParam.DownSample;
    
F=WaveParam.Freq;
wname=WaveParam.wname;
ntw=WaveParam.ntw;
nsw=WaveParam.nsw;
bin_width=1/samprateD;
fc = centfrq(wname);
scales=sort(fc./F.*samprateD);

flag_SMOOTH=true;
cfs_s1    = cwt(LFPtemp,scales,wname);
cfs_s10   = cfs_s1;
cfs_s1    = smoothCFS(abs(cfs_s1).^2,flag_SMOOTH,nsw,ntw);
S1Power    = sqrt(cfs_s1);



if WaveParam.Zscore
   S1Power=zscore(S1Power,0,2); %%%%%%%%%%%%%normalize over time
end

Time=([1:length(LFPtemp)]-1)/samprateD;
FPlot=F(end:-1:1);

temp=LFPphasetemp;
len_t=length(temp);

%%%%%%%detect theta trough
         thetamax_ts=[];
         for j=1:len_t
             if j==1&&temp(j)<-3.141
                 thetamax_ts=[thetamax_ts,j];
             elseif j==len_t&&temp(j)<-3.141
                 thetamax_ts=[thetamax_ts,j];
             elseif j~=1&&j~=len_t
                 if (temp(j-1)-temp(j))>6.0&&(temp(j+1)-temp(j))>0
                     thetamax_ts=[thetamax_ts,j-1+(pi-temp(j-1))/((pi-temp(j-1))+(pi+temp(j)))];

                 end
             else
                 
             end
         end
thetamax_ts=(thetamax_ts-1)/samprateD;
thetamax_ts(thetamax_ts>TInclude(2))=[];
thetamax_ts(thetamax_ts<TInclude(1))=[];
%%%%%%%detect theta trough

thetamax_ind=round(thetamax_ts*samprateD);
% thetamax_ts=round(thetamax_ts);
% % figure;
% % % plot(LFPthetatemp(:,2));hold on;plot(thetamax_ind,-pi,'r.')
% % plot(LFPthetatemp(:,1));hold on;plot(thetamax_ind,0,'r.')

% BackW=round(WaveParam.Range*samprateD);
% ForW=BackW;

if ii==1
    tempPhase=LFPphasetemp(:);
else
    tempPhase=[tempPhase;LFPphasetemp(:)];
end

thetamax_ts=thetamax_ts+TimeInclude(1,ii);
thetamax_ts=thetamax_ts(:);
for iii=1:(length(thetamax_ind)-1)
    s=round(max(thetamax_ind(iii),1));
    o=round(min(thetamax_ind(iii+1),len_t));
    ThetaSpec{end+1,1}=S1Power(:,s:o)';
    ThetaPhase{end+1,1}=LFPphasetemp(s:o);
    LFPOutput{end+1,1}=LFPtemp(s:o);
end
ThetaTimeStamps=[ThetaTimeStamps;thetamax_ts(1:(end-1))];


end
%%%%%%%%ii loops      multiple theta period in one recording file
   
   ThetaSpecN = aveY_discretizeX(ThetaSpec,ThetaPhase,PhaseBin);    
   for i=1:length(ThetaSpecN)
       ThetaSpecNew(:,:,i)=ThetaSpecN{i};
   end
   clear ThetaSpec;
   
   ThetaTS=ThetaTimeStamps; clear ThetaTimeStamps;
   
%     WaveParam.Freq=[20:2:120];
    Fplot=WaveParam.Freq(end:-1:1);
    Div=size(ThetaSpecNew);
    sample=reshape(ThetaSpecNew,Div(1)*Div(2),Div(3));
    Ivalid=find(sum(isnan(sample),1)>0);
    sample(:,Ivalid)=[];
    ThetaTS(Ivalid)=[];
    LFPOutput(Ivalid)=[];
    
    
for i=1:size(sample,2)
    sample(:,i)=sample(:,i)-min(sample(:,i));
    sample(:,i)=sample(:,i)/sum(sample(:,i))+0.00001;
    sample(:,i)=sample(:,i)/sum(sample(:,i));
end

ParamOut.sample=sample;
ParamOut.ThetaTS=ThetaTS;
ParamOut.LFPOutput=LFPOutput;
ParamOut.Fplot=Fplot;
ParamOut.Div=Div;
ParamOut.Div=Div;

toc

%----------------------------------------------------------------------
function CFS = smoothCFS(CFS,flag_SMOOTH,NSW,NTW)

if ~flag_SMOOTH , return; end
if ~isempty(NTW)
    len = NTW;
    F   = ones(1,len)/len;
    CFS = conv2(CFS,F,'same');
end
if ~isempty(NSW)
    len = NSW;
    F   = ones(1,len)/len;    
    CFS = conv2(CFS,F','same');
end

%----------------------------------------------------------------------
function yy = aveY_discretizeX(y,x,xbin,varargin)

%%%%%%%y is m*n matrix, x is m*1 vector
%%%%%%%divided x into xbin
%%%%%%%averaged the y in each x xbin.

if ~iscell(x)
    bins = discretize(x,xbin);    
    for ii=1:(length(xbin)-1)
        yy(ii,:)=mean(y(bins==ii,:),1);
    end
elseif iscell(x)
    yyy={};
    for i=1:size(x,1)
        for j=1:size(x,2)
            yyy{i,j}=aveY_discretizeX(y{i,j},x{i,j},xbin);
        end 
    end
    yy=yyy;
    clear yyy
elseif isempty(x)
    yy=[];  
else
    yy=[];  
  
end











