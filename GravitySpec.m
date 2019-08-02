function Output=GravitySpec(mapO,thRateRatio,PhasePlot,Fplot)

%%%%%%mapO: Spectrum map (Fre*Phase)
%%%%%%Gravity Threshold of Place Field was selected thRateRatio*PeakFiringRate
%%%%%%Place field was defined as >thNum of pixel with firing rate>thRateRatio*PeakFiringRate
PeakPhase=[];
PeakFre=[];
GravityPhase=[];
GravityFre=[];

if length(size(mapO))==2
    [PeakPhase,PeakFre,GravityPhase,GravityFre,GravNum,GravityPhaseStd,GravityFreStd]=GravitySpecMat(mapO,thRateRatio,PhasePlot,Fplot);   
elseif length(size(mapO))==3
    for j=1:size(mapO,3)
    [PeakPhase(j),PeakFre(j),GravityPhase(j),GravityFre(j),GravNum(j),GravityPhaseStd(j),GravityFreStd(j)]=GravitySpecMat(squeeze(mapO(:,:,j)),thRateRatio,PhasePlot,Fplot);   

    end
    elseif length(size(mapO))==4
    for j=1:size(mapO,3)
        for i=1:size(mapO,4)
            [PeakPhase(j,i),PeakFre(j,i),GravityPhase(j,i),GravityFre(j,i),GravNum(j,i),GravityPhaseStd(j,i),GravityFreStd(j,i)]=GravitySpecMat(squeeze(mapO(:,:,j,i)),thRateRatio,PhasePlot,Fplot);   
        end
    end


else

end
Output.PeakPhase=PeakPhase;
Output.PeakFre=PeakFre;
Output.GravityPhase=GravityPhase;
Output.GravityFre=GravityFre;
Output.GravNum=GravNum;
Output.GravityPhaseStd=GravityPhaseStd;
Output.GravityFreStd=GravityFreStd;



function [PeakPhase1,PeakFre1,GravityPhase1,GravityFre1,GravNum,GravityPhaseStd1,GravityFreStd1]=GravitySpecMat(map1,thRateRatio1,PhasePlot1,Fplot1)

   thNum=30;
   Cell=PlaceFieldFind2D(map1,thRateRatio1,thNum);
   GravX1=[];
   GravY1=[];
   for i=1:length(Cell)
       GravNum(i)=length(Cell(i).IndX); 
   end
   [GravNum,II1]=max(GravNum);
       GravX1=Cell(II1).IndX(:);
       GravY1=Cell(II1).IndY(:);

   
   PeakRate1=max(max(map1));
   [GravX1,GravY1,w1]=find(map1>=PeakRate1*thRateRatio1);
%    Weighted=map1(GravX1);
%    Weighted=Weighted-min(Weighted);
%    Weighted=Weighted/sum(Weighted);
   
   for i=1:length(GravX1)
       Weighted(i)=map1(GravX1(i),GravY1(i));
   end
   Weighted=Weighted-min(Weighted);
   Weighted=Weighted/sum(Weighted);

   GravNum=length(GravX1);
   [PeakX1,PeakY1]=find(map1>=PeakRate1);
   PeakPhase1=PhasePlot1(PeakX1);
   PeakFre1=Fplot1(PeakY1);
   
   kkk=Fplot1(GravY1);
   GravityFre1=sum(kkk(:).*Weighted(:));

   temp=PhasePlot1(GravX1);
   m=zeros(length(PhasePlot1),1);
   for i=1:length(PhasePlot1)
       kk=find(GravX1==i);
       if ~isempty(kk)
       m(i)=sum(Weighted(kk));
       end
   end
%    m=m.*Weighted;
%    GravityPhase1=circ_mean(temp(:));
    GravityPhase1 = circ_mean(PhasePlot1(:), m(:));
    GravityPhaseStd1 = circ_std(PhasePlot1(:), m(:));
    
    PopulationNum=100;
    SampleK=PopulationNum*Weighted(:);
    GravityFreStd1=sqrt(sum(SampleK(:).*(kkk(:)-GravityFre1).^2)/PopulationNum);
    



function Cell=PlaceFieldFind2D(mapO,thRateRatio,thNum)

%%%%%%map0: Firing Map
%%%%%%Firing Rate threshold of Place Field was selected thRateRatio*PeakFiringRate
%%%%%%Place field was defined as >thNum of pixel with firing rate>thRateRatio*PeakFiringRate
map=mapO;
PeakRate=max(max(map));

sigMap=zeros(size(map));

RateTH=PeakRate*thRateRatio;
sigMap(map>=RateTH)=1;

[m,n]=size(map);

sigCheck=1:n;
[index1,index2]=find(map==PeakRate);
index1=index1;
index2=index2;
StepMap=sigMap;

pixalM=1;
SigAll1=index1(:);
SigAll2=index2(:);
FieldNum=0;
Cell(1).IndX=[];
Cell(1).IndY=[];

while sum(sum(StepMap))>0
    for i=1:length(index1)
    StepMap(index1(i),index2(i))=0;
    end
    cal=[];
    cal1=[];
    cal2=[];
    tempMap=zeros(size(map));
    for i=1:length(index1)
        tcal1=[index1(i)-1 index1(i)+1 index1(i) index1(i)];
        tcal2=[index2(i) index2(i) index2(i)-1 index2(i)+1];
%           tcal1=[index1(i)-1 index1(i)-1 index1(i)-1 index1(i) index1(i) index1(i)+1 index1(i)+1 index1(i)+1];
%           tcal2=[index2(i)-1 index2(i) index2(i)+1 index2(i)-1 index2(i)+1 index2(i)-1 index2(i) index2(i)+1];

        invalid=union(find(tcal1<1|tcal1>m),find(tcal2<1|tcal2>n));
        
        tcal1(invalid)=[];
        tcal2(invalid)=[];
        cal1=[cal1;tcal1(:)];
        cal2=[cal2;tcal2(:)];
%         sigMap(cal1,cal2).*
    end
    if ~isempty(cal1)
    for i=1:length(cal1)
        tempMap(cal1(i),cal2(i))=1;
    end
    end
    
    for i=1:length(index1)
    tempMap(index1(i),index2(i))=0;
    end
    
    checkMap=tempMap.*StepMap;
    pixalM=pixalM+nansum(nansum(checkMap));

    [index1,index2]=find(checkMap==1);
    if ~isempty(index1)
        SigAll1=[SigAll1;index1(:)];
        SigAll2=[SigAll2;index2(:)];
        for i=1:length(index1)
            StepMap(index1(i),index2(i))=0;
        end
        if sum(sum(StepMap))==0
           if pixalM>=thNum
           FieldNum=FieldNum+1;
           Cell(FieldNum).IndX=SigAll1;
           Cell(FieldNum).IndY=SigAll2;
           end

        end
    else
        if pixalM>=thNum
           FieldNum=FieldNum+1;
           Cell(FieldNum).IndX=SigAll1;
           Cell(FieldNum).IndY=SigAll2;
        end
        map(find(StepMap==0))=0;
        PeakRate=max(max(map));
        [index1 index2]=find(map==PeakRate);
        pixalM=length(index1);
        SigAll1=index1;
        SigAll2=index2;
      

    end
        
    
end


for i=1:length(Cell)
    Rate=[];
    if ~isempty(Cell(i).IndX)
        for j=1:length(Cell(i).IndX)
            Rate(j)=mapO(Cell(i).IndX(j),Cell(i).IndY(j));
        end
        [Cell(i).PeakR,P]=max(Rate);
%         Cell(i).PeakI=[Cell(i).IndX(P) Cell(i).IndY(P)];
        Cell(i).PeakIX=Cell(i).IndX(P);
        Cell(i).PeakIY=Cell(i).IndY(P);
    else
        Cell(i).PeakR=[];
%         Cell(i).PeakI=[Cell(i).IndX(P) Cell(i).IndY(P)];
        Cell(i).PeakIX=[];
        Cell(i).PeakIY=[];
  

    end
end
