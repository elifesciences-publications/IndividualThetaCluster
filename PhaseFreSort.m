function SortI=PhaseFreSort(Fre,Phase)

[FreS,SortI]=sort(Fre);

%%%%%%%
%%%%%%%currently used 8/8-2018  Get low,median gamma by fre, then get early
%%%%%%%and late gamma based on phase refers to median gamma
Sort1=SortI(1:2);
Sort2=setdiff(1:length(Fre),Sort1);
Sort3=FindEarlyPhase(Phase(Sort2),Phase(Sort1(2)));
%%%%%%%currently used 8/8-2018
%%%%%%%



% if FreS(2)<=120||abs(FreS(2)-FreS(4))<40
%    Sort1=SortI(1:2);
% else
%    Sort1=SortI(1);
% end
% % if FreS(2)<=80
% else
%    Sort1=SortI(1);
% end

% % if abs(FreS(2)-FreS(4))<40
% % %    Sort1=SortI(1:2);
% % % else
% %    Sort1=SortI(1);
% % else
% %    Sort1=SortI(1:2);
% % end
% % Sort1=SortI(1:2);
% [~,Sort3]=sort(Phase(Sort2));

SortI=[Sort1 Sort2(Sort3)];

function SortI=FindEarlyPhase(Phase,PhaseRef)

I=find(Phase<PhaseRef);

if ~isempty(I)
   Phase(I)=Phase(I)+2*pi; 
end
   [~,SortI]=sort(Phase,'descend');

