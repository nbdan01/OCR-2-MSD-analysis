function checkparort(PAR,ORT)
%This function checks the sizes of PAR and ORT. If these are not of the
%same size, something has gone wrong. Most likely, a coumn was
%unintentionally labeled as para or ort.
sPAR=size(PAR);
sORT=size(ORT);
if min(sPAR==sORT)<1
        error('The number of columns with Para is unequal to Orto. There may be a Para or Orto on a location that is not correct.')
end
for i=1:length(PAR)
   p=PAR{i};
   o=ORT{i};
   if length(p)~=length(o)
       error('For some trajectory, the number of timepoints of par is unequal to the number of timepoints of ort.')
   end
end
end
