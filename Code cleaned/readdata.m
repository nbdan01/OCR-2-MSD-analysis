function [PAR,ORT]=readdata(loc)
%This function reads the PAR-ORT data from the excel file. The input is the
%a string with the excel location, the output are cells with the parallel
%and orthogonal data points
[~,~,R] = xlsread(loc); %Read the excel file. 
sr=size(R); 
col=1:sr(2); %tell how many collums there are

for i=col
T=R(:,i);
clear('Pc','Pnc')
Pc=zeros(1,length(T));
Pnc=zeros(1,length(T));
for j=1:length(T)
   Pc(j)=~isempty(strfind(T{j}, 'Para')); %The para data must start 
   %directly under a cell with 'para' or 'Para'. 
   Pnc(j)=(~isempty(strfind(T{j}, 'para')));
end
if ~isempty(find(Pc+Pnc,1))
   LP(i)=find(Pc+Pnc); %Get the position where the data starts
end

clear('Pc','Pnc') %The same, but then for orto
Oc=zeros(1,length(T));
Onc=zeros(1,length(T));
for j=1:length(T)
   Oc(j)=~isempty(strfind(T{j}, 'Orth'));
   Onc(j)=(~isempty(strfind(T{j}, 'orth')));
end
if ~isempty(find(Oc+Onc, 1))
   LO(i)=find(Oc+Onc);
end


clear('Oc','Onc') %The same, but then for orto
Tc=zeros(1,length(T));
Tnc=zeros(1,length(T));
for j=1:length(T)
   Tc(j)=~isempty(strfind(T{j}, 'frames'));
   Tnc(j)=(~isempty(strfind(T{j}, 'frames')));
end
if ~isempty(find(Tc+Tnc, 1))
%    TO(i)=find(Tc+Tnc);
end
end

O=find(LO); %All columns with orto data
P=find(LP); %All columns with para data
% TT=find(TO); %All columns with para data

for i=1:length(O)
C=O(i);
L=R(:,C);
sL=size(L);
nu=sL(1);
while max(isnan(R{nu,C}))==1
    nu=nu-1; %Get the size of the track (delete all the NaN's at the end)
end
if ~isempty(L(LO(C)+1:nu))
ORT{i}=cell2mat(L(LO(C)+1:nu)); %Make a cell with all the ortogonal trajectories
end
end

% for i=1:length(O)
% C=TT(i);
% L=R(:,C);
% sL=size(L);
% nu=sL(1);
% while max(isnan(R{nu,C}))==1
%     nu=nu-1; %Get the size of the track (delete all the NaN's at the end)
% end
% if ~isempty(L(TO(C)+1:nu))
% TIM{i}=cell2mat(L(TO(C)+1:nu)); %Make a cell with all the time trajectories
% end
% end

for i=1:length(P)
C=P(i);
L=R(:,C);
sL=size(L);
nu=sL(1);
while max(isnan(R{nu,C}))==1
    nu=nu-1;
end
if ~isempty(L(LP(C)+1:nu))
PAR{i}=cell2mat(L(LP(C)+1:nu)); %Make a cell with all the para trajectories
end
end

end
