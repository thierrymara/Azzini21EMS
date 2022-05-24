%% Copyright (C) 2020 Thierry A. MARA
%% 
%% This program is free software: you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation, either version 3 of the License, or
%% (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program.  If not, see
%% <https://www.gnu.org/licenses/>.

%% -*- texinfo -*- 
%% Usage: Out = MC_IA (X, y, SampleSize)
%% Compute the First- and Total-Order Sobol' indices and their bootstrap 95% confidence intervals
%% with the IA estimator (see Azzini I., Mara T., Rosati R., 2021).
%%
%% X is the input sample following the IA design of size Ntxd
%%
%% y is the associated vector of model responses
%%
%% d the number of inputs
%%
%% Nt = 2*Ns*(Ngr+1), with Ngr = the Number of groups (Ngr<=d)
%%
%% Ns is the sample size of each block [Xa;Xa_1;...Xa_Ngr;Xb_1;...,Xb_Ngr;Xb]
%%
%% Out.Si is a Ngrx3 matrix of  first-order Sobol' indices of each group
%%
%% Out.STi is a Ngrx3 matrix of  total-order Sobol' indices of each group
%%
%% Out.S(T)i = [95% bootstrap lower bound, Pointwise estimate, 95% bootstrap upper bound]
%%
%% Author: marathi <marathi@D01RI1602271>
%% Created: 2020-02-27

function FirstTotalOrder_SI = MC_IA (X, y, Nsample)
[Ns,Nvar]=size(X);
if (floor(Ns/Nsample)==ceil(Ns/Nsample))
else
  FirstTotalOrder_SI = [];
  disp("The input sample and the sample size are not consistent with the IA design!")
  return
end
From=Ns/2;
Gr=1:Nvar;
%%% Check Number of Groups %%%%%%%%%%
XA=X(1:Nsample,:);%The A
yA=y(1:Nsample);%y_A
XB=X((end-Nsample+1):end,:);%The B
yB=y((end-Nsample+1):end);%y_B
cnt=0;
for p=1:Nvar
  if Gr(p)==p
    cnt=cnt+1;
    Xp=X(cnt*Nsample+1:(cnt+1)*Nsample,:);
    DeltaX=sum(abs(XB-Xp));
    for pp=p:Nvar
      if DeltaX(pp)==0%True
        Gr(pp)=cnt;
      endif
    end
  endif
endfor
NumGroup=cnt;
fprintf('The number of groups is: %d\n',NumGroup);
%Pointwise estimate
cnt=0;
for Ngr=1:NumGroup
    cnt=cnt+1;
    InpNum=find(Gr==cnt);
    yBi=y(Ngr*Nsample+1:(Ngr+1)*Nsample);
    yAi=y(From+(Ngr-1)*Nsample+1:From+Ngr*Nsample);
    %First-order
    Si1 = mean((yA-yBi).*(yAi-yB));%
    %Total-order
    STi1=mean((yA-yAi).*(yBi-yB));%
    %%%%%%%%%%%%%%%%%%%%%%%
    Subscript{Ngr}=num2str(InpNum(1));
    for i=2:length(InpNum)
      Subscript{Ngr}=[Subscript{Ngr},'_',num2str(InpNum(i))];
    endfor
    Si(1,Ngr)=Si1/(mean((yA-yB).^2+(yAi-yBi).^2)/2);%The pointwise-estimate seems more accurate with Si2
    STi(1,Ngr)=1-STi1/(mean((yA-yB).^2+(yAi-yBi).^2)/2);
endfor
FirstTotalOrder_SI.Si=Si;
FirstTotalOrder_SI.STi=STi;
%Bootstrap estimate
NBoot=100;
for B=1:NBoot;
  BootSample=ceil(Nsample*rand(Nsample,1));
  yAA=yA(BootSample);%y_A
  yBB=yB(BootSample);%y_B
  for GrN=1:NumGroup
    yBi=y(GrN*Nsample+1:(GrN+1)*Nsample);
    yAi=y(From+(GrN-1)*Nsample+1:From+GrN*Nsample);
    yBi=yBi(BootSample);
    yAi=yAi(BootSample);
    %First-order
    Si1 = mean((yAA-yBi).*(yAi-yBB));
    %Total-order
    STi1=mean((yAA-yAi).*(yBi-yBB));
    %%%%%%%%%%%%%%%%%%%%%%%
    Si_Boot(B,GrN)=Si1/(mean((yAA-yBB).^2+(yAi-yBi).^2)/2);%
    STi_Boot(B,GrN)=1-STi1/(mean((yAA-yBB).^2+(yAi-yBi).^2)/2);%
  endfor
endfor
%Estimate
for GrN=1:size(Si_Boot,2)
  [Si_min,Si_med,Si_max]=Uncertainty(Si_Boot(:,GrN),0.05);
  [STi_min,STi_med,STi_max]=Uncertainty(STi_Boot(:,GrN),0.05);
  eval(['FirstTotalOrder_SI.S' Subscript{GrN} '= [Si_min;Si_med;Si_max];']);
  eval(['FirstTotalOrder_SI.ST' Subscript{GrN} '= [STi_min;STi_med;STi_max];']);
end
endfunction
