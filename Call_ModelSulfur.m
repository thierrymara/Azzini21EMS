%This program applies the IA-MC estimator to the Radiative Forcing Model 
%introduced in Tatang et. (1999) and like described in Azzini et al. (2021).
%References
%Tatang, M.A., Pan, W.W., Prin, R.G., McRae, G.J., 1997. 
%An efficient method for parametric uncertainty analysis of numerical geophisical model.
%J. Geophysics Research 102, 21925–21932.
% Azzini, I., Mara, T.A., Rosati, R., 2021.
%Comparison of two sets of Monte Carlo estimatiors of Sobol' indices
%Environ. Modell. and Softw. 144, 105167
warning('off')
pkg load statistics
struct_levels_to_print(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all
Nvar=9;%Number of Model Inputs;
%All the inputs of the model are log-normally distributed
for p=1:Nvar
  PDF.Type{p}='LogNormal';%'Gaussian';%
end
%Parameters of the Log-Normal distributions associated with the input model
PDF.Coeff=[0.76,1.2;0.39,1.1;0.85,1.1;0.3,1.3;5.0,1.4;1.7,1.2;71,1.15;0.5,1.5;5.5,1.5];
VarNames={'T';'Ac';'Rs';'Beta';'Psi';'f_psi';'Q';'Y';'L'};

Nsample=1000;
for rep=1:2%Do it twice for the A- and the B-sample resp.
  U=rand(Nsample,Nvar);%Note that in Azzini et al. 2021, Latin Hypercube Sampling was used
  for p=1:Nvar
    Sigma = log(PDF.Coeff(p,2));
    Mu = log(PDF.Coeff(p,1));
    XX(:,p)=exp(sqrt(2)*erfinv(2*U(:,p)-1)*Sigma + Mu);%LogNormal transformation
    if p==2
      XX(:,p) = 1-XX(:,p);%A_c
    elseif p==3
      XX(:,p) = 1-XX(:,p);%R_s
    end
  end
  if rep==1
    A = XX;
  else
    B = XX;
  end
end


yA=Sulfur_Model(A);
yB=Sulfur_Model(B);  
Dataset_A=[A,yA];
Dataset_B=[];
for i=1:Nvar,
  Ai=B;Ai(:,i)=A(:,i);
  Bi=A;Bi(:,i)=B(:,i);
  yAi(:,1)=Sulfur_Model(Ai);
  yBi=Sulfur_Model(Bi);
  Dataset_A=[Dataset_A;[Bi,yBi]];
  Dataset_B=[Dataset_B;[Ai,yAi]];
end
Dataset_B=[Dataset_B;[B,yB]];
Input = [Dataset_A(:,1:Nvar);Dataset_B(:,1:Nvar)];%This is how the dataset must be organised
Output = [Dataset_A(:,Nvar+1);Dataset_B(:,Nvar+1)];
SobolIndices = MC_IA(Input,Output,Nsample);
%First- and Total-order Sobol' indices (point estimate)
format short g
[SobolIndices.Si',SobolIndices.STi']
