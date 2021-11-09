%% CUP MODEL CODE WAS WRITTEN BY HOON CHOI, NATIONAL RENEWABLE ENERGY LABORATORY 
% 
% CUP MODEL CODING for ELUTION-EXTRUSION COUNTER-CURRENT CHROMATOGRAPHY (EECCC)
% CUP MODEL DERIVATION IS DESCRIBED IN THE PAPER, CHOI ET AL. SEP & PUR TECH(2022)
% THIS CODE ILLUSTRATES FOR AN EECCC EXAMPLE OF THE FIGURE 7C IN THE PAPER
% 
% ANY QUESTIONS TO AUTHOR: hoon.choi@nrel.gov
%
%The research was sponsored by the U.S. Department of Energy (DOE), Energy Efficiency and Renewable
%Energy Office, Bioenergy Technologies Office (BETO) under the BETO Bioprocessing Separations Consortium via Contract No.
%DE-AC36-08GO28308 with the National Renewable Energy Laboratory and via Contract No. DE-AC05?76RL01830 with Pacific Northwest National Laboratory.

%% 


clc
clear all


F = 1.0; %flow rate, ml/min
Vc = 27.5;  %column volume, ml


Sf = 0.75; %stationary phase retention factor; Sf = V[SP]/Vc

P = Sf/(1-Sf);% phase ratio 

Vinj = 1; % feed injection volume, mL
Vd = 3;  %extracolumn volume [reservoir-to-injection loop]

%Distribution Coefficient or Parition coefficient KD = C[SP]/C[MP]
KD = [ 0.35  0.6  0.92  1.28  2.1];
C0 = [5    5      5      10      5]; %feed concentration, g/L


Vcm = 24+Vd; %Vcm = end of classic elution, ml

Ncup = 724;  %column efficiency for all compounds 

%check number of species
n = length(KD);

%Calculate volumes
Vcup = Vc/Ncup;  %Cell volume calc
Vs = Vc*Sf; %V[SP]
Vm = Vc*(1-Sf); %V[MP]
vmcup = Vcup*(1-Sf); %cell MP volume

%Calculate timestep duration
dt_elution = vmcup/F; %this assumes that matrix is square

Turn_elution = round(Vcm/vmcup); %Number of timesteps that need to be run in order to hit switch volume (Vcm)


%calculation for elution peaks
%This is stage for Classical Elution

[Nturn Cout Ycm Xcm] = CupV3(Sf, KD, Vc, Ncup, Turn_elution, C0, Vinj);  

%Start Elution-extrusion CCC after Vcm
[Xtot, Ytot, Ct_out, Xe, Ye] = EECCC_V6(KD, Vc, Sf, Xcm, Ycm);


%This data is returned in the form where the rows are cells, and the
%columns are timesteps

telute1  = (dt_elution).*Nturn; %elution time for CM, min
Velute1 = vmcup.*Nturn;  %elution volume for CM, ml

Vex = Ct_out(:,1); 
Vext = Velute1(end)+ Vex';
Vtot = [ Velute1     Vext]; %overall elution volume, ml
Ttot = Vtot./F; %overall elution time, ml

%COMBINE ELUTION DATA FROM CM TO EECCC
for i = 1:n
    C_EECCC(:,i) = Ct_out(:,i+1); 
    Ctotal(i,:) = [Cout(i,:)  C_EECCC(:,i)'];
    
end


%find retention volume 

 export = [Ttot; Vtot; Ctotal]'; %to export the data to excel

 
 plot(Vtot , Ctotal, 'linewidth',2.0)
 set(gca,'FontWeight','bold','FontSize',14)
 title('Overall Elution Histories')
 xlabel('Elution Volume (mL)')
 ylabel('Concentration (mg/mL)')

 Velute1(end);
 Velute1(end)+Vex(Ncup);
 
 %use the code for if MATLAB version can read 'xline' function
%  xline(Velute1(end)/F,'-.r'); %extrusion boundary
%  xline((Velute1(end)+Vex(Ncup))/F,'--'); %sweep boundary

 
 figure()
 plot(Ttot , Ctotal, 'linewidth',2.0)
 set(gca,'FontWeight','bold','FontSize',14)
 title('Overall Elution Histories')
 xlabel('Elution Time (min)')
 ylabel('Concentration (mg/mL)')
 
 
