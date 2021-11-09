%% CUP MODEL CODE WAS WRITTEN BY HOON CHOI, NATIONAL RENEWABLE ENERGY LABORATORY 
% 
% CUP MODEL CODING for Dual mode of COUNTER-CURRENT CHROMATOGRAPHY (DualCCC)
% CUP MODEL DERIVATION IS DESCRIBED IN THE PAPER, CHOI ET AL. SEP & PUR TECH(2022)
% THIS CODE ILLUSTRATES FOR AN EECCC EXAMPLE OF THE FIGURE 7D IN THE PAPER
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
Vd = 3;  %extracolumn volume [flying leads]

%Distribution Coefficient or Parition coefficient KD = C[SP]/C[MP]
KD = [ 0.35  0.6  0.92  1.28  2.1];
C0 = [5    5      5      10      5]; %feed concentration, g/L


Vcm = 27-Vd; %Vcm = end of classic elution, ml

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
%This is stage for Classical Elution.
[Nturn Cout Ycm Xcm] = CupV3(Sf, KD, Vc, Ncup, Turn_elution, C0, Vinj);  

%input new variables for DUAL mode
F_dual = 1;
Sf_dual = 0.55;
time_dual = 50;
vscup_dual = Vc*(Sf_dual)/Ncup;

dt_dual = vscup_dual/F_dual;
Vol_dual = time_dual*F_dual;
Turn_dual = round(Vol_dual/vscup_dual);

%Calculate DUAL MODE
[Xtot, Ytot, Ct_out, X, Y] = Dual_V1(KD, Vc, Sf_dual, Xcm, Ycm, Turn_dual);


%COMBINE ELUTION DATA FROM ELUTION TO EECCC
for i = 1:n
    C_Dual(:,i) = Ct_out(:,i+1); 
    Ctotal(i,:) = [Cout(i,:)  C_Dual(:,i)'];
  end

Vdual = Ct_out(:,1)';
tdual = Vdual./F_dual;
telute1  = (dt_elution).*Nturn; %elution time (min)
Velute1 = vmcup.*Nturn;
Vext = Velute1(end)+ Vdual;

%total elution volume and time
Vtot = [ Velute1     Vext]; %ml
Ttot = Vtot./F; %min

 export = [Ttot; Vtot; Ctotal]'; %table for data export

 plot(Vtot , Ctotal, 'linewidth',2.0)
 set(gca,'FontWeight','bold','FontSize',14)
 title('Overall Elution Histories')
 xlabel('Elution Volume (mL)')
 ylabel('Concentration (mg/mL)')

 
 
 %use the code for if MATLAB version can read 'xline' function
%  Velute1(end);
%  Velute1(end)+Vex(Ncup);
%  xline(Velute1(end)/F,'-.r'); %extrusion boundary
%  xline((Velute1(end)+Vex(Ncup))/F,'--'); %sweep boundary

 
 figure()
 plot(Ttot , Ctotal, 'linewidth',2.0)
 set(gca,'FontWeight','bold','FontSize',14)
 title('Overall Elution Histories')
 xlabel('Elution Time (min)')
 ylabel('Concentration (mg/mL)')
 