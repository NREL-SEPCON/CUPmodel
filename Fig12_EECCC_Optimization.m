clc
clear all

% QF = [1  2  3  4  5  6  7  8  9   10]; %flow rate ml/min
% F = 1.0;
% n =length(F);
Vinj = 1;

Vc = 27.5; % column volume, ml
C0 = [1   1]; % injection concentration (g/ml)- assume C0 in mobile phase
Vinj = 1.0; %injection volume, mL  - assume no slug volume
KD = [2.5  3.0]; %partition coefficient = C[SP]/C[MP]; KD order from lowest to highest

% VCM = 35:5:90;
% n =length(VCM);

VCM = 35:10:95; 
m = length(VCM);

Qf = 0.5:0.35:2; %flowrate variation
n = length(Qf);

% Vin = 0.5:0.1:1;  %injection volume variation
% m = length(Vin);


%define variables
Y1 = zeros(n,m);
Y2 = zeros(n,m);
P1 = zeros(n,m);
P2 = zeros(n,m);
Pr = zeros(n,m);
Rs = zeros(n,m);


for i = 1:n  %loop for Qf
for j = 1:m  %loop for VCM or Vinj injection volume
    
F = Qf(i);
% Vinj = Vin(j);
Vcm = VCM(j);

    Sf = 0.9821-0.1426*F;  %analytical scale correlation

    P = Sf/(1-Sf);% phase ratio 



    Ncup = 50*F^2-298.22*F+556.2;  %analy scale correlation
    Ncup = round(Ncup);

    Vcup = Vc/Ncup;
        Vs = Vc*Sf;
        Vm = Vc*(1-Sf);
        vmcup = Vcup*(1-Sf);
        Tturn = vmcup/F*60; % mobile phase turnover time (s)

    Ninj = Vinj/vmcup;

    dt_elution = vmcup/F; % this assumes that matrix is square

    Turn_elution = round(Vcm/vmcup); %Number of timesteps that need to be run in order to hit switch volume (Vcm)
    [Nturn Cout Ycm Xcm] = CupV3(Sf, KD, Vc, Ncup, Turn_elution, C0, Vinj);  

    [Xtot, Ytot, Ct_out, Xe, Ye] = EECCC_V6(KD, Vc, Sf, Xcm, Ycm);

    telute  = (dt_elution).*Nturn; %elution time (min)
    Velute = vmcup.*Nturn;


    Vtot = Velute;

    Vex = Ct_out(:,1); 
    Vext = Velute(end)+ Vex';
    Vtot = [ Velute     Vext];
    Ttot = Vtot./F;

    %COMBINE ELUTION DATA FROM ELUTION TO EECCC
    for im = 1:2
        C_EECCC(:,im) = Ct_out(:,im+1); %collect conc. outlet histories at EECCC
        
    end
        Ctotal = [ Cout    C_EECCC'];

    sig = [Ttot;  Ctotal]';

    %given a purity spec for both components, return the cutoff times
    %and the resulting yield, purity
    % lowcut = 0.005; %indicates that peak will not be collected until signal is lowcut*maxpeakheight and collection will begin once signal is below lowcut*maxpeakheight

    [row col] = size(Ctotal);
     
    km = min (max(Ctotal(1,:)),max(Ctotal(2,:)));
    theta = 0.01;  % lowcut percentage from the lowest peaks e.g. 0.01 means 1% of peak height
    lowcut = km*theta;
    purity = 0.99; %the returned values are for comp 1 in the first peak eluted and comp 2 in the second peak

    [tbreak,yield,purity,res] = CalArea2(sig,lowcut,purity);
   

    Y1(i,j) = yield(1)*100;
    Y2(i,j) = yield(2)*100;
    P1(i,j) = purity(1)*100;
    P2(i,j) = purity(2)*100;
    Rs(i,j) = res;

  
    
    
    tEnd = tbreak(4)/60; % process time, hour

    pr = Vinj.*C0.*yield./(Vc/1000)./tEnd;  %productivity = input*purity*yield/column volume/process time

    Pr(i,j) = pr(1)+pr(2); %calculate total productivity
    
     if Pr(i,j) >= max(max(Pr))
        peak_F = i;
        peak_Vcm = j;
        
    end
    
    end
 clear C_EECCC  %need to clear cash 

end



FF = Qf./Vc;
VV = VCM./Vc;


figure (1)
surf(FF,VV,Y2')
xlabel('F/V_c (min^-^1)')
ylabel('V_C_M/V_c')
zlabel('Yield_2(%)')


figure (2)
surf(FF,VV,Pr')
xlabel('F/V_c (min^-^1)')
ylabel('V_C_M/V_c')
zlabel('Total Productivity(mg/L hr)')


figure (3)
contourf(VV,FF,Y2)
ylabel('F/V_c (min^-^1)')
xlabel('V_C_M/V_c')
zlabel('Yield_2(%)')
grid on

figure (4)
contourf(VV,FF,Pr)
ylabel('F/V_c (min^-^1)')
xlabel('V_C_M/V_c')
zlabel('Total Productivity(mg/L hr)')
grid on


