clc
clear all


F = 1.0; %flow rate ml/min
Vc = 27.5;

Sf = 0.83; % stationary retention factor (0~1)
P = Sf/(1-Sf);% phase ratio 

Vinj = 1; %mL
Kd = [0.7   1.0   1.5   2.5];
C0 = 1;
Vcm = 100; %mL  elution volume

Ncup = 312;

%check number of species
k = length(Kd);

%Calculate volumes
Vcup = Vc/Ncup;  %Cup volume calc
Vs = Vc*Sf;
Vm = Vc*(1-Sf);
vmcup = Vcup*(1-Sf);

%Calculate timestep duration
dt_elution = vmcup/F; % s, this assumes that matrix is square

Turn_elution = round(Vcm/vmcup); %Number of timesteps that need to be run in order to hit switch volume (Vcm)

for i = 1:k
    
    KD = Kd(i);
%calculation for elution peaks
%This is stage 1, Classical Elution!
% [Nturn Cout Ccup] = CupV2(Sf, KD, Vc, Ncup, Turn_elution, C0, Vinj);  

[Nturn Cout Y Ccup] = CupV3(Sf, KD, Vc, Ncup, Turn_elution, C0, Vinj);  

%This data is returned in the form where the rows are cells, and the
%columns are timesteps

telute1  = (dt_elution).*Nturn; %elution time (min)
Velute1 = vmcup.*Nturn;


Vtot = Velute1;
Ttot = Vtot./F-Vinj/F/2;
%  figure()
 plot(Ttot , Cout, 'linewidth',2.0)
 set(gca,'FontWeight','bold','FontSize',14)
 title('Overall Elution Histories')
 xlabel('Elution Time (min)')
 ylabel('Concentration (mg/mL)')
 
hold on
%analytical solution 
dt = 0.01;
et = Vcm/F; %elution time segment
% tau = linspace(0,et,et/0.01); %dimensionless time
% t = tau*Vc/F;      %elution time, min
t = linspace(0,et,et/dt);
tau = t.*(F/Vc);
n = Ncup;
p = 1./(1-Sf+Sf.*KD);
Cinj = C0*Vinj/Vc; %input volume

% C0 = C0/1000;
a1 = sqrt((n-1)/(2*pi));  %approximate factorial(n-1) to avoid double expression exceeding problem
a2 =p.*(n./(n-1)).^n;
Cn2 = Cinj.*a1.*a2.*(p.*tau).^(n-1).*exp(n-1-n.*p.*tau); %analytical solution of ODE (n-series CSTR)
 

a = -(1-p.*tau).^2;
xn = Cinj.*p./sqrt(2*pi/n).*exp(a.*(n/2));   % Apply to gaussian peak distribution 

% plot(t,Cn2,t,xn)
plot(t , xn,'-.', 'linewidth',2.0)
 set(gca,'FontWeight','bold','FontSize',14)
hold on
trapz(t,xn)
X = [t'  Cn2'];
 

end



% Sim	0.3698	0.273	0.181	0.1171
% analy	0.3412	0.256	0.189	0.1141
% error	0.077339102	0.062271062	-0.044198895	0.025619129
