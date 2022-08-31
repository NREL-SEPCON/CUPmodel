% Dual mode CCC by reversing the stationary while keep mobile phase
% written by Hoon Choi
% Loading New SP into a column with reverse flow

%Xtot: mobile phase concentration from classic elution
%Ytot: stationary phase concentration from classic elution
%Col_out: Column outlet proflies for MP and SP respectively
%X: Concentration in MP during dual mode
%Y: Concentration in SP during dual mode


function [Xtot, Ytot, Col_out, X, Y] = Dual_V1(KD, Vc, Sf, Ccup, CcupY, Tau)

P = Sf/(1-Sf);
Vm = Vc*(1-Sf);
Vs = Vc*Sf;
kA = 1./(1+P.*KD); %Retention Factor!



[Ncup el_time comp] = size(Ccup);   %Ccup profile (Ncup,time,species)
X_elution = zeros(Ncup,  el_time,  comp); %Initial MP Concentration g/L   (X , t, j)
Y_elution = zeros(Ncup, el_time, comp);

% tspan = fix(ceil(Ncup*(1+Sf)*1.1)); % timespan for MP (Ncup) + SP (Ncup*Sf);1.1 to give enough elution time
Y = zeros(Ncup, Tau,  comp); %SP Concentration g/L
X = zeros(Ncup, Tau, comp); % MP concentration

%take previous profiles from classic elution to set initial boundary
Xin = zeros(Ncup,1,comp);
Yin = zeros(Ncup,1,comp);


for comp = 1:length(KD)
    X_elution(:,:,comp) = Ccup(:,:,comp); % Keep use the same format (i,t,j)
    Y_elution(:,:,comp) =  CcupY(:,:,comp);               %KD(comp).*X_elution(:,:,comp);
    Xin(:,1,comp) = X_elution(:,el_time,comp);
    Yin(:,1,comp) = Y_elution(:,el_time,comp);       %KD(comp).*Xin(:,1,comp);
    
end

Vcup = Vc/Ncup;
Vmcup = Vm/Ncup;
Vscup = Vs/Ncup;

n = length(KD);

for j = 1:n  %comp
    

    t=1; %initial condition at t = 1
    
    X(Ncup,1,j) = kA(j)*Xin(Ncup,1,j);
    Y(Ncup,1,j) = KD(j)*X(Ncup,1,j);

    for i = Ncup:-1:2  %calculate from Ncup to 1st cell
        
        X(i-1,t,j) = kA(j)*(Xin(i-1,1,j)+P*Yin(i,1,j));
        Y(i-1,t,j) = KD(j)*X(i-1,t,j);
               
    end
    
              
     for t = 2:Tau % SP reverse elution (dT = Vscup) 
                
         
         % Boundary condition first for Ncup cell
         X(Ncup,t,j) = kA(j)*X(Ncup,t-1,j);
         Y(Ncup,t,j) = KD(j)*X(Ncup,t,j);
         
         %SP reverse elution from [Ncup-1] to 1st cell
         for i = Ncup:-1:2
                
             X(i-1,t,j) = kA(j)*(X(i-1,t-1,j)+P*Y(i,t-1,j));
             Y(i-1,t,j) = KD(j)*X(i-1,t,j);

         end
         
     end

     %Combine concentration profiles into Xtot & Ytot from CM to EECCC
          Xtot = [X_elution  X];
          Ytot = [Y_elution  Y];
     

          %make elution time and volume matrix for EECCC 
          tspan = linspace(1, Tau, Tau);
          Vspan = Vscup.*tspan;
          
              
          
          for i = 1:n
              
              M(i,:) = Y(1,:,i);   % column outlet is 1st cell in dual mode
              M(i+n,:)= X(1,:,i);  % Column outlet is 1st cell
                                         
          end
          
          Col_out = [Vspan; M]';  % [Vol yi  xi] column outlet profiles
                


end

          
     
  