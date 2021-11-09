% EECCC - extrusion mode of CCC by replacing mobile phase with SP
% written by Hoon Choi

% Assume the same flowrate is used
% Turnover time is automatically calculated 
% sweep state: turnover = Ncup, dT = vmcup 
% extrution state: turnover = Ncup*Sf, dT = vcup
%Xtot: total MP profiles (CM+EECCC)
%Ytot: total SP profiles (CM+EECCC)
%Ct_out: elution volume vs. column outlet concentration
%X: MP profiles during EECCC
%Y: SP profiles during EECCC
% Xcm: MP profiles from CM
% Ycm: SP profiles from CM


function [Xtot, Ytot, Ct_out, X, Y] = EECCC_V6(KD, Vc, Sf, Xcm, Ycm)

P = Sf/(1-Sf); %phase ratio
Vm = Vc*(1-Sf); % MP volume
% Vs = Vc*Sf;
kA = 1./(1+P.*KD); %Retention Factor!



[Ncup el_time comp] = size(Xcm);   %Ccup profile (Ncup,time,species)
X_elution = zeros(Ncup,  el_time,  comp); %Initial MP Concentration g/L   (column , time , species)
Y_elution = zeros(Ncup, el_time, comp);

tspan = fix(ceil(Ncup*(1+Sf)*1.3)); % timespan for MP (Ncup) + SP (Ncup*Sf);x1.3 to give enough elution time
Y = zeros(Ncup, tspan,  comp); %SP Concentration g/L
X = zeros(Ncup, tspan, comp); % MP concentration

%load previous concentration profiles from classic elution to set initial boundary
Xin = zeros(Ncup,1,comp);
Yin = zeros(Ncup,1,comp);


for comp = 1:length(KD)
    X_elution(:,:,comp) = Xcm(:,:,comp); % Keep use the same format (i,t,j)
    Y_elution(:,:,comp) =  Ycm(:,:,comp);               %KD(comp).*X_elution(:,:,comp);
    Xin(:,1,comp) = X_elution(:,el_time,comp);
    Yin(:,1,comp) = Y_elution(:,el_time,comp);       %KD(comp).*Xin(:,1,comp);
    
end

Vcup = Vc/Ncup;
Vmcup = Vm/Ncup;

n = length(KD);

for j = 1:n  %comp
    
%EECCC start - sweep state

    t=1; %initial condition at t = 1
    Y(1,1,j) = Sf*Yin(1,1,j);

    for i = Ncup:-1:2
        
        X(i,t,j) = kA(j)*(Xin(i-1,1,j)+P*Yin(i,1,j));
        Y(i,t,j) = KD(j)*X(i,t,j);
               
    end

              
     for t = 2:Ncup % MP elution end at Ncup time 
                
          for i = 2:Ncup      
           
                if i <= t
% one MP cell is displaced by SP after one turnover                    

                    Y(i,t,j) = (1-Sf)*Y(i-1,t-1,j)+Sf*Y(i,t-1,j);

                 elseif i > t 
                    X(i,t,j) = kA(j)*(X(i-1,t-1,j)+P*Y(i,t-1,j));
                    Y(i,t,j) = KD(j)*X(i,t,j);
                 end
                 
          end
          
         
     end
     
     %Sweep end  where dT = (1-Sf)Vc/Ncup = vmcup
     %%Extrusion start  where dT = Vc/Ncup =vcup
          
     for t = Ncup+1:tspan
     
         for i = 2:Ncup
                
             Y(i,t,j) = Y(i-1,t-1,j); %move forward to outlet without mixing

         end
         
     end

     %Combine concentration profiles into Xtot & Ytot from CM to EECCC
          Xtot = [X_elution  X];
          Ytot = [Y_elution  Y];
     
         % Combine elution volume (sweep +extrusion)
          
          tsweep = linspace(1, Ncup, Ncup);
          Vsweep = Vmcup.*tsweep;
          diff = tspan - Ncup;      %Sf*Ncup;
          tex = linspace(1, diff, diff);
          Vext = Vcup.*tex;
          Vext = Vsweep(end)+Vext;
        
           %make elution volume matrix for EECCC 
          V_EECCC = [Vsweep    Vext];
%           T_EECCC = [tsweep  tex];  
        
                 
          for i = 1:n  %Outlet concentrations for n speices
              
              for k = 1:tspan %adding total elution data from X(sweep)+ Y(extrusion)
                    if k <= Ncup
                        M(i,k) = X(Ncup,k,i); %sweep - MP elution
                    else
                        M(i,k) = Y(Ncup,k,i); %extrusion - SP elution
                    end
              end
%               
%               M(i+n,:) = Y(Ncup,:,i);
%               M(i+2*n,:)= X(Ncup,:,i);
              
          end
          
          Ct_out = [V_EECCC; M]';  % [Vol Ctot] column outlet profiles
                


end

          
     
  