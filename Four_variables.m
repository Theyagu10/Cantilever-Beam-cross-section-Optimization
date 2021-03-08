%% Cantilever Beam Cross Section Optimization %%
% Copyrights: Theyagarajan SELVAM, Sarvesh RAVICHANDRAN, Rajesh Kannan RAVINDRAN Centrale Nantes

% Cantilever Beam Using Four Variables
clear;
close all;
clc;

%% Initial parameters
E = 2.1*1e+11;  % Young modulus(pa)
rho  = 7800;    % Density(kg/m^3)
sa= 21*1e+8;    % Admissible stress(Pa)
L = 2;          % Length of beam(m)
F = 20000;      % Concentrated Force(N)
da = 0.01;      % Admissible deflection(m)
ma = 300;       % Admissible mass(Kg)


Lb = 0 ;        %Lower Bound
Ub = 0.4;       %Upper bound

x0 = [0.1 0.03 0.26 0.07]; % Guess values % for all defl


%% Selection Menu
Button = questdlg('Select the task to be performed',...
                   'Optimization',...
                   'Rectangle Section','I Section', '');
if(strcmp(Button,'Rectangle Section')==1) 
    recbox = questdlg('Select the minimization to be performed in Rectangle section',...
                      'Minimization',...
                      'Mass Min','Deflection Min','');
    if(strcmp(recbox,'Mass Min')==1)
        obj = @(x)Objectif(x,L,rho,F,E,Button,recbox); 
        const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,recbox);
    elseif(strcmp(recbox,'Deflection Min')==1)
        obj = @(x)Objectif(x,L,rho,F,E,Button,recbox);
        const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,recbox);
    end
elseif (strcmp(Button,'I Section')==1)
    Isectionbox = questdlg('Select the minimization to be performed in I Section',...
                       'Minimization',...
                       'Mass Min','Deflection Min','');
        if(strcmp(Isectionbox,'Mass Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Button,Isectionbox);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,Isectionbox);
        elseif(strcmp(Isectionbox,'Deflection Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Button,Isectionbox);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,Isectionbox);
        end
end

%% Computes the minimization
options = optimset('Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter',1000, 'MaxFunEvals', 1000);
x = fmincon(obj, x0, [], [], [], [], Lb, Ub, const, options);

x

%% Objective Function
function m = Objectif(x,L,rho,F,E,Button,Button2)
    if(strcmp(Button,'Rectangle Section')==1)
        if(strcmp(Button2,'Mass Min')==1)   %Rectangle mass minimization
             S = (x(1)*x(2))-(x(3)*x(4));
             m = rho*L*S;                    
        elseif(strcmp(Button2,'Deflection Min')==1)%Rectangle Deflection minimization
             I = ((x(1)*(x(2)^3))-((x(3)*x(4))^3))/12;
             m = F*L^3/(3*E*I); % (S= section)
        end        
    elseif (strcmp(Button,'I Section')==1)
        if(strcmp(Button2,'Mass Min')==1)
             S = (x(1)*x(2))-(x(3)*(x(1)-x(4)));%I-Section mass minimization
             m = rho*L*S;
        elseif(strcmp(Button2,'Deflection Min')==1) %I-Section Deflection minimization
             I = ((x(1)*(x(2)^3))-(x(3)^3*(x(1)-x(4))))/12;
             m = F*L^3/(3*E*I); % (S= section)
        end
    end
end

%% Constraint Function
function [c,ceq] = constraint(x,sa,F,L,E,da,ma,rho,Button,Button2)  
    
       c(1) = stress(x,F,L,Button)-sa;  %Maximumstress
    if(strcmp(Button2,'Mass Min')==1)   
       c(2) = def_min(x,F,L,E,Button)-da; % Maximum Deflection
    elseif(strcmp(Button2,'Deflection Min')==1)
         c(3) = Massmin(x,L,rho,Button)-ma; %Maximum Mass
    end
    ceq = [];
end

%% Moment of Inertia
function Ire = I_rect(x) % rectangle
         B = x(1);
         H = x(2);
         b = x(3);
         h = x(4);
         Ire =((B*(H^3))-((b*h)^3))/12;
end 
function Ii = I_I(x)   % I section
         B = x(1);
         H = x(2);
         h = x(3);
         b= x(4);
         Ii =((B*(H^3))-(h^3*(B-b)))/12;
end 
%% Admissible Stress
function sig = stress(x,F,L,Button)  
       if(strcmp(Button,'Rectangle Section')==1)
         v = (x(2)/2); %Shear force region
         sig = (F*L/I_rect(x))*v; % Admissible Stress 
       elseif(strcmp(Button,'I Section')==1)
          v = (x(2)/2); %shear force region 
          sig = (F*L/I_I(x))*v;  % Admissible Stress
       end
       
end
%% Admissible Deflection
function deflection = def_min(x,F,L,E,Button)
     if(strcmp(Button,'Rectangle Section')==1) 
        deflection = F*L^3/(3*E*I_rect(x));%Admissible Deflection
     elseif(strcmp(Button,'I Section')==1)
          deflection = F*L^3/(3*E*I_I(x)); %Admissible Deflection
     end
end


%% Mass 
function m = Massmin(x,L,rho,Button)
     if(strcmp(Button,'Rectangle Section')==1) 
         S = (x(1)*x(2))-(x(3)*x(4)); %(Section)
     elseif(strcmp(Button,'I Section')==1)
         S = (x(1)*x(2))-(x(3)*(x(1)-x(4))); %(Section)
     end
     m = rho*L*S; % Mass Formula
end


























