%% Cantilever Beam Cross Section Optimization %%
% Copyrights: Theyagarajan SELVAM, Sarvesh RAVICHANDRAN, Rajesh Kannan RAVINDRAN Centrale Nantes

% Cantilever Beam Using One Variable

clear;
close all;
clc;

% initial parameters
E = 2.1*1e+11;  % Young modulus(pa)
rho  = 7800;    % Density(kg/m^3)
sa= 21*1e+8;    % Admissible stress(Pa)
L = 2;          % Length of beam(m)
F = 20000;      % Concentrated Force(N)
da = 0.01;      % Admissible deflection(m)
ma = 300;       % Admissible mass(Kg)


Lb = 0 ;   % Lower Bound
Ub = 0.4;  % Upper Bound

x0 = 0.01;  %Guess Values
%% Selection Menu
Button = questdlg('Select the task to be performed',...
                   'Optimization',...
                   'Square Section','Disc Section', '');
if(strcmp(Button,'Square Section')==1)
    sqrbox = questdlg('Select the minimization to be performed in Square section',...
                      'Minimization',...
                      'Mass Min','Deflection Min','');
    if(strcmp(sqrbox,'Mass Min')==1)
        obj = @(x)Objectif(x,L,rho,F,E,Button,sqrbox);
        const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,sqrbox);
    elseif(strcmp(sqrbox,'Deflection Min')==1)
        obj = @(x)Objectif(x,L,rho,F,E,Button,sqrbox);
        const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,sqrbox);
    end
elseif(strcmp(Button,'Disc Section')==1)
    discbox = questdlg('Select the minimization to be performed in Disc section',...
                       'Minimization',...
                       'Mass Min','Deflection Min','');
        if(strcmp(discbox,'Mass Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Button,discbox);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,discbox);
        elseif(strcmp(discbox,'Deflection Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Button,discbox);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Button,discbox);
        end
end
 %% Computes the minimization   
options = optimset('Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter',1000, 'MaxFunEvals', 1000);
x = fmincon(obj, x0, [], [], [], [], Lb, Ub, const, options);
x

%% Graph Plotting
% Plotting the mass curve for square 
if(strcmp(Button,'Square Section')==1)
    if(strcmp(sqrbox,'Mass Min')==1)
        x = linspace(Lb,Ub);
        massplot = (rho*L*(x.^2));
        maxstress = ((6*F*L)/(sa))^(1/3);
        maxdefl = ((4*F*L^3)/(E*da))^(1/4);
        maxmass = ((ma)/(rho*L))^(1/2);
        plot(x,massplot);
        hold on;
        plot (maxstress,(((sa)/(rho*L))^(1/2)),'*');
        hold on;
        plot (maxdefl,(((da)/(rho*L))^(1/2)),'*');

        title('mass minimisation curve for square');
        xlabel('square side length(a)');
        ylabel('Mass');

    % Plotting the deflection curve for square
    elseif(strcmp(sqrbox,'Deflection Min')==1)
        x = linspace(Lb,Ub);
        deflplot = ((4*F*L^3)./(E*(x.^4)));
        maxstress = ((6*F*L)/(sa))^(1/3);
        maxdefl =   ((4*F*L^3)/(E*da))^(1/4);
        maxmass = ((ma)/(rho*L))^(1/2);
        plot(x,deflplot);
        hold on;
        plot (maxstress,((4*F*L^3)/(E*sa))^(1/4),'+');
        hold on;
        plot (maxmass,((4*F*L^3)/(E*ma))^(1/4),'+');

        title('deflection minimisation curve for square');
        xlabel('square side length(a)');
        ylabel('Deflection');
    end
elseif(strcmp(Button,'Disc Section')==1) %plotting mass curve for disc
    if(strcmp(discbox,'Mass Min')==1)
        x = linspace(Lb,Ub);
        massplot = (rho*L*pi*(x.^2)/4);
        maxstress = ((32*F*L)/(pi*(sa))^(1/3));
        maxdefl = ((64*F*L^3)/(3*E*pi*da))^(1/4);
        maxmass = ((ma)/(rho*L*pi))^(1/4);
        plot(x,massplot)
        xlim([-100 800]);
        ylim([-100 2000]);
        hold on;
        plot (maxstress,(((sa)/(rho*L*pi))^(1/4)),'*');
        hold on;
        plot (maxdefl,(((da)/(rho*L*pi))^(1/4)),'*');
        hold on;
        title('mass minimisation curve for Disc');
        xlabel('Diameter of Disc(D)');
        ylabel('Mass');

    % Plotting the deflection curve for disc 
      elseif(strcmp(discbox,'Deflection Min')==1)
        x = linspace(Lb,Ub);
        deflplot = ((64*F*L^3)./(3*E*pi*(x.^4)));
        maxstress = ((32*F*L)/(pi*(sa))^(1/3));
        maxdefl = ((64*F*L^3)/(3*E*pi*da))^(1/4);
        maxmass = ((ma)/(rho*L*pi))^(1/4);
        plot(x,deflplot);
        xlim([-100 800]);
        ylim([-100 2000]);
        hold on;
        plot (maxstress,((64*F*L^3)./(3*E*pi*sa))^(1/4),'+');
        hold on;
        plot (maxmass,((64*F*L^3)./(3*E*pi*ma))^(1/4),'+');

        title('deflection minimisation curve for Disc');
        xlabel('Diameter of Disc(D)');
        ylabel('Deflection');
    end
end

%% Objective Function
function m = Objectif(x,L,rho,F,E,Button,Button2)
   if(strcmp(Button,'Square Section')==1) %Square Section Mass Minimization 
        if(strcmp(Button2,'Mass Min')== 1)
             S = x^2; % (S= section)
             m = rho*L*S; %(Mass of Square)
        elseif(strcmp(Button2,'Deflection Min')== 1)
             I = (x^4)/12; %Moment of Inertia
             m = F*L^3/(3*E*I);% Deflection of Square
        end
    elseif(strcmp(Button,'Disc Section')==1)
         if(strcmp(Button2,'Mass Min')== 1)
             S = pi*(x^2)/4; % (S= section)
             m = rho*L*S;  % Mass
        elseif(strcmp(Button2,'Deflection Min')== 1)
             Ici = pi*(x^4)/64; %Moment of Inertia
             m = F*L^3/(3*E*Ici); % Deflection of Square
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
function Isq = I_sq(x) % Square Section
         a = x;
         Isq =a^4/12;
end 

function Ici = I_di(x) % Disc Section
         D = x;
         Ici = pi*(D^4)/64;
end 
%% Admissible Stress
function sig = stress(x,F,L,Button) 
        if(strcmp(Button,'Square Section')==1)
           v = (x/2);  %Shear force region
           sig = (F*L/I_sq(x))*v; % Admissible Stress 
        elseif(strcmp(Button,'Disc Section')==1)
           v = (x/2); %shear force region 
           sig = (F*L/I_di(x))*v; % Admissible Stress
        end
         
end

%% Admissible Deflection
function deflection = def_min(x,F,L,E,Button)
       if(strcmp(Button,'Square Section')==1)
        deflection = F*L^3/(3*E*I_sq(x)); %Admissible Deflection
       elseif(strcmp(Button,'Disc Section')==1)
         deflection =  F*L^3/(3*E*I_di(x)); %Admissible Deflection
       end
end


%% Mass 
function m = Massmin(x,L,rho,Button)
     if(strcmp(Button,'Square Section')==1)
         S = x^2;  %(Section)
     elseif(strcmp(Button,'Disc Section')==1)
         S = pi*(x^2)/4;  %(Section)
     end
     m = rho*L*S; %Mass  Formula
end





