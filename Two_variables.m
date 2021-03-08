%% Cantilever Beam Cross Section Optimization %%
% Copyrights: Theyagarajan SELVAM, Sarvesh RAVICHANDRAN, Rajesh Kannan RAVINDRAN Centrale Nantes

% Cantilever Beam Using Two Variables
clear;
close all;
clc;

%% initial parameters
E = 2.1*1e+11;  % Young modulus(pa)
rho  = 7800;    % Density(kg/m^3)
sa= 21*1e+8;    % Admissible stress(Pa)
L = 2;          % Length of beam(m)
F = 20000;      % Concentrated Force(N)
da = 0.01;      % Admissible deflection(m)
ma = 300;       % Admissible mass(Kg)


Lb = 0 ;        %Lower Bound
Ub = 0.3;       %Upper bound

x0 = [0.0001 0.001]; %Guess Values
%% Selection Menu
choice = menu('Select the task to be performed','Rectangular Section','Disc Section','Plus Section','Square Section');

if choice == 1
    recbox = questdlg('Select the minimization to be performed in Rectangular section',...
                      'Minimization',...
                      'Mass Min','Deflection Min','');
    if(strcmp(recbox,'Mass Min')==1)
        obj = @(x)Objectif(x,L,rho,F,E,recbox,choice);
        const = @(x)constraint(x,sa,F,L,E,da,ma,rho,recbox,choice);
    elseif(strcmp(recbox,'Deflection Min')==1)
        obj = @(x)Objectif(x,L,rho,F,E,recbox,choice);
        const = @(x)constraint(x,sa,F,L,E,da,ma,rho,recbox,choice);
    end
elseif choice == 2
        Discbox = questdlg('Select the minimization to be performed in Disc section',...
                       'Minimization',...
                       'Mass Min','Deflection Min','');
        if(strcmp(Discbox,'Mass Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Discbox,choice);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Discbox,choice);
        elseif(strcmp(Discbox,'Deflection Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Discbox,choice);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Discbox,choice);
        end
elseif choice == 3
        Plusbox = questdlg('Select the minimization to be performed in Plus Section',...
                       'Minimization',...
                       'Mass Min','Deflection Min','');
        if(strcmp(Plusbox,'Mass Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Plusbox,choice);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Plusbox,choice);
        elseif(strcmp(Plusbox,'Deflection Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Plusbox,choice);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Plusbox,choice);
        end
elseif choice == 4
     Squarebox = questdlg('Select the minimization to be performed in Square Section',...
                       'Minimization',...
                       'Mass Min','Deflection Min','');
        if(strcmp(Squarebox,'Mass Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Squarebox,choice);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Squarebox,choice);
        elseif(strcmp(Squarebox,'Deflection Min')==1)
            obj = @(x)Objectif(x,L,rho,F,E,Squarebox,choice);
            const = @(x)constraint(x,sa,F,L,E,da,ma,rho,Squarebox,choice);
        end
end
%% Computes the minimization
options = optimset('Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter',1000, 'MaxFunEvals', 1000);
[x] = fmincon(obj, x0, [], [], [], [], Lb, Ub, const, options);

x

%% Contour plotting 
    b = 0:0.001:0.3;
    h = 0:0.001:0.3;
 if choice == 1
   if(strcmp(recbox,'Mass Min')==1)
       [X,Y] = meshgrid(b,h);
       z = Mass_plot(X,Y);
       figure(1);
       meshc(X,Y,z);
       figure(2);
       contourf(X,Y,z);
       colorbar
       hold on;
       
     
       sigma = (6*F*L./(sa*Y.^2)); 
       plot(Y,sigma,'*');
       ylim([0 0.3]);
       hold on;
    elseif(strcmp(recbox,'Deflection Min')==1)
       [X,Y] = meshgrid(b,h);
       z = Mass_plot(X,Y);
       figure(1);
       meshc(X,Y,z);
       figure(2);
       contourf(X,Y,z);
       colorbar
       hold on;
       deflection = (4*F*L^3./(da*x.^3*E));
       plot(x,deflection);
       ylim([0 0.3]);
       hold on;
   end
 end
    
%% Graph Mass
  function m = Mass_plot(X,Y)
        S = X.*Y; % (S= section)
        m = 7800*2.*S; %mass Formula
  end


%% Objective Function

function m = Objectif(x,L,rho,F,E,Button,choice)
    if(strcmp(Button,'Mass Min')== 1)
        if choice == 1 %Rectangular Section
             S = x(1)*x(2); % (S= section)
        elseif  choice == 2 %Disc Section 
            S = pi*((x(1)^2)-(x(2)^2)/4); % (S= section)  
        elseif choice == 3  %plus Section
            S = (2*x(1)*x(2))-(x(1)^2);   % (S= section)
        elseif choice == 4   %square Section 
            S = ((x(1)^2)-(x(2)^2));     % (S= section)
        end
        m = rho*L*S; % Mass Formula
    elseif(strcmp(Button,'Deflection Min')== 1)
        if choice == 1 %rectangular Section
            I = x(1)*(x(2)^3)/12; 
        elseif  choice == 2 %Disc Section
           I = pi*((x(1)^4)-(x(2)^4))/64;   
        elseif choice == 3  %plus Section
            I = (x(1)*(x(2)^3)+((x(1)^3)*x(2))-(x(1)^4)/12);    
        elseif choice == 4   %square section
           I = ((x(1)^3)-((x(2)^2)*x(1)))/2;   
        end
           
        m = F*L^3/(3*E*I); % Mass Formula 
    end
end
%% Constraint Function

function [c,ceq] = constraint(x,sa,F,L,E,da,ma,rho,Button,choice)  
    
       c(1) = stress(x,F,L,choice)-sa;    %Maximumstress
    if(strcmp(Button,'Mass Min')==1)
       c(2) = def_min(x,F,L,E,choice)-da;  % Maximum Deflection
    elseif(strcmp(Button,'Deflection Min')==1)
         c(3) = Massmin(x,L,rho,choice)-ma; %Maximum Mass
    end
    ceq = [];
end

%% Moment of Inertia
function Irect = I_rect(x) %Rectangle Section
           b = x(1);
           h = x(2);      
           Irect = b*(h^3)/12;
end 
function Idisc = I_disc(x) % Disc section
         D = x(1);
         d = x(2);
         Idisc = pi*((D^4)-(d^4))/64;
end 
function Ip = I_plus(x) %Plus Section
         e = x(1);
         h = x(2);
         Ip = (e*(h^3)+((e^3)*h)-(e^4)/12);
end 
function Isq = I_sqr(x) %Square Section
         A = x(1);
         a = x(2);
         Isq =((A^3)-((a^2)*A))/2;
end 
%% Admissible Stress
function sigma = stress(x,F,L,choice) 
        if choice == 1%Rectangular Section
            v = ((x(2)/2));         %Shear force region
            sigma = (F*L/I_rect(x))*v; % Admissible Stress
        elseif  choice == 2 %Disc Section
           v = (x(1))/2;          %Shear force region
            sigma = (F*L/I_disc(x))*v; % Admissible Stress  
        elseif choice == 3  %plus Section
            v = (x(2))/2;        %Shear force region
            sigma = (F*L/I_plus(x))*v;  % Admissible Stress  
        elseif choice == 4  %square Section
            v = (x(1))/2;          %Shear force region
            sigma = (F*L/I_sqr(x))*v;  % Admissible Stress  
        end
             
end
%% Admissible Deflection
function deflection = def_min(x,F,L,E,choice)
     if choice == 1    %rectangukar Section
          deflection = F*L^3/(3*E*I_rect(x)); %Admissible Deflection  
     elseif  choice == 2 %Disc Section
          deflection =  F*L^3/(3*E*I_disc(x)); %Admissible Deflection   
     elseif choice == 3  %plus section
          deflection =  F*L^3/(3*E*I_plus(x));  %Admissible Deflection  
     elseif choice == 4  %square section
          deflection =  F*L^3/(3*E*I_sqr(x));    %Admissible Deflection
     end
        
end
%% Mass 
function m = Massmin(x,L,rho,choice)
     if choice == 1         %rectangukar Section
             S = x(1)*x(2); % (S= section)
        elseif  choice == 2    %disc Section
            S = pi*((x(1)^2)-(x(2)^2)/4);% (S= section)  
        elseif choice == 3       %plus Section
            S = (2*x(1)*x(2))-(x(1)^2); % (S= section)
        elseif choice == 4       %square Section
            S = ((x(1)^2)-(x(2)^2));  % (S= section)  
      end
        m = rho*L*S; % Mass Formula
             
end




