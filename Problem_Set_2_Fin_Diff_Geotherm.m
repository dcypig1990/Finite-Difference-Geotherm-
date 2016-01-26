% Finite Difference Model of Geotherm + Surface Temperature Warming
% Jan 25 2016

%% Data from Borehole

T_data=[-5.320000,-5.320000,-5.290000,-5.240000,-5.110000,-4.930000,...
   -4.670000,-4.440000,-4.140000,-3.840000,-3.560000,-3.240000,...
   -2.950000,-2.620000,-2.350000,-2.040000,-1.780000,-1.480000,...
   -1.180000,-0.9100000,-0.6000000,-0.2800000,7.0000000E-02];

z_data=[27.73685,42.97689,58.21692,73.45695,88.69698,103.9370,...      
119.1770,134.4171,149.6571,164.8971,180.1372,195.3772,...     
210.6172,225.8573,241.0973,256.3373,271.5773,286.8174,...      
302.0574,317.2974,332.5375,347.7775,363.0175];   

%% Initialize Finite Difference Function Using Geotherm Gradient

%--parameters for finite difference code

N = 23; %number of points
dz = 15.0; %spacing in m
k = 2; %W/mK
rho = 2000; %kg/m3
Cp = 2000; % J/kg K
kappa = k/(rho*Cp);
dt = (dz^2/kappa)*0.2; %s

%--parameters for geotherm gradient

z=0:dz:(dz*(N-1));
Ts_int=-6.75; %degrees C
k=2.5; %W/(m-K)
Qm=0.045; %W/m2

%--initial T equation+changing T equation

T_int = transpose(Ts_int+(Qm/k)*z);
T = transpose(Ts_int+(Qm/k)*z);

%--surface temperature change

%dT=0:0.25:2.5; %degrees C
dT=1.5; %degrees C
T(1) = Ts_int+dT; %degrees C

%--establish heat gradient

dTdz = zeros(N,1);
dTdz(N) = 0.025; %K/m

%--time in years of run
x=(dt*10)/(365*3600*24);

%% Run Finite Difference Code at 10 timesteps of dt

for t=1:10
    
    %--if I want to change surface temperature with every timestep
    
    %T(1) = Ts_int+dT(t); %degrees C
  
    %--calculate T gradient

    dTdz(1:N-1) = diff(T)/dz;

    %--calculate heat flux

    q = -k*dTdz;

    %--rate of change of T
    
    dqdz= diff(q)/dz;

    %--update T...top T is fixed
    
    T(2:N) = T(2:N) - (1/(rho*Cp))*dqdz*dt;
    
    %--if I want to plot for every timestep 
    
    %figure(1)
    %plot(T,z,'r','linewidth',2)
    %hold on
    
end

%% Plots

figure(1)

plot(T,z,'r','linewidth',2)

hold on

plot(T_int,z,'g','linewidth',2)
plot(T_data,z_data,'b','linewidth',2)

xlabel('Temperature (C)','fontname','arial','fontsize',21)
ylabel('Depth (m)','fontname','arial','fontsize',21)
set(gca,'fontsize',18,'fontname','arial')
set(gca,'YDIR','reverse')

legend('Finite Difference Model: 1.5 degrees warming over 28 years','Initial Geotherm Gradient','Observed Borehole Geotherm','Location','northeast'); 

hold off

%% Goodness of Fit

%--to get matrix dimensions to agree

T_data=transpose(T_data);

%--chi squared test

F=sum((T_data-T).^2);

disp(F)

%% Comments on Code

%To whomever is looking at my code, do you have any suggestions
%as to how to input specific z values into the finite difference
%model in order to compare the exact same z values between
%modeled and observed data?? I just estimated the observed
%depth values to be equally spaced at 15 m, and used that 
%value in the model.
%Also, do you have any suggestions for
%how to make my model better fit the observations?? I am not 
%sure I did the surface temperature change correctly...





