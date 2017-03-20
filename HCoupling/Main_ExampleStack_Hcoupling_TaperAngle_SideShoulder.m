% Behtash Behin-Aein
clear; clc
% Last Update: 9-08-2016
% Used for Extracting scaling laws

% --------------------------------------------------
% This code computes the diameter of the concentric 
% discs for a given input taper angle and desired thickness
%
% -------------------------------------------------
global R t
% CV0 stack
tic
Theta = 60 * pi/180;   % [72 65 58 51]
t_free = 2.58;
R_free = 70/2; % 
f = 1;
t_AP1 = 4.16 * 1.0;   % For %10 higher, multiply Ms_AP1 and t_AP1 by 1.1
% ---------
t_AP2 = 1.7;
t_TL = 0.35;  
t_RL = 0.86; 
% ---------
% Ms in units of 10^6 A/m
V2_POR_mag = 2.36;
V2_POR_tAP1 = 4.16;
V2_POR_Mag_o_tot = V2_POR_mag / V2_POR_tAP1;
%
%
Ms_AP1 = 366 * 1.0 * 1e-6/(t_AP1*1e-7) * (1000 * 1e-6);    
Ms_AP2 = 255 * 1e-6/(t_AP2*1e-7) * (1000 * 1e-6);
Ms_RL =  111 * (95/95) * 1e-6/(t_RL*1e-7)  * (1000 * 1e-6);

t_Ru = 0.2;
t_TB1 = 1.2; 
t_seed = 0;
H = t_TB1 + t_RL + t_TL + t_AP2 + t_Ru + t_AP1;
W = H / tan(Theta);
SS = 1 *(20 - (t_AP2+t_RL+t_TL)/tan(Theta));   % Width of the side shoulder 
%------------
res_nm = 0.3;
%------------
EDR_RL  = 1;% 0.83 with 10; % Etch damage ratio  (lower is more damage)
EDR_AP2 = 1;% 0.83 with 10 % Etch damage ratio  (lower is more damage)
EDR_AP1 = 1;% 1.05 with 10 % Etch damage ratio  (lower is more damage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for FL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_FL = round(t_free/res_nm);   % Number of discs in a magnet
y_a = linspace(t_seed+t_AP1+t_Ru+t_AP2+t_RL+t_TL+t_TB1 , t_seed+t_AP1+t_Ru+t_AP2+t_RL+t_TL+t_TB1+t_free, N_FL);
start = t_seed+t_AP1+t_Ru+t_AP2+t_TL +t_RL+t_TB1;
R_FL = zeros(1,N_FL);
delta_FL = zeros(1,N_FL);
for yy=1:length(y_a)
    t = t_free/N_FL;
    xi = (W/H)*y_a(yy);
    Wi = W - xi;
    R_FL(yy) = R_free + Wi*0; R = R_FL(yy);
    delta_FL(yy) = t*(yy-1);
end
D_FL = 2*R_FL;
sV_FL = sum(R_FL)*t;   % Scaled volume
for cc = 1:N_FL
    figure(1)
    plot(linspace(-R_FL(cc),R_FL(cc),21),start+(cc-1)*t*ones(1,21),'g-','linewidth',8);hold on 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for RL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = round(t_RL/res_nm);   % Number of discs in a magnet
t = t_RL/N;
y_a = linspace(t_seed+t_AP1+t_Ru+t_AP2+t_TL, t_seed+t_AP1+t_Ru+t_AP2+t_TL+t_RL, N);
start = t_seed+t_AP1+t_Ru+t_AP2+t_TL;
R_RL = zeros(1,N);
delta_RL = zeros(1,N);
Hcpl_RL = zeros(1,N);
crt=0;
for yy=1:length(y_a)
    xi = (W/H)*y_a(yy);
    Wi = W - xi;
    R_RL(yy) = R_free + Wi + EDR_RL *SS; R = R_RL(yy);
    delta_RL(yy) = t*(N-yy) + t_TB1;
    Hcpl_RLyFLy = zeros(1,N_FL);
    % --------------------------------------------------------
    for nn=1:N_FL
        R_free_nn = R_FL(nn);
        sV_FLy= R_free_nn*t_free/N_FL;
        Hcpl_RLyFLy(nn) = Hcpl_Disc_Offset(delta_RL(yy)+delta_FL(nn),Ms_RL,R_free_nn,t_free/N_FL,crt) * (sV_FLy/sV_FL);
    end
    Hcpl_RL(yy) = sum(Hcpl_RLyFLy);
    % --------------------------------------------------------
end
Area_RL = ( pi * (R_RL(1).^2) * 1e-14);
mu_RL = Ms_RL * 1e3 * 1e6 * sum( pi * (R_RL.^2) * t *1e-21);
H_RL = sum(Hcpl_RL);
D_RL = 2*R_RL;
for cc = 1:N
    figure(1)
    plot(linspace(-R_RL(cc),R_RL(cc),21),start+(cc-1)*t*ones(1,21),'b-','linewidth',8);hold on 
end 
%%%% ----------------------------------------
%%%% ----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for AP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = round(t_AP2/res_nm);   % Number of discs in a magnet
y_a = linspace(t_seed+t_AP1+t_Ru,t_seed+t_AP1+t_Ru+t_AP2, N);
start = t_seed+t_AP1+t_Ru;
RAP2 = zeros(1,N);
delta_AP2 = zeros(1,N);
Hcpl_AP2 = zeros(1,N);
crt=0;
t = t_AP2/N;
for yy=1:length(y_a)
    xi = (W/H)*y_a(yy);
    Wi = W - xi;
    RAP2(yy) = R_free + Wi + EDR_AP2 *SS; R = RAP2(yy);
    delta_AP2(yy) = t*(N-yy) + t_TB1 + t_RL;
    Hcpl_AP2yFLy = zeros(1,N_FL);
    % --------------------------------------------------------
    for nn=1:N_FL
        R_free_nn = R_FL(nn);
        sV_FLy= R_free_nn*t_free/N_FL;
        Hcpl_AP2yFLy(nn) = Hcpl_Disc_Offset(delta_AP2(yy)+delta_FL(nn),Ms_AP2,R_free_nn,t_free/N_FL,crt) * (sV_FLy/sV_FL);
    end
    Hcpl_AP2(yy) = sum(Hcpl_AP2yFLy); 
end
Area_AP2 = ( pi * (RAP2(1).^2) * 1e-14);
mu_AP2 = Ms_AP2 * 1e3 * 1e6 * sum( pi * (RAP2.^2) * t *1e-21); 
H_AP2 = sum(Hcpl_AP2);
D_AP2 = 2*RAP2; 
for cc = 1:N
    figure(1)
    plot(linspace(-RAP2(cc),RAP2(cc),21),start+(cc-1)*t*ones(1,21),'b-','linewidth',8);hold on 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for AP1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = round(t_AP1/res_nm);   % Number of discs in a magnet
y_a = linspace(t_seed,t_seed+t_AP1,N);
start = t_seed;
RAP1 = zeros(1,N);
delta_AP1 = zeros(1,N);
Hcpl_AP1 = zeros(1,N);
crt=0;
t = t_AP1/N;
for yy=1:length(y_a);
    xi = (W/H)*y_a(yy);
    Wi = W - xi;
    RAP1(yy) = R_free + Wi + EDR_AP1 *SS; R = RAP1(yy);2*R
    delta_AP1(yy) = t*(N-yy) + t_Ru + t_AP2 + t_RL + t_TB1; 
    Hcpl_AP1yFLy = zeros(1,N_FL);
    % --------------------------------------------------------
    for nn=1:N_FL
        R_free_nn = R_FL(nn);
        sV_FLy= R_free_nn*t_free/N_FL;
        Hcpl_AP1yFLy(nn) = Hcpl_Disc_Offset(delta_AP1(yy)+delta_FL(nn),Ms_AP1,R_free_nn,t_free/N_FL,crt) * (sV_FLy/sV_FL);
    end
    Hcpl_AP1(yy) = sum(Hcpl_AP1yFLy);
end
Area_AP1 = ( pi * (RAP1(1).^2) * 1e-14);
mu_AP1 = Ms_AP1 * 1e3 * 1e6 * sum( pi * (RAP1.^2) * t *1e-21); 
H_AP1 = sum(Hcpl_AP1);
D_AP1 = 2*RAP1;
for cc = 1:N
    figure(1)
    plot(linspace(-RAP1(cc),RAP1(cc),21),start+(cc-1)*t*ones(1,21),'r-','linewidth',8);hold on 
end
figure(1)
ylim([-1 15])
xlim([-70 70])
set(gcf,'color','w')
set(gca,'fontsize',24)
xlabel('Width  [nm]','fontsize',24)
ylabel('Height [nm]','fontsize',24)
%
%%%%%%%%%%%%%%%%%%%%
BitCount = 7.4e8;
moment_AP2 = (mu_AP2 + mu_RL) * BitCount
mu_AP2_Film = moment_AP2 / (BitCount * Area_AP2);  % mu-emu/cm^2
%----------------------------------------
moment_AP1 = mu_AP1 * BitCount 
mu_AP1_Film = moment_AP1 / (BitCount * Area_AP1);  % mu-emu/cm^2
%----------------------------------------
moment_ratio_2o1 = moment_AP2 / moment_AP1
H_tot_cpl = H_RL + H_AP2 - H_AP1

toc