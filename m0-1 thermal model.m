%% 
clear all
%% Dust mental
lo=1.25e3;       %% Mass density (kg/m^3)
lo=730; % new
kk=6e-4;         %% thermal conductivity. (W/m K)
Cp=670;           %%specific heat [Joule/kg.K]
Cp=1300;  % new

%% Rock
% lo=3;          %% density (g/cm**3)
% kk=6e-3;       %% thermal cond. (W/cm K)
% Cp=670          %%specific heat [Joule/kg.K]
% Ti=1000;        %% Thermal inertia (MKS) (J/K m**2s**0.5)
%% Constants CERES_FANALE
ffd=72;
ftd=36;
dthi=pi/ftd;                   % thita的間隔
dfi=2*pi/ffd;                    % fi的間隔
z=12.3;                         % 地軸傾角[度]
z=z/180*pi;
P=9.074*60*60;               % 週期 [s]
n=5;                        % 每個單位格點的熱傳導次數(每個單位格點經過幾個dt)
dt=P/ffd/n;                       % 每次熱傳導經過時間[sec]


courant = 0.1;

Ti = sqrt(lo*kk*Cp);            %% Thermal inertia (MKS)(sqrt(lo*kk*Cp)) (J/(K m^2 s^0.5))
dr = sqrt(courant*2*dt*kk/lo/Cp)        % 每個格點的深度[m]       

courant=dr^2/dt/2/kk;
if(courant>=1),disp('Courant OK! : D');else disp('Courant failed!: (');end% Courant-Friedrichs-Levy condition
disp(courant);

depth = 0.05;                  % 溫度會變化的深度範圍[m]
h = ceil(depth/dr)                       % 深度軸的格點數
Ri=1.436e6;                  % Radius of Iapetus[m]
sigma=5.6704*1e-8;           % Stefan-Boltzman constant [W/s/m^2/K^4]
De=1;                        % Distance from Earth to Sun [AU]
Dc=2.6;                   % Distance from Ceres to Sun [AU]
Si=1368*((De/Dc)^2);         % Solar constant at Ceres [J/s/m^2]

% rr=linspace(0,depth,h);

Torr2Pa = 133.32237;
calj=4.2;          %cal=4.2 Joule
mole=6.02e23;
K=8.3144/mole;   % ideal gas constant
mmass = 18e-3/mole;      % molar weight
p_mat(1) = -2445.5646;
%------ parameter of H2O --------------
        p_mat(2) = -6.757169;
        p_mat(3) = -0.01677006;
        p_mat(4) = 1.20514e-5;
        p_mat(5) = 8.2312;
        tran_unit_p = Torr2Pa;
        LH_mat(1) = 0;
        LH_mat(2) = 12420;
        LH_mat(3) = -4.8;
        LH_mat(4) = 0;
        LH_mat(5) = 0;
        tran_unit_LH = calj;
%--------------------------------------


%% Comet thermal Constant
%bn=0.7;                      % Knudsen Number
%AA=3.56e13;                  % dynes/cm^2
%BB=6141.667;                 % Kelvin
%p=1;                         % porosity
%Pv=AA*exp(-BB/T(typ,test))
%ep=                          %  infrared emissitivity
%dg=                          % Intergrain distance of the smallest
                              % residual grain in the mantle
%xai=                         % the mass of dust to ice
%loi=                         % density of ice
%lod=                         % density of silicate
%ep=1;                % infrared emissivity
%r0=10e-6;            % pore radius ( capillary radius ) [m]
%f=0.5;               % surface fraction of ice
%taw=2;               % tortuosity 
%los=2.79e3;          % density of  asteroid regolith [ kg/cm^3 ]
%bn=0.7;              % Knudsen Number
%xai=1.25/0.8;        % mass of dust to ice
%loi=0.8;              % density of ice
%lod=2.0;              % density of silicate
%fa=(1/(1+xai*loi/lod))^(2/3);  % areal factor of 
%C1=4*sigma*ep*r0;
%C2=bn*r0*AA*((8*Ru/pi/M)^0.5);
%C3=3/2*fa*AA*((Ru/2/pi/M)^0.5);
%% constants(0714)
A=0.02;
dF=Si*(1-A);                  % Flux per unit area
a=zeros(h,h+2);
b=zeros(h,h+2);
R=kk*dt/lo/2/Cp/((dr)^2)
CC=1/dr;
            for i=1:h
                a(i,i)=-R;
                a(i,i+1)=1+2*R;
                a(i,i+2)=-R;
            end
            a(:,1)=[];
            a(:,h+1)=[]; 
            a(1,2)=0;
            a(1,1)=1;
            a(h,h)=CC;
            a(h,h-1)=-CC;
            for i=1:h
                b(i,i)=R;
                b(i,i+1)=1-2*R;
                b(i,i+2)=R;
            end                 
            b(:,1)=[];
            b(:,h+1)=[];       
            b(1,1)=1;
            b(1,2)=0;
            b(h,h)=CC;
            b(h,h-1)=-CC;
            a=sparse(a);
            b=sparse(b);
            p_cond=inv(a)*b;   
clear a b
% Solor Flux
td=1;
for td=1:ftd
    for fd=1:ffd
        thi=(td-0.5)/ftd*pi-pi/2;
        fi=(fd-0.5)/ffd*2*pi;
        if(fi>=pi)
            Flux(td,fd)=0;
            continue;
        end
        Flux(td,fd)=dF*cos(thi)*sin(fi);
    end
end

imagesc(Flux)
%% %%%%%%%%%%%%%%%%%%% Main Program%%%%%%%%%%%%%%%%%%%%%%

T=220*ones(h,ftd,ffd);  %Temperature(depth,degree of thita,d.. fi)
for td = 1 : ftd
    for fd = 1 : ffd
        if(Flux(td, fd)>0)
        T(:,td, fd) = sqrt(sqrt(Flux(td, fd)/sigma));
        oT = T(:,td, fd);
        else
        T(:,td, fd) = oT;
        end
    end
end
plotT(:,:)=T(1,:,:);figure;imagesc(plotT);
colorbar;
qT=zeros(h,1);
Flux=abs(Flux);                       %計算表面溫度
%% calaulate for many times
hold on;
timeT=zeros(ftd,ffd);


for td = 1 : ftd
    qT(:) = T(:,td,1);
    for turn=1:1
    for fd = 1 : ffd
    F = Flux(td, fd);
for hh = 1 : n
    oT=T(2,td,fd);
    % 2-divided method for finding roots
    T0=50;
    sb_T = T0;
    p = p_mat(1)/sb_T + p_mat(2)+ p_mat(3)*sb_T + p_mat(4)*(sb_T ^ 2) + p_mat(5)*log10(sb_T);
    LH = LH_mat(1)/sb_T + LH_mat(2)+ LH_mat(3)*sb_T + LH_mat(4)*(sb_T ^ 2) + LH_mat(5)*log10(sb_T);
    fT0 = F + kk/dr*oT - kk/dr*T0 - sigma*(T0^4) - tran_unit_LH * LH * sqrt(mmass/(2*pi*K*T0)) * tran_unit_p * (10^p);
    T1=300;
    sb_T = T1;
    p = p_mat(1)/sb_T + p_mat(2)+ p_mat(3)*sb_T + p_mat(4)*(sb_T ^ 2) + p_mat(5)*log10(sb_T);
    LH = LH_mat(1)/sb_T + LH_mat(2)+ LH_mat(3)*sb_T + LH_mat(4)*(sb_T ^ 2) + LH_mat(5)*log10(sb_T);
    fT1 = F + kk/dr*oT - kk/dr*T1 - sigma*(T1^4) - tran_unit_LH * LH * sqrt(mmass/(2*pi*K*T1)) * tran_unit_p * (10^p);
for m = 1:30
    if(T1>=10000)
        disp('wrong')
        break;
    end
    if(fT0*fT1>0)
        T1=T1*5;
        fprintf('again in time:%f,fd:%f,T1:%f,fT1:%f,fT0:%f\n',time,fd,T1,fT1,fT0)
        continue;
    end
    T2=(T0+T1)/2;
    sb_T = T2;
    p = p_mat(1)/sb_T + p_mat(2)+ p_mat(3)*sb_T + p_mat(4)*(sb_T ^ 2) + p_mat(5)*log10(sb_T);
    LH = LH_mat(1)/sb_T + LH_mat(2)+ LH_mat(3)*sb_T + LH_mat(4)*(sb_T ^ 2) + LH_mat(5)*log10(sb_T);
    fT2 = F + kk/dr*oT - kk/dr*T2 - sigma*(T2^4) - tran_unit_LH * LH * sqrt(mmass/(2*pi*K*T2)) * tran_unit_p * (10^p);
    if(fT0*fT2<=0)
        T1=T2;
        fT1=fT2;
    else
        T0=T2;
        fT0=fT2;
    end
    if(abs(T0-T1)<=0.01)
        break;
    end
end

qT(1) = T2;
qT = p_cond * qT;
T(:, td, fd) = qT(:);
end
timeT(td, fd) = T(1, td, fd);

    end
if(td==1),plot(qT);hold on;end
end

fprintf('%d\n', td);

end
fprintf('done!\n');

T(:,:,ffd+1)=T(:,:,1);
T(:,:,1)=[];
%%
figure;
imagesc(timeT)
colorbar

%%
td = 1;
figure;
for fd = 1 : ffd
    qT(:) = T(:, td, fd);
    plot(qT);
    hold on;
end
%% Write it on file
plotT=zeros(ftd,ffd);
plotT(:,:)=T(1,:,:);
fid=fopen('Ceres_T_VP.dat','w');
fprintf(fid,'%d\n',ftd);
fprintf(fid,'%d\n',ffd);
for td=1:ftd
    for fd=1:ffd
        fprintf(fid,'%8.6f\n',plotT(td,fd));
    end
end
fclose(fid);
%% To make a beautiful plot
plotT=zeros(ftd,ffd);
plotT(:,:)=T(1,:,:);
imagesc(plotT)
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'xtick',[gcapoint1]);
set(gca,'ytick',[gcapoint2]);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'-90','0','90'});
title('P=2[hours]')
xlabel('Local Time');
ylabel('Longitude');
zlabel('Temperature[K]');
colorbar
%% To make a beautiful sphere
td=1;
for thi=linspace(-pi/2,pi/2,ftd)
    td=td+1;
    fd=1;
    for fi=linspace(0,2*pi,ffd)
        fd=fd+1;
        x(td,fd)=sin(fi)*cos(thi);
        y(td,fd)=cos(fi)*cos(thi);
        z(td,fd)=sin(thi);
    end
end

figure
surf(x,y,z,plotT)
shading flat
hidden off;
grid off
colorbar
