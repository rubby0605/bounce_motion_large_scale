clear all;
%% Read the atmospheric density
fid=fopen('den_0415_1e4.dat','r');
for ii=1:3
    line = fgetl(fid);
    %str=sscanf(line,'%d');
    str=sscanf(line, '%f');
    if(ii==1),ay=str;end
    if(ii==2),by=str;end
    if(ii==3),cy=str;end
end
ftd=ay;
ffd=by;
hd = cy;
den_atm=zeros(hd,ftd,ffd);
for ht=1:hd
for ii=1:ftd
    for jj=1:ffd
    line = fgetl(fid);
    str=sscanf(line, '%f');
    den_atm(ht,ii,jj)=str;
    end
end
end
fclose(fid);
copy_atm=den_atm;

%% Read the temperature
fid=fopen('Ceres_T_VP.dat','r');
for ii=1:2
    line = fgetl(fid);
    %str=sscanf(line,'%d');
    str=sscanf(line, '%f');
    if(ii==1),ay=str;end
    if(ii==2),by=str;end
end
ftd=ay;
ffd=by;
copyT=zeros(ffd,ftd,ffd);
copyT2=zeros(ftd,ffd);

for ii=1:ftd
    for jj=1:ffd
    line = fgetl(fid);
    str=sscanf(line, '%f');
    copyT2(ii,jj)=str;
    end
end

for ii = 1 : ffd
    copyT(ii,:,:) = copyT2(:,:);
    copyT2(:,ffd + 1) = copyT2(:,1);
    copyT2(:,1)=[];
end



fclose(fid);
figure(3)
%%
plotT=zeros(ftd,ffd);
for tt=1:18
plotT(:,:)=copyT(tt,:,:);
for td=1:ftd
    plotT(td,:)=copyT(tt,ftd-td+1,:);
end
imagesc(plotT(:,:));
pause(0.01)
end
colorbar
imagesc(plotT(:,:));
hold on;
[C,h] = contour(plotT,'w');
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
colormap jet
colorbar
hold on;
colorbar
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint2);
set(gca,'yticklabel',{'Northern Cap','0','Southern Cap'});
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'xtick',gcapoint1);
set(gca,'xticklabel',{'6','12','18','24','6'});
xlabel('Local Time');
ylabel('latitude');
plotT(:,:)=copyT(1,:,:);
clear copyT
%% RK try
% constant
totaln=1e4;
G=6.67e-11;                  % m^3/kg/s^2
K=1.38e-23;                  % Boltzmann constant
M=9.39e20;                   % kg
mole=6.02e23;                % moleculer 
distance=2.8;                % [AU]
lifetime=1/1.02e-5*(2.8)^2           % [sec]
Radii=476.2e3;               % [m] 
mmass=18e-3;                 % H2O
ve=sqrt(2*G*M/Radii);        % escape speed [m/s]
fprintf('Escape Velocity= %f\n',ve);
P=9.07417*60*60;             % 週期 [s]
h=5;    %% delta time
op=0;
cc=0;
%toP=zeros(ftd,ffd);
%toZ=zeros(ftd,ffd);
%copy=zeros(ftd,1);
%copy2=zeros(ftd,1);
alph=0.5; % alph between 0~1, 1 means sticking. alph is the function of temperature.
 % 總共有1.9e17顆水分子，每次模擬代表1e12顆水分子
dt=P/ffd;  %% sec
%wio=ones(ftd,ffd);  %% with water ice or not?
%wa=zeros(ftd,ffd);
%modeln=5e11;
%wv=zeros(ftd,ffd);
%Mr=0.1;    % radius of a meteorite [m]
%water=0.9406*(2*Mr)^2*Mr/5*8.9*1000/18*mole;     % crater 蒸發的water [g]
k=1.38e-23;      % Boltzmann constant
G=6.67e-11;      % M.K.S
mole=6.02*1e23;  % moleculer 
mmass=18e-3;             % M of atmosphere per 1 mole

% Scale height
T_scale = mean(mean(plotT));
g_value = G * M / (Radii ^ 2);

scale_height = k * T_scale / mmass * mole / g_value;

%% Calculate just for the time scale histogram
solar = zeros(ffd, 3);
cal = zeros(1, 3);
horizon = zeros(ftd,ffd);

max_height = 3e6*sqrt(2);
hd = 50;
h_line=zeros(hd, 1);

alpha2 = 0.04;
unih = (max_height-Radii)/((1+alpha2)^(hd)-1)*alpha2;
scale_height / (1+alpha2)^(hd-1)
unih
if( unih >= scale_height * 0.8 / (1+alpha2)^(hd-1)),
    fprintf('Use larger hd\n');
else
    fprintf('hd is %d, unih is %f scale height is %f\n', hd, unih, scale_height);
end
%% collision parameters

dev_td = 3;
dev_fd = dev_td;
v_record = zeros(1e5, 6);
count_nn = zeros(ftd/dev_td, ffd/dev_fd, hd);

%%

for num = 0 : hd - 1
    h_line(num + 1) = ((1 + alpha2)^(num + 1)-1)/(alpha2)*unih+Radii;
end

Ri = max_height-1;
ht = ceil(log(alpha2*(Ri-Radii)/unih+1)/log(1+alpha2))


za = 3/180*pi;%%
disp(za)
thi=-za;
for time=1:ffd
    fi=(time-0.5)/ffd*2*pi+pi;
    solar(time,:)=[cos(thi)*sin(fi) cos(thi)*cos(fi) sin(thi)];
end
%aa=1;
s=zeros(3,1);
%% velocity matrix for different T(0226)!!!!
nnn=1;
num=zeros(23,1);
p_max=zeros(23,1e3);
for MIV_T=linspace(20,240,23)
    n_v=0;
for vv=linspace(0,1000,101)
aa=4*pi*vv^2*sqrt((1/2/pi*mmass/mole/k/MIV_T)^3)*exp(-mmass/mole*3*vv^2/2/k/MIV_T);
aa=aa*5e5;
if(floor(aa)==0),continue;end
p_max(nnn,n_v+1:n_v+1+floor(aa))=vv;
n_v=n_v+floor(aa);
end
n_v=n_v+1;
num(nnn)=n_v;
nnn=nnn+1;
end
disp(num)
num2=num;
%% check
ii = 23;
ppp = p_max(ii, 1 : num(ii));

x = 10 : 10 : 1000; figure;

hist(ppp,x)
xlabel('thermal velocity [m/s]');
ylabel('test numbers');

hold on
den_atm = zeros(hd, ftd, ffd);
%% 

max_nn = 0;
max_rot = 0;
totaln = 1e4;
record_nn = 1;
ouou = 0;
ttime = zeros(8, totaln);
%%
nt = 1;


figure;
%sphere;
hold on;
for n = 1 : totaln
   rest_time = 0;
   ou=0;
   rot=0;
   t_n=0;
%   record_io=0;
   
   nt = 1;
   %if(nt==0),nt=nt+ffd;end
   
   time=0;
   disp(n)
   W=1;
   solar_time=0;
   thi=rand*pi-pi/2;
   ttime(1,n)=thi;
   fi=rand*pi*2;
   flight_time=0;
   s=[Radii*cos(thi)*sin(fi) Radii*cos(thi)*cos(fi) Radii*sin(thi)];
   %traj=zeros(1,3,5e4);
   del_time=0;
 %  vout_first=p_Max(ceil(rand*n_v1));
   v_io=0;
    while(W >= 0.01)
    KF=2*pi*Radii*abs(cos(thi))/P;
    td=ceil((thi+pi/2)/pi*ftd);
    fd=ceil(fi/2/pi*ffd);
        if(fd==0),fd=1;end%if(td==0),td=ftd;end
    local_T=plotT(td,fd);
        if(local_T<100),
            for add=1:ffd
            rot=rot+1;
            %nt=nt+1;
            time=time+dt;
            rest_time = rest_time + dt;
            %disp('1');
            %pause(0.1);
            %sfd=fd-add;
            %if(sfd<=0),sfd=sfd+ffd;end
            %fi=sfd/ffd*pi*2;
            %s(:)=[Radii*cos(thi)*sin(fi) Radii*cos(thi)*cos(fi) Radii*sin(thi)];
            %plot3(s(1),s(2),s(3),'rp');pause(0.1);
            %if(abs(s(1)-os(1))>=1e4),disp(s(1));end
            %os=s;
            %hold on;
           %t_n=t_n+1;
           %traj(1,:,t_n)=s(:);
            %if(nt>ffd),nt=nt-ffd;end
            if(plotT(td,fd)>=100),
                if(fd<=0),fd=fd+ffd;end
                rot_fi = (add-fd)/ffd*2*pi;
                s(:)=[s(1)*cos(rot_fi)-sin(rot_fi)*s(2) sin(rot_fi)*s(1)+cos(rot_fi)*s(2) s(3)];
                fd = add;
                fi=fd/ffd*pi*2;
                break;
            end
            end
            if(max(plotT(td,:)) < 100),horizon(td,fd)=horizon(td,fd)+W;break;end
        end
   flyd=rand*pi+thi-pi/2;
   xi=rand*pi+fi-pi/2;
   %if(abs(flyd)<=0.1),continue;end  ??
   v(1:3)=0;
   local_T=plotT(td,fd);
    %if(v_io==0),v_io=1;vout=vout_first;
    %else
        if(local_T>=220),vr=p_max(21,ceil(rand*num(21)));
        else
        vr=p_max(ceil(local_T/10)-1,ceil(rand*num(ceil(local_T/10)-1)));
        end
        vout=sqrt(alph*vr^2+(1-alph)*sum(v(:).^2));
    %end
   v(1)=cos(flyd)*sin(xi);
   v(2)=cos(flyd)*cos(xi);
   v(3)=sin(flyd);
   v(:)=v(:).*vout;
   v(1)=v(1)+KF*sin(fi+pi/2);  %% velocity induced by rotational force 
   v(2)=v(2)+KF*cos(fi+pi/2);
    if(sqrt(sum(v(:).^2))>=ve),ou=1;break;end
   Ri=sqrt(sum(s(:).^2));%%
   oa=G*M*s(:)/(Ri^3);
   %disp('!')
    for ii=1:1e7
        %disp('2')
        flight_time=flight_time+h;
        time=time+h;
        %if(time>=3e3),break;end
        %pause(0.1)
        for i=1:3
            ns(i)=s(i)+v(i)*h-G*M*s(i)/(Ri^3)*(1.5*(h^2))-oa(i)*(h^2)/6;
            v(i)=v(i)-G*M*s(i)/(Ri^3)*h;
            oa(i)=G*M*s(i)/(Ri^3);
        end
        s(:)=ns(:);
        Ri = sqrt(sum(s(1:3).^2));
        if(Ri <= Radii),break;end
        if(Ri >= max_height),ou=1;break;end
       
        %t_n=t_n+1;traj(1,:,t_n)=s(:);
        
        thi=asin(s(3)/Ri);
        if(s(1)==0 && s(2)>=0), fi = 0 ;
        elseif(s(1)==0 && s(2)<0), fi = pi;
        else
            if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
            else fi=-atan(s(2)/s(1))+pi/2;
            end
        end
        fd=ceil(fi/2/pi*ffd);
        td=ceil((thi+pi/2)/pi*ftd);
        if(td==0),td=1;end
        
        if( fd <= ffd/2),
            W = W * exp(-h/lifetime);
            solar_time = solar_time+h;
        end
        %if(ceil(ii/100)~=ii/100),continue;end
    %if(abs(s(1)-os(1))>=1e4),disp(s(1)*1e4);end
    %os=s;
    
    %if(ii/5e1==ceil(ii/5e1)),plot3(Ri/Radii*cos(thi)*sin(fi),Ri/Radii*cos(thi)*cos(fi),Ri/Radii*sin(thi),'ko');hold on;end
    %if(ii/1e3==ceil(ii/1e3)),pause(0.1);end
    
    
    if(ii==1),continue;end
    ht = ceil(log(alpha2*(Ri-Radii)/unih+1)/log(1+alpha2));
    
        %if(ht>=hd),ou=1;disp(ht);break;end
    td2 = ceil(td / dev_td);
    fd2 = ceil(fd / dev_fd);
    if(max_nn~=td2+fd2+ht && count_nn(td2, fd2, ht) <= 50),
        max_nn = td2+fd2+ht;
        v_record(record_nn, 1 : 3) = v(1:3);
        v_record(record_nn, 4 : 6) = [td2 fd2 ht];
        record_nn = record_nn+1;
        count_nn(td2, fd2, ht) = count_nn(td2, fd2, ht) + 1;
        if(record_nn/1e4 == ceil(record_nn/1e4)),disp(record_nn),end
    end
        den_atm(ht, td, fd) = den_atm(ht, td, fd) + W;
    end
        if(s(1)==0 && s(2)>=0), fi = 0 ;
        elseif(s(1)==0 && s(2)<= 0), fi = pi;
        else
            if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
            else fi=-atan(s(2)/s(1))+pi/2;
            end
        end
        fd=ceil(fi/2/pi*ffd);
        td=ceil((thi+pi/2)/pi*ftd);
        if(td==0),td=1;end
    
    if(ou==1),break;end
    del_time=del_time+h*ii;
    if(td<=2||td>=35),horizon(td,fd)=horizon(td,fd)+0.05*W;W=W*0.95;end
    if(del_time>=dt),del_n=floor(del_time/dt);rot=rot+del_n;fd=fd-del_n;
        if(fd<=0),fd=fd+ffd;end,
        del_time=del_time-dt*del_n;
        continue;
    end
    %if(time>=1e4),ou=1;break;end
    end
    
   
    
   if(ou==1),ouou=ouou+1;continue;
   
   %else break;
   end
   if(rot>=max_rot),max_rot=rot;end
ttime(2,n)=thi;
ttime(3,n)=time;
ttime(4,n)=W;
ttime(5,n)=flight_time;
ttime(8,n)=rest_time;
   %endingpoint(1)=t_n;
%   disp(time)
  
end




%% write "v_record" on file

fid = fopen('v_record_6.dat','w');
nnn = record_nn - 1;
fprintf(fid,'%d\n', record_nn - 1);

for num = 1 : nnn
fprintf(fid,'%f %f %f\n',v_record(num, 1), v_record(num, 2), v_record(num, 3));
fprintf(fid,'%d\n', v_record(num, 4));
fprintf(fid,'%d\n', v_record(num, 5));
fprintf(fid,'%d\n', v_record(num, 6));
end

%% Normalized
%
%den_atm=copy_den;
%
%totaln=totaln;
prn= 1e24;



dfi = 2 * pi / ffd;
dthi = pi / ftd;
N=1e29/sum(sum(sum(den_atm)));


Volum=ones(hd,ftd,ffd);
oRi=Radii;
Ri=Radii;
for ht=1:hd
    Ri=Ri+unih*(1+alpha2)^(ht-1);
    oth=-pi/2;
    for td=1:ftd
        if(td~=1)
        oth=thi;
        end
        thi=td/ftd*pi-pi/2;
        Volum(ht,td,:)=abs((oRi^3+Ri^3))/6*dfi*abs(sin(oth)-sin(thi));
    end
    oRi=Ri;
end
%for =1:hd
%    for td = ceil(ftd/2) : ftd
%Volum(ht,td,:)=Volum(ht,ftd-td+1,:);
%    end
%end
copy_atm=den_atm;
for ht=1:hd
    for td=1:ftd
   copy_atm(ht,td,:) = copy_atm(ht,td,:).*N./Volum(ht,td,:);
    end
end
% 2

copy_atm3=den_atm;

copy_atm2=den_atm;
% Area 
dthi = pi/ftd;
Area = zeros(ftd,1);
for td = 1 : ftd
    othi = td/ftd*pi-pi/2-dthi;
    thi = td/ftd*pi-pi/2;
    dx = Radii * dthi;
    oRR = Radii * abs(sin(othi));
    RR = Radii * abs(sin(thi));
    Area(td) = (oRR+RR)/2*dx;
end


f = prn/totaln;
vth = zeros(ftd,ffd);
% thermal speed
for td = 1 : ftd
    for fd = 1 : ffd
        TTT = plotT(td,fd);
        vth(td,fd) = sqrt(2*k*TTT/mmass*mole);
    end
end


for ht = 1 : hd
    for td = 1 : ftd
        for fd =1  :ffd
            copy_atm2(ht,td,fd) = f * copy_atm2(ht,td,fd) / Area(td) / vth(td,fd);
        end
    end
end


% copy_atm=copy_atm2;

%% Running average

add_td = 1 ;
copy_atm3=copy_atm;
for td = 1 : ftd
    
    for ht = 1 : hd
        
        for fd = 1 : ffd
        mean_value = 0;    
        for adtd = td-add_td : td+add_td
            fd2 = fd;
            td2 = adtd;
            if(td2<=0)
                fd2=fd+ffd/2;
                td2=abs(adtd)+1;
                if(fd2>ffd),fd2=fd2-ffd;end
            end
            if(td2>ftd)
                fd2=fd+ffd/2;
                td2=ftd-(adtd-ftd)+1;
                if(fd2>ffd),fd2=fd2-ffd;end
            end
            
            mean_value = mean_value+copy_atm(ht,td2,fd2);
        end
            copy_atm3(ht,td,fd) = mean_value;
            
        end
    end
    end

copy_atm=copy_atm3;
%%
figure
plot_den=zeros(hd,ftd);
add_fd = 2;
plot_den(:,:)=copy_atm(:,:,ceil(3*ffd/4-add_fd+1));
plotT1=plot_den;
for ht=1:hd
    for td=1:ftd
        if(plotT1(ht,td)<=0),plotT12(ht,td)=0;
        else plotT12(ht,td)=log10(plotT1(ht,td));
        end        
    end
end
imagesc(plotT12);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
RR=(max_height-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 24:00     ' 'Number desity by log_1_0(N) in [#^.m^-^3]'];
title(aa);
%% 
figure
plot_den=zeros(ht,ftd);

plot_den(:,:)=copy_atm(:,:,ceil(3*ffd/4));
plotT13=plot_den;

for ht=1:hd
    for td=1:ftd
        if(plot_den(ht,td)==0),plotT13(ht,td) = 0;
            else
            plotT13(ht, td) = log10(plot_den(ht, td));
        end
    end
end


imagesc(plotT13);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
RR=(max_height-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 12:00     ' 'Number desity by log_1_0(N) in [#^.m^-^3]'];
title(aa);

%% Nightside

figure(7)
%subplot(1,2,1)
plot_den=zeros(hd,ftd);
plot_den2=zeros(hd,ftd);
add_fd=2;
for fd_ii = ffd*3/4-add_fd : ffd*3/4+add_fd-1
plot_den(:,:) = plot_den(:,:) + copy_atm(:,:,fd_ii);
plot_den2(:,:) = plot_den2(:,:) + copy_atm3(:,:,fd_ii);
end
plot_den(:,:) = plot_den(:,:) ./ 2/add_fd;
plot_den2(:,:) = plot_den2(:,:) ./ 2/add_fd;

plotT1=plot_den;
for ht=1:hd
    for td=1:ftd
        if(plotT1(ht,td)<=0),plotT12(ht,td)=0;
        else plotT12(ht,td)=log10(plotT1(ht,td));
        end        
    end
end
imagesc(plotT12);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
RR=(max_height-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 24:00 '];
title(aa);
hold on;
%% Mean
%subplot(1,2,2)
figure;
plot_den=zeros(hd,ftd);
plot_den3=zeros(hd,ftd);

for fd_ii = ffd/4-add_fd : ffd/4+add_fd-1
plot_den(:,:) = plot_den(:,:) + copy_atm(:,:,fd_ii);
plot_den3(:,:) = plot_den3(:,:) + copy_atm3(:,:,fd_ii);
end
plot_den(:,:) = plot_den(:,:) ./ 2./add_fd;
plot_den3(:,:) = plot_den3(:,:) ./ 2./add_fd;

plotT1=plot_den;
for ht=1:hd
    for td=1:ftd
        if(plot_den(ht,td)<=0),plotT1(ht,td)=0;
        else plotT1(ht,td)=log10(plotT1(ht,td));
        end        
    end
end
imagesc(plotT1);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
RR=(max_height-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 12:00   '];
title(aa);
hold on;
%%
figure
td=ftd/2;
thi=ceil(((td/ftd)-0.5)*180);
plot_den=zeros(hd,ffd);


for td_ii = ftd/2-1 : ftd/2 +1
    for ht = 1 : hd
        for fd = 1 : ffd
plot_den(ht,fd) = plot_den(ht,fd) + copy_atm(ht,td_ii,fd);
%plot_den3(:,:) = plot_den3(:,:) + copy_atm3(:,:,fd_ii);
        end
    end
end

plotT2=plot_den./3;

for ht=1:hd
    for fd=1:ffd
        if(plotT2(ht,fd)<=0),plotT2(ht,fd)=0;
        else plotT2(ht,fd)=log10(plotT2(ht,fd));
        end        
    end
end


copyT2=plotT2;

%for tt=1:hd
%    for fd=1:ffd
%        if(plot_den(tt,fd)==0),copyT2(tt,fd)=NaN;end
%    end
%end
plot_den(:,:)=copyT2;

imagesc(plot_den);
plotT2=plot_den;
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('Local time');
ylabel('height [km]');
aaa=['latitude=' num2str(thi) '.0^o' '    Number desity by log_1_0(N)'];
title(aaa)
%% H2O cap
Area=zeros(ftd,ffd);
thi=-pi/2;
for td=1:ftd
    oth=thi;
    thi=td/ftd*pi-pi/2;
    Area(td,:)=Radii^2*dfi*abs(cos(thi)-cos(oth));
end

figure(3)
copy=zeros(ftd,ffd);
copy2=horizon;
for td=1:ftd
    copy2(td,:)=horizon(td,:)./Area(td,:).*N;
end
for td=1:ftd
    copy(td,:)=copy2(ftd-td+1,:);
end
% copy(:,:)=log10(copy(:,:));
imagesc(copy);
plotT3=copy2;

hold on;
colorbar
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint2);
set(gca,'yticklabel',{'Northern Cap','0','Southern Cap'});
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'xtick',gcapoint1);
set(gca,'xticklabel',{'6','12','18','24','6'});
xlabel('Local Time');
ylabel('latitude');
title('horizontal CO2')
%% write it on file
fid=fopen('den_0417_1e4_col1e24.dat','w');
fprintf(fid,'%d\n',ftd);
fprintf(fid,'%d\n',ffd);
fprintf(fid,'%d\n',hd);
for tt=1:hd
for td=1:ftd
    for fd=1:ffd
        fprintf(fid,'%9.6f\n',den_atm(tt,td,fd));
    end
end
end
fclose(fid);
%% write it on file
fid=fopen('horizon_1e5_right_0125.dat','w');
fprintf(fid,'%d\n',ftd);
fprintf(fid,'%d\n',ffd);
for td=1:ftd
    for fd=1:ffd
        fprintf(fid,'%9.6f\n',horizon(td,fd));
    end
end
fclose(fid);
%% Write the ttime
fid=fopen('ttime_1e4_1e24_0417.dat','w');
fprintf(fid,'ouou is %d, total num is %d\n',ouou,totaln);
for nn=1:9
    for n=1:totaln
        fprintf(fid,'%9.5f\n',ttime(nn,n));
    end
end
%% 3D model again
unithi=5/180*pi;
hold on;
figure(12)
for ht = 1:hd
            Ri=unih*ht+Radii;
    for td=1:ftd
            thi=(td-ftd/2)/ftd*pi;
    %    for fd=1:ffd
            fd=ffd/4;
            fi=fd/ffd*2*pi;
            
            if(copy_atm(ht,td,fd)<=6),continue,end
            for n=1:copy_atm(ht,td,fd)
                rr=Ri+unih*rand;
                th=thi+(rand-0.5)*unithi;
                ff=fi+(rand-0.5)*unithi;
                plot(rr*cos(th),rr*sin(th),'r');
                hold on;
            end
    %    end
    end
    disp(ht)
end
max(max(max(copy_atm)))

%% read velocity distribution

fid = fopen('v_record_6.dat','r');
num_record = 0;

line = fgetl(fid);
    str = sscanf(line, '%d ');
td_list=zeros(str, 1);
fd_list=zeros(str, 1);
ht_list=zeros(str, 1);
v_record=zeros(str, 3);
    fprintf('%d\n',str);
while ~feof(fid)
    num_record = num_record + 1;
    line = fgetl(fid);
    str = sscanf(line, '%f %f %f ');
    v_record(num_record, 1 : 3) = str(1:3);
    
    line = fgetl(fid);
    str = sscanf(line, '%d ');
    td_list(num_record) = str;
    
    line = fgetl(fid);
    str = sscanf(line, '%d ');
    fd_list(num_record) = str;
    
    line = fgetl(fid);
    str = sscanf(line, '%d');
    ht_list(num_record) = str;
    if(num_record/1e5==ceil(num_record/1e5)),
        fprintf('%d\n',num_record/1e5),
    end
end

%% Arrange the v_record


num_record = num_record - 1;
new_v_record = zeros(ftd, ffd, hd, 50);
count_num= zeros(ftd,ffd,hd);
fprintf('Estimeted --> Elapsed time is 15.364213 seconds.\n');
for num_list = 1 : num_record
    td = td_list(num_list);
    fd = fd_list(num_list);
    ht = ht_list(num_list);
    count_num(td, fd, ht) = count_num(td, fd, ht) + 1;
    new_v_record(td, fd, ht, count_num(td, fd, ht)) = sqrt(sum(v_record(num_list, :).^2));
    if(num_list/1e5==ceil(num_list/1e5)),disp(num_list),end
end

fprintf('%d\n',sum(sum(sum(count_num))));
fprintf('%d\n',num_record);

%% parameters
clear v_record td_list fd_list ht_list Volum
%% New find velocity distribution using v_record
tic
num_poly = 3;
poly_mat = zeros(ftd,ffd,hd,num_poly+1);
min_v_mat = zeros(ftd,ffd,hd);
max_v_mat = zeros(ftd,ffd,hd);
check_count = zeros(ftd,ffd,hd);
for td = 1 : ftd
    for fd = 1 : ffd
        for ht = 1 : hd
            if(count_num(td,fd,ht) <= 8),
                check_count(td,fd,ht) = 1;
                continue;
            end
            v_dist = zeros(count_num(td, fd, ht), 1);
            nn = count_num(td,fd,ht);
            v_dist(1:nn) =  new_v_record(td, fd, ht, 1:nn);
             min_v = min(v_dist(:));
             max_v = max(v_dist(:));
            for ii = 1 : nn
                for jj = 1 : nn
                    if(v_dist(jj) > v_dist(ii))
                        copy = v_dist(jj);
                        v_dist(jj) = v_dist(ii);
                        v_dist(ii) = copy;
                        
                    end
                end
            end
            yy = hist(v_dist, ceil(nn/2));
            xx = zeros(ceil(nn/2), 1);
            xx(:) = linspace(min_v, max_v, ceil(nn/2));
            zz = zeros(nn, 1);
            zz(:) = 1 : nn;
            yy2=zeros(ceil(nn/2), 1);
            for num_yy = 1 : ceil(nn/2)
                yy2(num_yy) = sum(yy(num_yy:ceil(nn/2)));
            end
            plot(xx, yy2, 'kp');hold on;
            %p= polyfit(zz, v_dist(zz), num_poly);
            %plot(zz,v_dist(zz),'rp');hold on;
            %ff1 = polyval(p, zz);
            %plot(ff1, zz)
            p2= polyfit(xx, yy2, num_poly);
            ff = polyval(p2, xx);
            %if(min(size(find(ff<0)))~=0),
            %num_od_ff = find(ff<0);
            %if(num_od_ff~=ceil(nn/2)),disp('!!!!!!');end
            %fprintf('%d %d %d\n', td, fd, ht);break;
            %end
            plot(xx,ff);hold on;
            min_v_mat(td,fd,ht) = min_v;
            max_v_mat(td,fd,ht) = max_v;
            poly_mat(td,fd,ht,1:num_poly+1) = p2(1:num_poly+1);
        end
    end
    disp(td)
end
toc



%% velocity matrix for different T(0226)!!!!
nnn=1;
num=zeros(23,1);
p_max=zeros(23,1e3);
for MIV_T=linspace(20,240,23)
    n_v=0;
for vv=linspace(0,1000,101)
aa=4*pi*vv^2*sqrt((1/2/pi*mmass/mole/k/MIV_T)^3)*exp(-mmass/mole*3*vv^2/2/k/MIV_T);
aa=aa*5e5;
if(floor(aa)==0),continue;end
p_max(nnn,n_v+1:n_v+1+floor(aa))=vv;
n_v=n_v+floor(aa);
end
n_v=n_v+1;
num(nnn)=n_v;
nnn=nnn+1;
end
disp(num)
num2=num;
%% check
ii = 23;
ppp = p_max(ii, 1 : num(ii));

x = 10 : 10 : 1000; figure;

hist(ppp,x)
xlabel('thermal velocity [m/s]');
ylabel('test numbers');

hold on

%% run 2 setting first

totaln = 1e4;


h=0.5;
num_Bong = 0;
den_atm = zeros(hd, ftd, ffd);
max_nn = 0;
max_rot = 0;
record_nn = 1;
ouou = 0;
cross_section = 1e-19;
v2 = zeros(3,1);
bomb_line = zeros(totaln, 1);
ttime = zeros(9, totaln);
del_v=zeros(3,1);
vcm=del_v;

%%
%plotT=zeros(ftd,ffd);
%figure
%sphere;
%hold on;
for n = 1 : totaln
    rest_time = 0;
   del_s=0;
   num_col = 0;
   real_col_p = 0;
   col_p = rand;
   ou=0;
   rot=0;
   t_n=0;
%   record_io=0;
   
   time=0;
   disp(n)
   W=1;
   solar_time=0;
   thi=rand*pi-pi/2;
   ttime(1,n) = thi;
   
   fi=rand*pi*2;
   ttime(6,n) = fi;
   flight_time=0;
   s(:)=[Radii*cos(thi)*sin(fi) Radii*cos(thi)*cos(fi) Radii*sin(thi)];   
   %traj=zeros(1,3,5e4);
   del_time=0;
 %  vout_first=p_Max(ceil(rand*n_v1));
    v_io=0;
    td=ceil((thi+pi/2)/pi*ftd);
    fd=ceil(fi/2/pi*ffd);
    ttime(7,n) = plotT(td,fd);
    while(W>=0.01)
    KF=2*pi*Radii*abs(cos(thi))/P;
    td=ceil((thi+pi/2)/pi*ftd);
    fd=ceil(fi/2/pi*ffd);
        if(fd==0),fd=ffd;end
    local_T=plotT(td,fd);
        if(local_T<100),
            horizon(td,fd) = horizon(td,fd) +W;
            for add=1:ffd
            rot=rot+1;
            
            time=time+dt;
            rest_time = rest_time + dt;
            %disp('1');
            %pause(0.1);
            %sfd=fd-add;
            %if(sfd<=0),sfd=sfd+ffd;end
            %fi=sfd/ffd*pi*2;
            %s(:)=[Radii*cos(thi)*sin(fi) Radii*cos(thi)*cos(fi) Radii*sin(thi)];
            %plot3(s(1),s(2),s(3),'rp');pause(0.1);
            %if(abs(s(1)-os(1))>=1e4),disp(s(1));end
            %os=s;
            %hold on;
           %t_n=t_n+1;
           %traj(1,:,t_n)=s(:);
            %if(nt>ffd),nt=nt-ffd;end
            if(plotT(td,add)>=100),
                %if(fd<=0),fd=fd+ffd;end
                rot_fi = (add-fd)/ffd*2*pi;
                s(:)=[s(1)*cos(rot_fi)-sin(rot_fi)*s(2) sin(rot_fi)*s(1)+cos(rot_fi)*s(2) s(3)];
                fd = add;
                fi=fd/ffd*pi*2;
                break;
            end
            end
            if(max(plotT(td,:))<100),horizon(td,fd)=horizon(td,fd)+W;break;end
        end
   flyd=rand*pi+thi-pi/2;
   xi=rand*pi+fi-pi/2;
   %if(abs(flyd)<=0.1),continue;end  ??
   v(1:3)=0;
   local_T = plotT(td,fd);
    %if(v_io==0),v_io=1;vout=vout_first;
    %else
        if(local_T>=240),vr=p_max(23,ceil(rand*num(21)));
        else
        vr=p_max(ceil(local_T/10)-1,ceil(rand*num(ceil(local_T/10)-1)));
        end
        vout=sqrt(alph*vr^2+(1-alph)*sum(v(:).^2));
    %end
   v(1)=cos(flyd)*sin(xi);
   v(2)=cos(flyd)*cos(xi);
   v(3)=sin(flyd);
   v(:)=v(:).*vout;
   v(1)=v(1)+KF*sin(fi+pi/2);  %% velocity induced by rotational force 
   v(2)=v(2)+KF*cos(fi+pi/2);
    if(sqrt(sum(v(:).^2))>=ve),ou=1;break;end
   Ri = sqrt(sum(s(:).^2));%%
   oa = G*M*s(:)/(Ri^3);
   %disp('!')
    for ii=1:1e7
        %disp('2')
        flight_time=flight_time+h;
        time=time+h;
        %if(time>=3e3),break;end
        %pause(0.1)
        for i=1:3
            ns(i)=s(i)+v(i)*h-G*M*s(i)/(Ri^3)*(1.5*(h^2))-oa(i)*(h^2)/6;
            v(i)=v(i)-G*M*s(i)/(Ri^3)*h;
            oa(i)=G*M*s(i)/(Ri^3);
            now_del_s(i)=s(i)-ns(i);
        end
        
        s(:)=ns(:);
        Ri = sqrt(sum(s(1:3).^2));
        if(Ri <= Radii),break;end
        if(Ri >= max_height),ou=1;break;end
       
        %t_n=t_n+1;traj(1,:,t_n)=s(:);
        
        thi=asin(s(3)/Ri);
        if(s(1)==0 && s(2)>0), fi = 0 ;
        elseif(s(1)==0 && s(2)<0), fi = pi;
        else
            if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
            else fi=-atan(s(2)/s(1))+pi/2;
            end
        end
        fd=ceil(fi/2/pi*ffd);
        td=ceil((thi+pi/2)/pi*ftd);
        if(td==0),td=1;end
   
        if( fd <= ffd/2),
            W = W * exp(-h/lifetime);
            solar_time = solar_time+h;
        end
        %if(ceil(ii/100)~=ii/100),continue;end
    %if(abs(s(1)-os(1))>=1e4),disp(s(1)*1e4);end
    %os=s;
    
    %if(ii/5e1==ceil(ii/5e1)),plot3(Ri/Radii*cos(thi)*sin(fi),Ri/Radii*cos(thi)*cos(fi),Ri/Radii*sin(thi),'kp');hold on;end
    %if(ii/1e3==ceil(ii/1e3)),pause(0.1);end
    
    td=ceil((thi+pi/2)/pi*ftd);
    fd=ceil(fi/2/pi*ffd);
    if(ii==1),continue;end
    ht = ceil(log(alpha2*(Ri-Radii)/unih+1)/log(1+alpha2));
    %---------------- collision ----------------
    real_col_p = real_col_p + sqrt(sum(now_del_s(:).^2)) * copy_atm(ht, td, fd) * cross_section;
    
    if(col_p < 1 - exp(-real_col_p) && count_num(td, fd, ht) > 8),
        real_col_p = 0;
        disp('*');
        max_v = max_v_mat(td, fd, ht);
        min_v = min_v_mat(td, fd, ht);
        nn = count_num(td,fd,ht);
        num_col = num_col + 1;
        p2(1:num_poly+1) = poly_mat(td,fd,ht,1:num_poly+1);
        %xx1 = linspace(min_v, max_v, ceil(nn/2));
        %vv2 = polyval(p2, xx1);
        df = polyval(p2, max_v);
        ranking = rand * (nn-df) + df;
        f1 =  polyval(p2, max_v) - ranking;
        v1 = max_v;
        f0 = polyval(p2, min_v) - ranking;
        v0 = min_v;
                for num_dev = 1 : 1e2
                    vv2 = (v1 + v0)/2;
                    f2 = polyval(p2, vv2) - ranking;
                    if(f2 * f1 < 0),f0 = f2; v0 = vv2;
                    else f1 = f2; v1 = vv2;
                    end
                    if(abs(v1-v0) < 0.01)
                        break;
                    end
                end
                
           % plot(xx1, ones(size(xx1)).*ranking,'r-');
        %plot(xx1, vv2, 'kp');hold on;
        %disp('Bong!');
        num_Bong = num_Bong + 1;
        thi2 = pi * rand - pi / 2;
        fi2 = 2 * pi * rand;
        v2(:) = vv2 .* [cos(thi2)*cos(fi2) cos(thi2)*sin(fi2) sin(thi2)];
        for rp_vn = 1 : 3
        del_v(rp_vn) = v(rp_vn) - v2(rp_vn);
        vcm(rp_vn) = v(rp_vn) + v2(rp_vn);
        end
        for rp_vn = 1 : 3
            v(rp_vn) = vcm(rp_vn) + 0.5 * del_v(rp_vn);
        end
        continue
    end
    
    %-------------------------------------------
        %if(ht>=hd),ou=1;disp(ht);break;end
    
    %if(Ri <= max_height/3 && max_nn~=td+fd+ht && count_nn(td, fd, ht) <= 200),
    %    v_record(record_nn, 1:3) = v(1:3);
    %    v_record(record_nn, 4:6) = [td fd ht];
    %    record_nn = record_nn+1;
    %    count_nn(td, fd, ht) = count_nn(td, fd, ht) + 1;
    %    max_nn = td+fd+ht;
    %end
    
        den_atm(ht,td,fd)=den_atm(ht,td,fd)+W;
    end
    
    thi=asin(s(3)/Ri);
        if(s(1)==0 && s(2)>0), fi = 0 ;
        elseif(s(1)==0 && s(2)<0), fi = pi;
        else
            if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
            else fi=-atan(s(2)/s(1))+pi/2;
            end
        end
        fd=ceil(fi/2/pi*ffd);
        td=ceil((thi+pi/2)/pi*ftd);
        if(td==0),td=1;end
   
    if(ou==1),break;end
    del_time=del_time+h*ii;
    if(td<=2||td>=35),horizon(td,fd)=horizon(td,fd)+0.05*W;W=W*0.95;end
    if(del_time>=dt),del_n=floor(del_time/dt);rot=rot+del_n;fd=fd-del_n;if(fd<=0),fd=fd+ffd;end,del_time=del_time-dt*del_n;continue;end
   %if(time>=1e4),ou=1;break;end
    end
   
   bomb_line(n) = num_col; 
   ttime(3,n) = time;   
   ttime(4,n) = W;
   ttime(5,n) = flight_time;
   ttime(8,n) = rest_time;
   ttime(9,n) = bomb_line(n);

   if(ou==1),ouou=ouou+1;continue;
   %else break;
   end
   if(rot>=max_rot),max_rot=rot;end
ttime(2,n)=thi;

   %endingpoint(1)=t_n;
%   disp(time)
  
end
if(n~=totaln)
bomb_line(n:totaln) = [];
end
%% check
td = 1;
fd = 1;
ht = 1;
str = zeros(3, 1);
figure
for num_check = 1 : count_num(td,fd,ht)
str(1)=new_v_record(td,fd,ht,num_check,1);
plot(str(1)),hold on;
end

%%
sum_atm = 0;
for td = 1 : ftd
    
    for fd = 1 : ffd
        
        for ht = 1 : hd
            sum_atm = sum_atm + Volum(ht,td,fd) * copy_atm(ht,td,fd);
        end
    end
end

%% Get the histogram of bomb numbers

figure

num_hist_final = 100;

copy_x = bomb_line;
copy_x(find(bomb_line==0))=[];
x = linspace(0, max(bomb_line), num_hist_final);

hist(copy_x, x) 
xlabel('number of collision happening');
ylabel('histogram');
%%
figure;

x = linspace(0, max(bomb_line), num_hist_final);
num_w_hist = zeros(size(x));
time_hist = zeros(size(x));
thi_hist =  zeros(size(x));
ii=7;
for num_hist = 1 : 1e3
    
    nn = ceil(bomb_line(num_hist)/max(bomb_line)*num_hist);
    if(nn==0),continue;end
    if(ttime(3,num_hist) ~= 0)
    %time_hist(nn) = time_hist(nn) + ttime(3, num_hist) - ttime(ii, num_hist);
    thi_hist(nn) = thi_hist(nn) + ttime(ii, num_hist);
    plot(nn, ttime(ii, num_hist),'p');hold on;
    num_w_hist(nn) = num_w_hist(nn) + 1;
    end
end
for num_hist = 1 : num_hist_final
    if(num_w_hist(num_hist)~=0)
        time_hist(num_hist) = time_hist(num_hist)./num_w_hist(num_hist);
        
        thi_hist(num_hist) = thi_hist(num_hist)./num_w_hist(num_hist);
    end
end
for num_hist = 1 : num_hist_final
x2(num_hist) = x(num_hist) - 2;
end
%%
plot(x(:),time_hist(:),'r')
legend('total time','flight time');

%%
new_v_record2 = zeros(ftd, hd);
for td  =  1 : ftd
    for hd = 1 : ht
for ii = 1 : 200
        new_v_record2(td, hd)=new_v_record2(td, hd) + new_v_record(td, hd, ii);
end
    end
end
new_v_record2(:,:) = new_v_record2(:,:)./80567;

figure;imagesc(new_v_record2)


%% Dayside & Nightside histogram

%n = n - 1;
T_line = zeros(n, 1);
region_line = zeros(n, 1);
T_line(:) = ttime(7, 1:n);
region_1 = zeros(n, 1);
region_2 = zeros(n, 1);
aa = 0;
bb = 0;
judge_T = 100;
totaln_bomb = 0;
for nn = 1 : n
    if(bomb_line(nn)==0),
        totaln_bomb = totaln_bomb + 1;
        continue;end
    if(T_line(nn)>judge_T)
        aa = aa + 1;
        region_line(nn) = 1;
        region_1(aa) = nn;
    else
        bb = bb + 1;
        region_line(nn) = 2;
        region_2(bb) = nn;
        
    end 
end
region_1(aa + 1 : n) = [];
region_2(bb + 1 : n) = [];

copy_bomb = bomb_line;
bomb1=zeros(aa, 1);
for nn = 1 : aa
bomb1(nn) = copy_bomb(region_1(nn));
end
figure;

subplot(1, 2, 1);
str1_1 = ['Day side(T>' num2str(judge_T) 'K)'];

hist(bomb1, 20)
title(str1_1)

%
bomb2=zeros(bb, 1);
for nn = 1 : bb
bomb2(nn) = copy_bomb(region_2(nn));
end

subplot(1, 2, 2);
str1_2 = ['Night side(T<' num2str(judge_T) 'K)'];

hist(bomb2, 20)
title(str1_2)

%% plot the minimum temperature
imagesc(plotT)
hold on;
for thi = linspace(1, ftd, 7)
    td = ceil(thi);
    min_T = min(plotT(td,:));
    plot(linspace(1,ffd,1000),linspace(td,td,1000));hold on;
    text(td/ftd*ffd,td,num2str(min_T),'FontSize',18)
end
%% contour
%contour(plotT)
plotT2=zeros(ftd,ffd);
for td  = 1 : ftd
    for fd = 1 : ffd
        if(plotT(td,fd) > 100),
            plotT2(td,fd) = 1;
        else plotT2(td,fd) = 0;
        end
    end
end
imagesc(plotT2)
%% test
ss = [1 2 3 ];
rot_fi = 0.03;
for ii = 1 : 40
ss(:)=[ss(1)*cos(rot_fi)-sin(rot_fi)*ss(2) sin(rot_fi)*ss(1)+cos(rot_fi)*ss(2) ss(3)];
plot3(ss(1),ss(2),ss(3),'rp');hold on;
end



%%
G = 6.67e-11;
n_mat = 1:5;
M_mat = zeros(3, 2);
M_mat(1,:)=[2 16];
M_mat(2,:)=[0.8 1.2];
M_mat(3,:)=[0.075 0.5];
T_mat(1,:)=[ ];
T_mat(2,:)=[ ];
T_mat(3,:)=[ ];

figure(2)

for num = 1 : 3;
    subplot(3,3,num);
R=1;M=1;n=1;
Pc2=((4*pi*G)^(1/n))/(n+1)*(G*M/Mn)^((n-1)/n)*(R/Rn)^((3-n)/n)*(M/R^3)^((n+1)/n);
figure(2)
for nc = 1 : 100
ran  = rand;
M=M_mat(num,1).*ran+(num-ran)*M_mat(num,2);

R = 1;
Mn = 1 ;
Rn = 1 ;
ii=0;
for n = 1 : 0.05 : 5
ii=ii+1;

Pc(ii) = 1/Pc2*((4*pi*G)^(1/n))/(n+1)*(G*M/Mn)^((n-1)/n)*(R/Rn)^((3-n)/n)*(M/R^3)^((n+1)/n);
end
hold on;
plot((1:0.05:5),Pc(:),'g')
end
end
%%
%xx = linspace(1, ftd, ftd);
%yy = ones(ftd, 1) .* (scale_height/unih);
%plot(xx, yy, 'k-');
%hold on;
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:')

exp_line = zeros(ftd, 1);
for td = 1 : ftd
    den_most = plot_den(1, td);
    for ht = 1 : hd
        if(plot_den(ht, td) <= 1e10) %den_most * exp(-1)
            exp_line(td) = ht;
            break;
        end
    end
end
plot(xx, exp_line,'LineWidth',2);
hold all;

exp_line = zeros(ftd, 1);
for td = 1 : ftd
    den_most = plot_den(1, td);
    for ht = 1 : hd
        if(plot_den(ht, td) <= 1e8) %den_most * exp(-1)
            exp_line(td) = ht;
            break;
        end
    end
end
plot(xx, exp_line,'LineWidth',2);
hold all;


exp_line = zeros(ftd, 1);
for td = 1 : ftd
    den_most = plot_den(1, td);
    for ht = 1 : hd
        if(plot_den(ht, td) <= 1e7) %den_most * exp(-1)
            exp_line(td) = ht;
            break;
        end
    end
end
plot(xx, exp_line,'LineWidth',2);
hold all;

legend(3,'Scale height','1e10','1e8','1e7');
