clear all;
%% Read the atmospheric density
fid=fopen('den_0211.dat','r');
for ii=1:2
    line = fgetl(fid);
    %str=sscanf(line,'%d');
    str=sscanf(line, '%f');
    if(ii==1),ay=str;end
    if(ii==2),by=str;end
end
ftd=ay;
ffd=by;
den_atm=zeros(100,ftd,ffd);
for time=1:100
for ii=1:ftd
    for jj=1:ffd
    line = fgetl(fid);
    str=sscanf(line, '%f');
    den_atm(time,ii,jj)=str;
    end
end
end
fclose(fid);
copy_atm=den_atm;

%% Read the temperature
fid=fopen('time_temperature_VP.dat','r');
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
for time=1:ffd
for ii=1:ftd
    for jj=1:ffd
    line = fgetl(fid);
    str=sscanf(line, '%f');
    copyT(time,ii,jj)=str;
    end
end
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
hd=50;   %% height grid points
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
%% collision parameters

v_record = zeros(1e5, 6);
count_nn = zeros(ftd, hd);

%% Calculate just for the time scale histogram
solar = zeros(ffd,3);
cal = zeros(1,3);
horizon = zeros(ftd,ffd);
hd=100;   %% height grid points
unih = Radii*0.1;
maxheight = unih * hd + Radii;  %% [meter]
den_atm = zeros(hd, ftd, ffd);
za = 12.3/180*pi;%%
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
%% 
max_nn = 0;
max_rot = 0;
totaln = 1e4;
record_nn = 1;
ouou = 0;
%%
for n=1:totaln
   ou=0;
   rot=0;
   t_n=0;
%   record_io=0;
   nt=ceil(rand*ffd);if(nt==0),nt=nt+ffd;end
   time=0;
   disp(n)
   W=1;
   solar_time=0;
   thi=rand*pi-pi/2;
   ttime(1,n)=thi;
   fi=rand*pi*2;
   flight_time=0;
   s(:)=[Radii*cos(thi)*sin(fi) Radii*cos(thi)*cos(fi) Radii*sin(thi)];   
   %traj=zeros(1,3,5e4);
   del_time=0;
 %  vout_first=p_Max(ceil(rand*n_v1));
   v_io=0;
    while(W>=0.01)
    KF=2*pi*Radii*abs(cos(thi))/P;
    td=ceil((thi+pi/2)/pi*ftd);
    fd=ceil(fi/2/pi*ffd);
        if(fd==0),fd=ffd;end,if(td==0),td=ftd;end
    local_T=copyT(nt,td,fd);
        if(local_T<100),
            for add=1:ffd
            rot=rot+1;
            nt=nt+1;
            time=time+dt;
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
            if(nt>ffd),nt=nt-ffd;end
            if(copyT(nt,td,fd)>=100),fd=fd-add;
                if(fd<=0),fd=fd+ffd;end
                fi=fd/ffd*pi*2;
                break;
            end
            end
            if(max(copyT(:,td,fd))<100),horizon(td,fd)=horizon(td,fd)+W;break;end
        end
   flyd=rand*pi+thi-pi/2;
   xi=rand*pi+fi-pi/2;
   %if(abs(flyd)<=0.1),continue;end  ??
   v(1:3)=0;
   local_T=copyT(nt,td,fd);
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
        if(Ri >= maxheight),ou=1;break;end
       
        %t_n=t_n+1;traj(1,:,t_n)=s(:);
        
   thi=asin(s(3)/Ri);
    if(s(1)==0), fi=0;fd=0;
    else
        if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
        else fi=-atan(s(2)/s(1))+pi/2;
        end
        fd=ceil(fi/2/pi*ffd);
    end
   td=ceil((thi+pi/2)/pi*ftd);
   if(td==0),td=1;end
        
        if(fd <= ffd/2),
            W = W * exp(-h/lifetime);
            solar_time = solar_time+h;
        end
        %if(ceil(ii/100)~=ii/100),continue;end
    %if(abs(s(1)-os(1))>=1e4),disp(s(1)*1e4);end
    %os=s;
    %plot3(Ri*cos(thi)*sin(fi),Ri*cos(thi)*cos(fi),Ri*sin(thi),'kp');hold on;pause(0.1);

    if(ii==1),continue;end
    ht=ceil((Ri-Radii)/unih);
        %if(ht>=hd),ou=1;disp(ht);break;end
    
    if(max_nn~=td+fd+ht && count_nn(td, ht) <= 200),
        v_record(record_nn, 1 : 3) = v(1:3);
        v_record(record_nn, 4 : 6) = [td fd ht];
        record_nn = record_nn+1;
        count_nn(td, ht) = count_nn(td, ht) + 1;
        max_nn = td+fd+ht;
    end
    if(record_nn>=100*ftd*hd),disp('!!!');end
        den_atm(ht,td,fd)=den_atm(ht,td,fd)+W;
    end
   thi=asin(s(3)/Ri);
    if(s(1)==0), fi=0;fd=0;
    else
        if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
        else fi=-atan(s(2)/s(1))+pi/2;
        end
        fd=ceil(fi/2/pi*ffd);
    end
   td=ceil((thi+pi/2)/pi*ftd);
   if(td==0),td=1;end
    if(ou==1),break;end
   del_time=del_time+h*ii;
    if(td<=2||td>=35),horizon(td,fd)=horizon(td,fd)+0.05*W;W=W*0.95;end
    if(del_time>=dt),del_n=floor(del_time/dt);rot=rot+del_n;fd=fd-del_n;if(fd<=0),fd=fd+ffd;end,del_time=del_time-dt*del_n;continue;end
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
   %endingpoint(1)=t_n;
%   disp(time)
  
end
%% write "v_record" on file

fid = fopen('v_record.dat','w');
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
dfi = 2 * pi / ffd;
dthi = pi / ftd;
N=1e29/sum(sum(sum(den_atm)));
Volum=ones(hd,ftd,ffd);
oRi=0;
Ri=Radii;
for ht=1:hd
    Ri=Ri+unih;
    oth=-pi/2;
    for td=1:ftd
        oth=thi;
        thi=td/ftd*pi-pi/2;
        Volum(ht,td,:)=abs((oRi^3+Ri^3))/6*dfi*abs(sin(oth)-sin(thi));
    end
    oRi=Ri;
end
copy_atm=den_atm;
for ht=1:hd
    for td=1:ftd
   copy_atm(ht,td,:) = copy_atm(ht,td,:).*N./Volum(ht,td,:);
    end
end
%%
figure
plot_den=zeros(hd,ftd);
plot_den(:,:)=copy_atm(:,:,ceil(3*ffd/4));
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
RR=(maxheight-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 12:00     ' 'Number desity by log_1_0(N) in [#^.m^-^3]'];
title(aa);
%% 
figure
plot_den=zeros(ht,ftd);
plot_den(:,:)=copy_atm(:,:,ceil(ffd/4));
plotT1=plot_den;
for ht=1:hd
    for td=1:ftd
        if(plotT1(ht,td)==0),plotT1(ht,td)=NaN;end
    end
end
plotT11=zeros(size(plotT1));
plotT11(:,:)=log10(plotT1(:,:));
imagesc(plotT11);
hold on;
colorbar
gcapoint1=[1 hd/2 hd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'ytick',gcapoint1);
set(gca,'xtick',gcapoint2);
set(gca,'xticklabel',{'-90','0','90'});
RR=(maxheight-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 12:00     ' 'Number desity by log_1_0(N) in [#^.m^-^3]'];
title(aa);

%% Mean
figure
plot_den=zeros(hd,ftd);
for fd_ii = 1 : ffd
plot_den(:,:) = plot_den(:,:) + copy_atm(:,:,fd_ii);
end
plot_den(:,:) = plot_den(:,:) ./ ffd;

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
RR=(maxheight-Radii)/1000;
RR2=RR/2;
str=num2str(RR);
str2=num2str(RR2);
set(gca,'yticklabel',{'horizon',str2,str});
xlabel('latitude');
ylabel('distance [km]');
aa=['Local Time = 12:00     ' 'Number desity by log_1_0(N) in [#^.m^-^3]'];
title(aa);
hold on;

xx = linspace(1, ftd, ftd);
yy = ones(ftd, 1) .* (scale_height/unih);
plot(xx, yy, 'k-');
hold on;
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
        if(plot_den(ht, td) <= 1e6) %den_most * exp(-1)
            exp_line(td) = ht;
            break;
        end
    end
end
plot(xx, exp_line,'LineWidth',2);
hold all;

legend(3,'Scale height','1e10','1e8','1e6');
%%
figure
td=ftd/2;
thi=ceil(((td/ftd)-0.5)*180);
plot_den=zeros(hd,ffd);
plot_den(:,:)=copy_atm(:,ceil(td),:);
copyT2=plot_den;
%for tt=1:hd
%    for fd=1:ffd
%        if(plot_den(tt,fd)==0),copyT2(tt,fd)=NaN;end
%    end
%end
plot_den(:,:)=log10(copyT2);

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
fid=fopen('den_0213.dat','w');
fprintf(fid,'%d\n',ftd);
fprintf(fid,'%d\n',ffd);
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
fid=fopen('ttime_1e3.dat','w');
for nn=1:5
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

fid = fopen('v_record.dat','r');
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
clear new_v_record

num_record = num_record - 1;
new_v_record = zeros(ftd, hd, 200);
count_num= zeros(ftd,hd);
fprintf('Estimeted --> Elapsed time is 15.364213 seconds.\n');
for num_list = 1 : num_record
    td = td_list(num_list);
    ht = ht_list(num_list);
    count_num(td, ht) = count_num(td, ht) + 1;
    new_v_record(td, ht, count_num(td, ht)) = sqrt(sum(v_record(num_list, :).^2));
    
end

fprintf('%d\n',sum(sum(sum(count_num))));
fprintf('%d\n',num_record);

%% parameters
clear v_record td_list fd_list ht_list Volum

% Scale height
T_scale = mean(mean(plotT));
g_value = G * M / (Radii ^ 2);

scale_height = k * T_scale / mmass * mole / g_value;


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
ttime = zeros(7, totaln);
%%
figure
for n =1: totaln
   del_s=0;
   num_col = 0;
   real_col_p = 0;
   col_p = rand;
   ou=0;
   rot=0;
   t_n=0;
%   record_io=0;
   nt=ceil(rand*ffd);if(nt==0),nt=nt+ffd;end
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
    ttime(7,n) = copyT(1,td,fd);
    while(W>=0.01)
    KF=2*pi*Radii*abs(cos(thi))/P;
    td=ceil((thi+pi/2)/pi*ftd);
    fd=ceil(fi/2/pi*ffd);
        if(fd==0),fd=ffd;end,if(td==0),td=ftd;end
    local_T=copyT(nt,td,fd);
        if(local_T<100),
            for add=1:ffd
            rot=rot+1;
            nt=nt+1;
            time=time+dt;
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
            if(nt>ffd),nt=nt-ffd;end
            if(copyT(nt,td,fd)>=100),fd=fd-add;
                if(fd<=0),fd=fd+ffd;end
                fi=fd/ffd*pi*2;
                break;
            end
            end
            if(max(copyT(:,td,fd))<100),horizon(td,fd)=horizon(td,fd)+W;break;end
        end
   flyd=rand*pi+thi-pi/2;
   xi=rand*pi+fi-pi/2;
   %if(abs(flyd)<=0.1),continue;end  ??
   v(1:3)=0;
   local_T=copyT(nt,td,fd);
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
    for ii=1:1e6
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
        if(Ri >= maxheight),ou=1;break;end
       
        
        thi=asin(s(3)/Ri);
            if(s(1)==0), fi=0;fd=0;
            else
                if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
                else fi=-atan(s(2)/s(1))+pi/2;
                end
            fd=ceil(fi/2/pi*ffd);
            end
        td=ceil((thi+pi/2)/pi*ftd);
        if(td==0),td=1;end
        if(fd <= ffd/2),
            W = W * exp(-h/lifetime);
            solar_time = solar_time+h;
        end
        %if(ceil(ii/100)~=ii/100),continue;end
    %if(abs(s(1)-os(1))>=1e4),disp(s(1)*1e4);end
    %os=s;
    
    %plot3(Ri*cos(thi)*sin(fi),Ri*cos(thi)*cos(fi),Ri*sin(thi),'kp');hold on;
    %if(ii/1e2==ceil(ii/1e2)),pause(0.1);end
    
    if(ii==1),continue;end
    ht = ceil((Ri-Radii)/unih);
    %---------------- collision ----------------
    real_col_p = real_col_p + sqrt(sum(now_del_s(:).^2)) * copy_atm(ht, td, fd) * cross_section;
    
    if(col_p < 1 - exp(-real_col_p)&& ht <= fix(hd/3) && count_num(td, ht)~=0),
        disp('*');
        del_s = 0;
        col_cand = randi(count_num(td, ht), 1);% collision candidatefrom the last time atmospheric velocity distribution
        num_col = num_col + 1;
        %disp('Bong!');
        num_Bong = num_Bong + 1;
        vv2 = new_v_record(td,ht,col_cand);
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
    
    %if(Ri <= maxheight/3 && max_nn~=td+fd+ht && count_nn(td, ht) <= 200),
    %    v_record(record_nn, 1:3) = v(1:3);
    %    v_record(record_nn, 4:6) = [td fd ht];
    %    record_nn = record_nn+1;
    %    count_nn(td, ht) = count_nn(td, ht) + 1;
    %    max_nn = td+fd+ht;
    %end
    
        den_atm(ht,td,fd)=den_atm(ht,td,fd)+W;
    end
    
   thi=asin(s(3)/Ri);
    if(s(1)==0), fi=0;fd=0;
    else
        if(s(1)<0), fi=-atan(s(2)/s(1))+3*pi/2;
        else fi=-atan(s(2)/s(1))+pi/2;
        end
        fd=ceil(fi/2/pi*ffd);
    end
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

num_hist_final = 3561;

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
new_v_record2 = zeros(ftd, ht);
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

n = n - 1;
T_line = zeros(n, 1);
region_line = zeros(n, 1);
T_line(:) = ttime(7, 1:n);
region_1 = zeros(n, 1);
region_2 = zeros(n, 1);
aa = 0;
bb = 0;
judge_T = 100;
for nn = 1 : n
    if(bomb_line(nn)==0),continue;end
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
