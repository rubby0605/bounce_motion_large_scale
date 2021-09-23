%% Meridian Plane
%figure;
hold on;
plotT11=log10(plot_den);
xe=zeros(hd,2*ftd);
ye=zeros(hd,2*ftd);
ze=zeros(hd,2*ftd);
fi=pi/2;
th=linspace(-pi/2,pi/2,ftd);
Ri=1;
for ht=1:hd
    td=0;
    Ri=Ri+unih/Radii*(1+alpha2)^(ht-1);
for thi=linspace(0,2*pi,2*ftd)
    td=td+1;
    xe(ht,td)=Ri*cos(thi);
    ze(ht,td)=Ri*sin(thi);

end
end
    ye(:,:)=0;
hold on;
cp_pT11=zeros(ht,ftd);
for td=1:ftd
    cp_pT11(:,ftd-td+1)=plotT1(:,td);
end
a1=plotT12(:,1:ftd/2);
a2=plotT12(:,ftd/2+1:ftd);
b1=cp_pT11(:,1:ftd/2);
b2=cp_pT11(:,ftd/2+1:ftd);
CC1=[b2 a1 a2 b1];

surf(xe,ye,ze,CC1)
%axis([-3.5 3.5 -3.5 3.5 -3.5 3.5 ])
%axis([-10 10 -10 10 -10 10  ])
shading interp
hidden off;
colormap(jet)
colorbar
set(gca,'XTick',linspace(-3e6/Radii,3e6/Radii,7),...
    'XTickLabel',{'-300km','-200km','-100km','0km','100km','200km','300km'})
set(gca,'ZTick',linspace(-3e6/Radii,3e6/Radii,7),...
    'ZTickLabel',{'-300km','-200km','-100km','0km','100km','200km','300km'},...
    'FontSize',16)
Y_Label=8:13;
set(gca,'YTickLabel',Y_Label)
hei = 3e6/Radii;
axis([-hei hei -hei hei -hei hei]);
caxis([7 12])
colorbar('YTick',[7:12],'YTickLabel',...
    {'10^7','10^8','10^9','10^10','10^11',...
     '10^12'},...
     'FontSize',16)
title('Number density of exosphere of Ceres [molecules-m^-^2-s^-^1]')
 
%colorbar('YTick',yyy);
%% Meridian Plane

figure;
hold on;
%plotT11=log10(plot_den);
%plotT12 = plot_den2;
%plotT11 = plot_den3;
xe=zeros(hd,2*ftd);
ye=zeros(hd,2*ftd);
ze=zeros(hd,2*ftd);
fi=pi/2;
th=linspace(-pi/2,pi/2,ftd);
Ri=1;
for ht=1:hd
    td=0;
    Ri=Ri+0.1;
for thi=linspace(0,2*pi,2*ftd)
    td=td+1;
    xe(ht,td)=Ri*cos(thi);
    ze(ht,td)=Ri*-sin(thi);

end
end
    ye(:,:)=0;
hold on;
cp_pT11=zeros(ht,ftd);
for td=1:ftd
    cp_pT11(:,ftd-td+1)=plotT11(:,td);
end
a1=plotT12(:,1:ftd/2);
a2=plotT12(:,ftd/2+1:ftd);
b1=cp_pT11(:,1:ftd/2);
b2=cp_pT11(:,ftd/2+1:ftd);
CC1=[b2 a1 a2 b1];

surf(xe,ye,ze,CC1)
%axis([-3.5 3.5 -3.5 3.5 -3.5 3.5 ])
%axis([-10 10 -10 10 -10 10  ])
shading interp
hidden off;
colormap hot
colorbar

hei = 3e6/Radii;
axis([-hei hei -hei hei -hei hei]);
for thi=linspace(0,pi*2,100)
rr=2;plot3(rr*cos(thi),-0.1,rr*sin(thi),'.k'),hold on;
end


%% for test==>equatorial [plotT2=>hd,ffd]
hold on;
xe=zeros(hd,ffd);
ye=zeros(hd,ffd);
ze=zeros(hd,ffd);
uni=unih/Radii;
ff=linspace(0,2*pi,ffd);
Ri=1;
for ht=1:hd
    Ri=Ri+unih/Radii*(1+alpha2)^(ht-1);
    for fd=1:ffd
        fi=ff(fd);
        xe(ht,fd)=Ri*cos(fi);
        ye(ht,fd)=Ri*sin(fi);
        ze(ht,fd)=0;
    end
end
hold on;surf(xe,ye,ze,plotT2)

shading interp
hidden off;
colorbar
%('YTickLabel',...
%{'10^8','10^9','10^10','10^11','10^12','10^13','10^14','10^15'})
title('Number Density distribution [#/m^3]')
%xlabel('X [R_C_e_r_e_s]',...
%    'FontSize',16)
%ylabel('Y [R_C_e_r_e_s]',...
%    'FontSize',16)
shading interp
hidden off;
colormap(jet)
colorbar
set(gca,'XTick',linspace(-3e6/Radii,3e6/Radii,7),...
    'XTickLabel',{'-300km','-200km','-100km','0km','100km','200km','300km'})
set(gca,'YTick',linspace(-3e6/Radii,3e6/Radii,7),...
    'YTickLabel',{'-300km','-200km','-100km','0km','100km','200km','300km'},...
    'FontSize',16)
Y_Label=8:13;
hei = 3e6/Radii;
axis([-hei hei -hei hei -hei hei]);
caxis([7 12])
colorbar('YTick',[7:12],'YTickLabel',...
    {'10^7','10^8','10^9','10^10','10^11',...
     '10^12'},...
     'FontSize',16)
title('Number density of exosphere of Ceres [molecules-m^-^2-s^-^1]')
 

%% Density=1e23 Surface
figure
clear x y z
s=zeros(3,1);
x=zeros(ftd,ffd);
y=x;
z=x;

for td=1:ftd
for fd=1:ffd
    for ht=hd:-1:1
        if(log10(copy_atm(ht,td,fd))>=12)
            Ri=ht*unih+Radii;
            thi=(td-0.5)/ftd*pi-pi/2;
            fi=fd-0.5/ffd*2*pi;
            x(td,fd)=Ri*cos(thi)*sin(fi);
            y(td,fd)=Ri*cos(thi)*cos(fi);
            z(td,fd)=Ri*sin(thi);
            plot3(x,y,z,'r')
            hold on;
            break;
        end
    end
end
end
surf(x,y,z)
%% Ice polar Cap

%figure(10)
hold on
x=zeros(ftd,ffd);
y=zeros(ftd,ffd);
z=zeros(ftd,ffd);
td=0;

for thi=linspace(-pi/2,pi/2,ftd)
    td=td+1;
    fd=0;
    for fi=linspace(0,2*pi,ffd)
        fd=fd+1;
            x(td,fd)=cos(thi)*(-sin(fi));
            y(td,fd)=cos(thi)*(-cos(fi));
            z(td,fd)=sin(thi);
    end
end
figure;
surf(x,y,z,log10(copy_horizon))
shading interp
% shading faceted
%% 
figure;imagesc(horizon)



copy_horizon=horizon;
copy=copy_horizon(:,1:ffd/2);
copy_horizon(:,ffd+1:ffd/2*3)=copy;
copy_horizon(:,1:ffd/2)=[];
imagesc(copy_horizon)
whos horizon
colorbar

%% 
for fd=1:ffd
figure(1)
    plotT12(:,:)=copy_atm(:,:,ceil(fd/2));
	plotT11(:,:)=copy_atm(:,:,ceil(fd/2+ffd/2));
    a1=plotT12(:,1:ftd/2);
    a2=plotT12(:,ftd/2+1:ftd);
    b1=cp_pT11(:,1:ftd/2);
    b2=cp_pT11(:,ftd/2+1:ftd);
    CC1=[b2 a1 a2 b1];
    surf(xe,ye,ze,CC1);axis tight;view(0,90);
    filename=['E' num2str(fd) '.png']
	print('-dpng',filename)
figure(2)
    surf(xt,yt,zt,plot00);
    view(fd,35)
    pause(0.1)
    filename=['T' num2str(fd) '.png']
	print('-dpng',filename)
    disp(fd)
end
cd('c:\test')












