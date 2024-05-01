clear all;
close all;
clc;
% Исходные данные
H=2000;
L=1;
l1=0.5*L;
bv=0.1*L;
bn=0.05*L;
V=850;
% Параметры стандартной атмосферы
P=79501.4; % статическое давление
T=275.154; % температура в Кельвинах
R=287.053; % газовая постоянная воздуха
ro=1.00655; % плотность
a_inf=332.532; % скорость звука в невозмущённом набегающем потоке
k=1.4; % коэффициент адиабаты воздуха
AO=l1;
OC=L-l1;
BO=bv;
DO=bn;
% Определение абсолютных значений углов между гранями профиля и его хордой
beta1=atan2(BO,AO);
beta2=atan2(DO,AO);
beta3=atan2(BO,OC);
beta4=atan2(DO,OC);
M1=V/a_inf; % число Маха перед скачком уплотнения
P0=P*(1+(k-1)*M1^2/2)^(k/(k-1)); % давление торможения набегающего потока
M2=1; % число Маха за скачком уплотнения
beta_min=asin(1/M1); % минимальный угол наклона скачка уплотнения
% Определение угла наклона скачка уплотнения, при котором угол поворота потока максимален
A=2*k*M1^4;
B=M1^2*(4-(k+1)*M1^2);
C=-(M1^2*(k+1)+2);
beta_zv=asin(sqrt((-B+sqrt(B^2-4*A*C))/2/A));
% Метод половинного деления для определения угла наклона скачка уплотнения, при котором угол поворота потока максимален для заданного числа Маха
a=beta_min;
b=beta_zv;
err=10^(-5);
while abs(b-a)>=err
beta=(a+b)/2;
vsb1=2*(M1^2*(sin(beta))^2-1);
vsb2=tan(beta)*(M1^2*(k+cos(2*beta))+2);
Z=M2^2*(sin(beta-atan2(vsb1,vsb2)))^2;
Km=(1+(k-1)*M1^2*(sin(beta))^2/2)/(k*M1^2*(sin(beta))^2-(k-1)/2);
if Km<Z
b=beta;
elseif Km>Z
a=beta;
end
end
% Максимальный угол поворота потока
Teta=atan2(2*(M1^2*(sin(beta))^2-1),tan(beta)*(M1^2*(k+cos(2*beta))+2));
% Определение предельных углов атаки для заданного числа Маха
lowest=-(Teta-beta1);
highest=Teta-beta2;
% Выбор нескольких углов атаки для расчёта
step=(highest-lowest)/9;
alpha=[lowest:step:highest];
nomer=0;
fname='Results.txt';
fid=fopen(fname,'w');
z=0;
for i=1:length(alpha)
if (alpha(i)>0) || (alpha(i)<=0 && nomer==0)
z=z+1;
if alpha(i)<=0
alpha(i)=0;
nomer=1;
end
alpha(i)=alpha(i)*180/pi;
alpha(i)=fix(alpha(i));
alpha(i)=alpha(i)*pi/180;
% Реализация первой схемы обтекания с присоединённым скачком уплотнения
if beta1-alpha(i)>0
% Расчёт косых скачков уплотнения AB, AD
ab_Teta=beta1-alpha(i);
ad_Teta=alpha(i)+beta2;
for j=1:2
if j==1
T=ab_Teta;
elseif j==2
T=ad_Teta;
end
a=beta_min;
b=beta_zv;
err=10^(-5);
while abs(b-a)>=err
beta=(a+b)/2;
vsb1=2*(M1^2*(sin(beta))^2-1);
vsb2=tan(beta)*(M1^2*(k+cos(2*beta))+2);
K=atan2(vsb1,vsb2);
if K>T
b=beta;
elseif K<T
a=beta;
end
end
if j==1
beta_ab=beta;
elseif j==2
beta_ad=beta;
end
end
Pab=(2*k*M1^2*(sin(beta_ab))^2/(k+1)-(k-1)/(k+1))*P;
vsb1=(1+(k-1)*M1^2*(sin(beta_ab))^2/2)/(k*M1^2*(sin(beta_ab))^2-(k-1)/2);
Mab=sqrt(vsb1/(sin(beta_ab-ab_Teta))^2);
P0ab=Pab*(1+(k-1)*Mab^2/2)^(k/(k-1));
Pad=(2*k*M1^2*(sin(beta_ad))^2/(k+1)-(k-1)/(k+1))*P;
vsb1=(1+(k-1)*M1^2*(sin(beta_ad))^2/2)/(k*M1^2*(sin(beta_ad))^2-(k-1)/2);
Mad=sqrt(vsb1/(sin(beta_ad-ad_Teta))^2);
P0ad=Pad*(1+(k-1)*Mad^2/2)^(k/(k-1));
else
% Реализация второй схемы обтекания с присоединённым скачком уплотнения (в точке A вверх отходит волна разрежения, вниз – косой скачок уплотнения)
% Расчёт косого скачка уплотнения AD
ad_Teta=alpha(i)+beta2;
a=beta_min;
b=beta_zv;
err=10^(-5);
while abs(b-a)>=err
beta=(a+b)/2;
vsb1=2*(M1^2*(sin(beta))^2-1);
vsb2=tan(beta)*(M1^2*(k+cos(2*beta))+2);
K=atan2(vsb1,vsb2);
if K>ad_Teta
b=beta;
elseif K<ad_Teta
a=beta;
end
end
beta_ad=beta;
Pad=(2*k*M1^2*(sin(beta_ad))^2/(k+1)-(k-1)/(k+1))*P;
vsb1=(1+(k-1)*M1^2*(sin(beta_ad))^2/2)/(k*M1^2*(sin(beta_ad))^2-(k-1)/2);
Mad=sqrt(vsb1/(sin(beta_ad-ad_Teta))^2);
P0ad=Pad*(1+(k-1)*Mad^2/2)^(k/(k-1));
% Расчёт течения разрежения AB
beta_ab=alpha(i)-beta1;
omega=sqrt((k+1)/(k-1))*atan(sqrt((k-1)/(k+1))*sqrt(M1^2-1))-atan(sqrt(M1^2-1));
omegaAB=omega+beta_ab;
a=1; b=10;
err=10^(-5);
while abs(b-a)>=err
M=(a+b)/2;
K=sqrt((k+1)/(k-1))*atan(sqrt((k-1)/(k+1))*sqrt(M^2-1))-atan(sqrt(M^2-1));
if K>omegaAB
b=M;
elseif K<omegaAB
a=M;
end
end
omegaAB=omegaAB*180/pi;
Mab=M;
P0ab=P0;
Pab=P0ab*(1+((k-1)*Mab^2)/2)^(-k/(k-1));
end
% Расчёт течений разрежения BC, DC
beta_bcdc=beta1+beta3;
for j=1:2
if j==1
omega=sqrt((k+1)/(k-1))*atan(sqrt((k-1)/(k+1))*sqrt(Mab^2-1))-atan(sqrt(Mab^2-1));
omega2=omega+beta_bcdc;
omegaBC=omega2*180/pi;
elseif j==2
omega=sqrt((k+1)/(k-1))*atan(sqrt((k-1)/(k+1))*sqrt(Mad^2-1))-atan(sqrt(Mad^2-1));
omega2=omega+beta_bcdc;
omegaDC=omega2*180/pi;
end
a=1; b=10;
err=10^(-5);
while abs(b-a)>=err
M=(a+b)/2;
K=sqrt((k+1)/(k-1))*atan(sqrt((k-1)/(k+1))*sqrt(M^2-1))-atan(sqrt(M^2-1));
if K>omega2
b=M;
elseif K<omega2
a=M;
end
end
if j==1
Mbc=M;
elseif j==2
Mdc=M;
end
end
P0bc=P0ab;
Pbc=P0bc*(1+((k-1)*Mbc^2)/2)^(-k/(k-1));
P0dc=P0ad;
Pdc=P0dc*(1+((k-1)*Mdc^2)/2)^(-k/(k-1));
% Расчёт аэродинамических сил и коэффициентов
X1=Pab*BO+Pad*DO-Pbc*BO-Pdc*DO;
Y1=Pad*AO+Pdc*OC-Pab*AO-Pbc*OC;
X=X1*cos(alpha(i))+Y1*sin(alpha(i));
Y=-X1*sin(alpha(i))+Y1*cos(alpha(i));
ctau=X1/(0.5*P*M1^2*L*k);
cn=Y1/(0.5*P*M1^2*L*k);
cx=ctau*cos(alpha(i))+cn*sin(alpha(i));
cy=-ctau*sin(alpha(i))+cn*cos(alpha(i));
% Приближенный расчёт
if alpha(i)==0
phi=sqrt(cx*sqrt(M1^2-1))/2;
end
cxappr=4*(alpha(i)^2+phi^2)/sqrt(M1^2-1);
cyappr=4*alpha(i)/sqrt(M1^2-1);
Cx(z)=cx;
Cy(z)=cy;
Cxappr(z)=cxappr;
Cyappr(z)=cyappr;
% Запись результатов в файл
alfa=alpha(i)*180/pi;
beta_ab=beta_ab*180/pi;
beta_ad=beta_ad*180/pi;
fprintf(fid,'%s\r\n',' ');
fprintf(fid,'%s',' Alpha=');
fprintf(fid,'%.0f',alfa);
if beta1-alpha(i)>0
fprintf(fid,'%s',' beta_ab=');
fprintf(fid,'%.3f',beta_ab);
else
fprintf(fid,'%s',' omega_ab=');
fprintf(fid,'%.3f',omegaAB);
end
fprintf(fid,'%s',' beta_ad=');
fprintf(fid,'%.3f',beta_ad);
fprintf(fid,'%s',' Pab=');
fprintf(fid,'%.3e',Pab);
fprintf(fid,'%s',' P0ab=');
fprintf(fid,'%.3e',P0ab);
fprintf(fid,'%s',' Mab=');
fprintf(fid,'%.4f',Mab);
fprintf(fid,'%s',' Pbc=');
fprintf(fid,'%.3e',Pbc);
fprintf(fid,'%s',' P0bc=');
fprintf(fid,'%.3e',P0bc);
fprintf(fid,'%s',' Mbc=');
fprintf(fid,'%.4f',Mbc);
fprintf(fid,'%s',' Pad=');
fprintf(fid,'%.3e',Pad);
fprintf(fid,'%s',' P0ad=');
fprintf(fid,'%.3e',P0ad);
fprintf(fid,'%s',' Mad=');
fprintf(fid,'%.4f',Mad);
fprintf(fid,'%s',' Pdc=');
fprintf(fid,'%.3e',Pdc);
fprintf(fid,'%s',' P0dc=');
fprintf(fid,'%.3e',P0dc);
fprintf(fid,'%s',' Mdc=');
fprintf(fid,'%.4f',Mdc);
fprintf(fid,'%s',' omega_bc=');
fprintf(fid,'%.3f',omegaBC);
fprintf(fid,'%s',' omega_dc=');
fprintf(fid,'%.3f',omegaDC);
fprintf(fid,'%s',' X=');
fprintf(fid,'%.3e',X);
fprintf(fid,'%s',' Y=');
fprintf(fid,'%.3e',Y);
fprintf(fid,'%s',' Cx=');
fprintf(fid,'%.3f',cx);
fprintf(fid,'%s',' Cy=');
fprintf(fid,'%.3f',cy);
fprintf(fid,'%s',' Cx_приближенное=');
fprintf(fid,'%.3f',cxappr);
fprintf(fid,'%s',' Cy_приближенное=');
fprintf(fid,'%.3f',cyappr);
% Построение картины обтекания
beta_ab=beta_ab*pi/180;
beta_ad=beta_ad*pi/180;
figure;
hold on;
xlim=([-0.2,1.2]);
ylim=([-0.2,0.2]);
axis equal;
% Построение геометрии
% Проекции на ось Ox узловых точек профиля
x=[0,cos(beta1-alpha(i))*sqrt(AO^2+BO^2),(AO+OC)*cos(alpha(i)),cos(beta2+alpha(i))*sqrt(AO^2+DO^2),0];
% Проекции на ось Oy узловых точек профиля
y=[0,sin(beta1-alpha(i))*sqrt(AO^2+BO^2),-(AO+OC)*sin(alpha(i)),-sin(beta2+alpha(i))*sqrt(AO^2+DO^2),0];
plot(x,y,'LineWidth',3,'Color','black');
str=strcat('Угол атаки= ',num2str(alpha(i)*180/pi),' градусов');
title(str);
s=0.3;
% Косой скачок уплотнения A и C
if beta1-alpha(i)>0
plot([x(1),x(1)+s*cos(beta_ab)],[y(1),y(1)+s*sin(beta_ab)],'LineWidth',2,'Color','red');
plot([x(3),x(3)+s*cos(-beta_bcdc)],[y(3),y(3)+s*sin(-beta_bcdc)],'LineWidth',2,'Color','red');
else
% Волна разрежения AB (веер характеристик)
amab=asin(1/Mab);
b=beta1+alpha(i);
b2=beta1-alpha(i);
am=[amab-b2:5*pi/180:amab+b];
for nn=1:length(am)
plot([x(1),x(1)+s*cos(am(nn))],[y(1),y(1)+s*sin(am(nn))],'--');
end
plot([x(1),x(1)+s*cos(amab+b)],[y(1),y(1)+s*sin(amab+b)],'--','LineWidth',2);
plot([x(1),x(1)+s*cos(amab-b2)],[y(1),y(1)+s*sin(amab-b2)],'--','LineWidth',2);
% Волна разрежения вблизи задней кромки, в нижней части профиля (веер характеристик)
amdc=asin(1/Mdc);
b=beta4+alpha(i);
b2=beta4-alpha(i);
am=[amdc-b2:5*pi/180:amdc+b];
for nn=1:length(am)
plot([x(3),x(3)+s*cos(am(nn))],[y(3),y(3)-s*sin(am(nn))],'--');
end
plot([x(3),x(3)+s*cos(amdc+b)],[y(3),y(3)-s*sin(amdc+b)],'--','LineWidth',2);
plot([x(3),x(3)+s*cos(amdc-b2)],[y(3),y(3)-s*sin(amdc-b2)],'--','LineWidth',2);
end
% Косой скачок уплотнения AD и скачок вблизи задней кромки, в нижней части профиля
plot([x(3),x(3)+s*cos(beta_bcdc)],[y(3),y(3)+s*sin(beta_bcdc)],'LineWidth',2,'Color','red');
plot([x(1),x(1)+s*cos(-beta_ad)],[y(1),y(1)+s*sin(-beta_ad)],'LineWidth',2,'Color','red');
% Волны разрежения AD, BC, DC (пучки характеристик)
amab=asin(1/Mab);
ambc=asin(1/Mbc);
amad=asin(1/Mad);
amdc=asin(1/Mdc);
b=beta3+alpha(i);
b2=beta3-alpha(i);
am=[ambc-b2:5*pi/180:amab+b];
for nn=1:length(am)
plot([x(2),x(2)+s*cos(am(nn))],[y(2),y(2)+s*sin(am(nn))],'--');
end
plot([x(2),x(2)+s*cos(amab+b)],[y(2),y(2)+s*sin(amab+b)],'--','LineWidth',2);
plot([x(2),x(2)+s*cos(ambc-b2)],[y(2),y(2)+s*sin(ambc-b2)],'--','LineWidth',2);
am=[amdc-b2:10*pi/180:amad+b];
for nn=1:length(am)
plot([x(4),x(4)+s*cos(am(nn))],[y(4),y(4)-s*sin(am(nn))],'--');
end
plot([x(4),x(4)+s*cos(amad+b)],[y(4),y(4)-s*sin(amad+b)],'--','LineWidth',2);
plot([x(4),x(4)+s*cos(amdc-b2)],[y(4),y(4)-s*sin(amdc-b2)],'--','LineWidth',2);
str2=strcat('картина обтекания ',num2str(z));
saveas(gcf,str2,'jpg');
end
end
% Построение поляры
figure;
hold on;
grid on;
plot(Cx,Cy,'-*','LineWidth',2);
plot(Cxappr,Cyappr,'r -*','LineWidth',2);
legend('Точное','Приближенное',0);
xlabel('Cx');
ylabel('Cy');
saveas(gcf,'Поляра','jpg');
fclose(fid);