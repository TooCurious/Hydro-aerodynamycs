clear all;
clc;
alpha = [0, 3, 6, 9, 12];
Xc=[142, 165, 200, 230, 386];
Y1c=[0, 31, 64, 91, 85];
Y2c=[0, 48, 176, 331, 539];
dH=[150, 150, 150, 150, 150];
Kx=2.47;
Ky=2.5;
d=0.116;
g=9.81;
AB=0.2;
R=0.058;
H1=0.216;
H2=0.116;
%Расчеты
x=Kx*Xc-200;
y1=Ky*Y1c;
y2=Ky*Y2c;
y=y1+y2;
xn=x*g/1000;
yn=y*g/1000;
y2n=y2*g/1000;
u=4.111*(dH).^(1/2);
q=g*dH;
S=pi/4*d^2;
z=q*S;
Cx=xn./z;
Cy=yn./z;
Cy2=y2n./z;
attack_angle1 = alpha*pi/180;
AC=AB*Cy2.*cos(attack_angle1)./(Cy.*cos(attack_angle1)+Cx.*sin(attack_angle1));

%
Cx0 =     0.09083;
a1 =       8.069;
a2 =       3.168;
a3 =       36.05;
a4 =       2.062;
a5 =        46.4;

AC(1)=AB*a4/(a2+Cx0);
%
V1=2/3*pi*R^3;
xx1=R*3/8;
X1=-xx1-0.008;
V2=pi*R^2*H1;
xx2=H1/2;
X2=xx2-0.008;
r=0.0435;
V3=pi*H2/3*(R^2+R*r+r^2);
xx3=H2/4*(3-(R*r+2*R^2)/(R^2+R*r+r^2));
X3=xx3+0.216-0.008;
x0=(V1*X1+V2*X2+V3*X3)/(V1+V2+V3)

AO=x0;
OC0=AC(1)-AO;
OC=AC-AO;



