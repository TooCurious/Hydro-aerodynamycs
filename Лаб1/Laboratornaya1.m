clear all
clc
alpha=[0;3;6;9;12];
Xi=[142;165;200;230;386];
Yi1=[0;31;64;91;85];
Yi2=[0;48;176;331;539];
dH=[150;150;150;150;150];
g=9.81;
d=0.116;
Kx=2.47;
Ky=2.5;
AB=0.2;

%Вычисления состовляющих силы R, cкорости и скоростного напора%
X=(Kx*Xi-200)*g*0.001;
Y1=(Ky*Yi1)*g*0.001;
Y2=(Ky*Yi2)*g*0.001;
Y=Y1+Y2;
T=table(alpha,X,Y1,Y2,Y);
Hsr=mean(dH);
v=4.111*Hsr.^(1/2);
q=g*Hsr;
Sm=(pi*d^2)/4;
Cx=X/(q*Sm);
Cy=Y/(q*Sm);
Cy2=Y2/(q*Sm);
attack_angle=alpha*pi/180;
T2=table(alpha,attack_angle,Cx,Cy,Cy2);

%коэффициенты аппроксимации Сx,Cy,Cy1%
Cx0 =     0.09083;
a1 =       8.069;
a2 =       3.171;
a3 =       36.08;
a4 =       2.062;
a5 =        46.4;

%Нахождение положения центра давления%

AC=AB*Cy2.*cos(attack_angle)./(Cy.*cos(attack_angle)+Cx.*sin(attack_angle));
AC(1)=AB*(a4/(a2+Cx0));
T3=table(attack_angle,AC)

%определение центра тяжести модели%
R=0.058;
r=0.0435;
H1=0.216;
H2=0.116;
V1=(2/3)*pi*R^3;
Z1=(3/8)*R;
V2=pi*(R^2)*H1;
Z2=H1/2
V3=(pi*H2/3)*(R^2+r*R+r^2);
Z3=(H2/4)*(3-((r*R+2*R^2)/(R^2+r*R+r^2)))
X1=-Z1-0.008;
X2=Z2-0.008;
X3=Z3+H1-0.008;
x0=(V1*X1+V2*X2+V3*X3)/(V1+V2+V3)
Xd=x0-AC
T4=table(alpha,attack_angle,Xd)

%Момент тангажа%
Mz=(Cy.*cos(attack_angle)+Cx.*sin(attack_angle)).*Xd.*q*Sm
Mz0=(Cy.*cos(attack_angle)+Cx.*sin(attack_angle)).*Xd.*q*Sm

%Коэффициент момента%
l=R+H1+H2
Cmz=(Cy.*cos(attack_angle)+Cx.*sin(attack_angle)).*Xd/l
Cmzpr=-(a2+Cx0)*Xd(1)/l*attack_angle
T5=table(alpha,attack_angle,Cx,Cy,Cy2,Xd,AC,Mz,Cmz)
