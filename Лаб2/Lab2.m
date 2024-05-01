clc
clear all 
%вычисление коэффициента давления%
hi=[147;80;32;-160;-186;-113;-47;-10;64;50;-20;5;-2;8;15];
hp=148;
hi1=(hi(1)+hi(9))/2;
hi3=(hi(1)+hi(2))/2;
hifull=[hi1;147;80;hi3;32;-160;-186;-113;-47;-10;64;50;-20;5;-2;8;15];
x=[0;0.2;1.1;1.5;2.6;13;25;47;91;120;2;3;4;5;44;89;120];
y=[0;2;3.5;5;6.6;13.2;16.3;18;13;7.5;-3;-4;-4.6;-5;-6.8;-4;-2];
Cpi=hifull/hp;
b=150;
x1=x/b;
y1=y/b;
Numberofpoint=[1;2;3;4;5;5.1;5.2;6;7;8;9;10;11;12;13;14;15];
T1=table(Numberofpoint,Cpi,x1,y1);

%Построение диаграммы Cpi(x)
Cpiup=[Cpi(1);Cpi(2);Cpi(3);Cpi(4);Cpi(5);Cpi(6);Cpi(7);Cpi(8);Cpi(9);Cpi(10);0];
Cpidown=[Cpi(1);Cpi(11);Cpi(12);Cpi(13);Cpi(14);Cpi(15);Cpi(16);Cpi(17);0];
Xup=[x1(1);x1(2);x1(3);x1(4);x1(5);x1(6);x1(7);x1(8);x1(9);x1(10);1];
Xdown=[x1(1);x1(11);x1(12);x1(13);x1(14);x1(15);x1(16);x1(17);1];
h=0:0.001:1;
lineup=pchip(Xup,Cpiup,h);
figure(1)
gr11=plot(h,lineup);
grid on
hold on
linedown=pchip(Xdown,Cpidown,h);
gr12=plot(h,linedown);
plot(x1(1),Cpi(1),'x', x1(2),Cpi(2),'+',x1(3),Cpi(3),'o',x1(4),Cpi(4),'s',x1(5),Cpi(5),'>')
plot(x1(6),Cpi(6),'<',x1(7),Cpi(7),'^',x1(8),Cpi(8),'*',x1(9),Cpi(9),'h',x1(10),Cpi(10),'p')
plot(x1(11),Cpi(11),'v',x1(12),Cpi(12),'h',x1(13),Cpi(13),'d',x1(14),Cpi(14),'s')
plot(x1(15),Cpi(15),'d',x1(16),Cpi(16),'x',x1(17),Cpi(17),'o')
legend('gr11','gr12','1','2','3','4','5','5a','5b','6','7','8','9','10','11','12','13','14','15')

%Построение сетки для нахождения значений для подстановки в формулу
%Симпсона для вычисления коэффициента cn(cnorm)
Xmax=1; Xmin=0; N=11;
hsimp=(Xmax-Xmin)/(N-1);
L1=line([Xmin Xmin],[-1.5;1.5],'LineStyle','--');
L2=line([Xmin+hsimp Xmin+hsimp],[-1.5;1.5],'LineStyle','--');
L3=line([Xmin+2*hsimp Xmin+2*hsimp],[-1.5;1.5],'LineStyle','--');
L4=line([Xmin+3*hsimp Xmin+3*hsimp],[-1.5;1.5],'LineStyle','--');
L5=line([Xmin+4*hsimp Xmin+4*hsimp],[-1.5;1.5],'LineStyle','--');
L6=line([Xmin+5*hsimp Xmin+5*hsimp],[-1.5;1.5],'LineStyle','--');
L7=line([Xmin+6*hsimp Xmin+6*hsimp],[-1.5;1.5],'LineStyle','--');
L8=line([Xmin+7*hsimp Xmin+7*hsimp],[-1.5;1.5],'LineStyle','--');
L9=line([Xmin+8*hsimp Xmin+8*hsimp],[-1.5;1.5],'LineStyle','--');
L10=line([Xmin+9*hsimp Xmin+9*hsimp],[-1.5;1.5],'LineStyle','--');
L11=line([Xmin+10*hsimp Xmin+10*hsimp],[-1.5;1.5],'LineStyle','--');
%находим по графику точки пересечения сетки с Срxup и Срxdown и получаем значения
delta_Cpxup=[0.7128;-1.131;-1.205;-0.7991;-0.6008;-0.4577;-0.3268;-0.1769;-0.0677;-0.0207;0];
delta_Cpxdown=[0.7128;0.0261;0.0001;-0.0134;0.0005;0.0284;0.0056;0.0845;0.1014;0.0719;0];
delta_Cnorm=delta_Cpxdown-delta_Cpxup;
%Находим Cnorm по формуле Симпсона
Cnorm1=4*(delta_Cnorm(2)+delta_Cnorm(4)+delta_Cnorm(6)+delta_Cnorm(8)+delta_Cnorm(10));
Cnorm2=2*(delta_Cnorm(3)+delta_Cnorm(5)+delta_Cnorm(7)+delta_Cnorm(9));
%delta_Cnorm(1)=delta_Cnorm(11)=0
Cnorm=(Cnorm1+Cnorm2)/30;

%Построение диаграммы Cpi(y)
Cpilob=[Cpi(15);Cpi(14);Cpi(13);Cpi(12);Cpi(11);Cpi(1);Cpi(2);Cpi(3);Cpi(4);Cpi(5);Cpi(6);Cpi(7);Cpi(8)];
Cpikorm=[Cpi(15);Cpi(16);Cpi(17);0;Cpi(10);Cpi(9);Cpi(8)];
Yup=[y1(15);y1(14);y1(13);y1(12);y1(11);y1(1);y1(2);y1(3);y1(4);y1(5);y1(6);y1(7);y1(8)];
Ydown=[y1(15);y1(16);y1(17);0;y1(10);y1(9);y1(8)];
Ymin=min(y1);
Ymax=max(y1);
hy=Ymin:0.001:Ymax;
linelob=pchip(Yup,Cpilob,hy);
figure(2)
gr21=plot(hy,linelob);
grid on
hold on
linekorm=spline(Ydown,Cpikorm,hy);
gr22=plot(hy,linekorm);
plot(y1(1),Cpi(1),'x', y1(2),Cpi(2),'+',y1(3),Cpi(3),'o',y1(4),Cpi(4),'s',y1(5),Cpi(5),'>')
plot(y1(6),Cpi(6),'<',y1(7),Cpi(7),'^',y1(8),Cpi(8),'*',y1(9),Cpi(9),'h',y1(10),Cpi(10),'p')
plot(y1(11),Cpi(11),'v',y1(12),Cpi(12),'h',y1(13),Cpi(13),'d',y1(14),Cpi(14),'s')
plot(y1(15),Cpi(15),'d',y1(16),Cpi(16),'x',y1(17),Cpi(17),'o')
legend('gr21','gr22','1','2','3','4','5','5a','5b','6','7','8','9','10','11','12','13','14','15')
hold on
%Построение сетки для нахождения значений для подстановки в формулу
%Симпсона для вычисления коэффициента ct(ctau)
N=11;
hysimp=(Ymax-Ymin)/(N-1);
L21=line([Ymin Ymin], [-2;2],'LineStyle','--');
L22=line([hysimp+Ymin hysimp+Ymin], [-2;2],'LineStyle','--');
L23=line([2*hysimp+Ymin 2*hysimp+Ymin], [-2;2],'LineStyle','--');
L24=line([3*hysimp+Ymin 3*hysimp+Ymin], [-2;2],'LineStyle','--');
L25=line([4*hysimp+Ymin 4*hysimp+Ymin], [-2;2],'LineStyle','--');
L26=line([5*hysimp+Ymin 5*hysimp+Ymin], [-2;2],'LineStyle','--');
L27=line([6*hysimp+Ymin 6*hysimp+Ymin], [-2;2],'LineStyle','--');
L28=line([7*hysimp+Ymin 7*hysimp+Ymin], [-2;2],'LineStyle','--');
L29=line([8*hysimp+Ymin 8*hysimp+Ymin], [-2;2],'LineStyle','--');
L210=line([9*hysimp+Ymin 9*hysimp+Ymin], [-2;2],'LineStyle','--');
L211=line([10*hysimp+Ymin 10*hysimp+Ymin], [-2;2],'LineStyle','--');
%находим по графику точки пересечения сетки с Срilob и Срikorm и получаем значения
delta_Cpilob=[-0.0135;0.0675;0.5347;0.8188;0.6099;0.6312;-0.1625;-0.7114;-1.0682;-1.224;-0.7635];
delta_Cpikorm=[-0.0135;0.0408;0.0979;-0.0318;-0.0733;-0.0630;-0.0765;-0.1679;-0.3201;-0.5245;-0.7635];
delta_Ct=delta_Cpilob-delta_Cpikorm
%Находим Ct по формуле Симпсона
Ctau1=4*(delta_Ct(2)+delta_Ct(4)+delta_Ct(6)+delta_Ct(8)+delta_Ct(10));
Ctau2=2*(delta_Ct(3)+delta_Ct(5)+delta_Ct(7)+delta_Ct(9));
%delta_Ct(1)=delta_Ct(11)=0
Ctau=((Ymax-Ymin)*(Ctau1+Ctau2))/30;

%Нахождение cx,cy
Cnorm=0.5249;
Ctau=0.0104;
cx=-Ctau*cosd(8)+Cnorm*sind(8)
cy=-Ctau*sind(8)+Cnorm*cosd(8)

%Для нахождения cmz по формуле Симпсона можно взять значения delta_Cpxup,
%delta_Cpxdown и добавить только:
delta_x=[0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1];
delta_Cmz=delta_Cpxdown.*delta_x-delta_Cpxup.*delta_x
%Тогда cmz по формуле Симпсона
Cmz1=4*(delta_Cmz(2)+delta_Cmz(4)+delta_Cmz(6)+delta_Cmz(8)+delta_Cmz(10));
Cmz2=2*(delta_Cmz(3)+delta_Cmz(5)+delta_Cmz(7)+delta_Cmz(9));
Cmz=-(Cmz1+Cmz2)/30;
%Нахождение центра давления
Xd=-Cmz/Cnorm