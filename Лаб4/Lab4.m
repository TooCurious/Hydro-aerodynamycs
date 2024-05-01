clc
clear all
%Необходимо определить tkr
y=1.4;
pn=1;
p0kr=pn*((y+1)/2)^(y/(y-1));
t=[0;0.2;0.4;0.6;0.8;1.0;1.2;1.4;1.6];
p=[6;4.7;3.5;2.5;2.0;1.6;1.4;1.2;1];
p00=6;
P=p/p00;
h=0:0.01:1.6;
gr1=pchip(t,P,h);
figure(1)
plot(t,P,'x',h,gr1)
grid on
xlabel('t,c');
ylabel('p0/p00');
hold on
pkr=p0kr/p00;
pkrline=line([0 1.6], [pkr pkr],'LineStyle','--');
tkr=0.8476;
tkrline=line([tkr tkr], [0 1],'LineStyle','--');
hold off

B=0.13;
T00=293;
V=0.009;
d=0.008;
R=287;
Fa=pi*(d/2)^2;
M=1;

%Сверхкритический режим
tsverh=[0;0.3;0.4;0.6;0.8];
for i=1:1:5
    pteor(i)=(1+B*tsverh(i)).^(-2*y/(y-1));%теоретическое изменение давлений p0/p00
    ro(i)=pteor(i).^(1/y);%теоретическое изменение плотности ro0/ro00
    T(i)=pteor(i)^((y-1)/y);%теоретическое изменение температуры T0/T00
    a(i)=pteor(i)^((y-1)/2*y);%теоретическое изменение скорости звука a0/a00 
end

for i=1:1:5
    p0(i)=pteor(i)*p00*10^5;
    ro00=(p00*10^5)/(R*T00);
    ro0(i)=ro00*ro(i);
    T0(i)=T(i)*T00;
    a00=sqrt(y*R*T00);
    a0(i)=a(i)*a00;
end

for i=1:1:5
    pa(i)=p0(i)*(1+((y-1)/2)*M^2)^(-y/(y-1));
    roa(i)=ro0(i)*(1+((y-1)/2)*M^2)^(-1/(y-1));
    Ta(i)=T0(i)*(1+((y-1)/2)*M^2)^(-1);
    aa(i)=a0(i)*(1+((y-1)/2)*M^2)^(-1/2);
    Va(i)=M*aa(i);
    Q(i)=pa(i)*roa(i)*Fa;
end
tsverh=[0;0.2;0.4;0.6;0.8];
p0sverh=[600000;443860;332470;251900;192890];
ro0sverh=[7.1351;5.7531;4.6801;3.8386;3.1722];
T0sverh=[293;268.8;247.5;228.7;211.9];
a0sverh=[343;318;310;296;282];
Vasverh=[313;290;283;270;258];
Qsverh=[0.0720668;0.0459862;0.0391932;0.0292771;0.0213002]
T1=table(tsverh,p0sverh,ro0sverh,T0sverh,a0sverh,Vasverh,Qsverh)

%Докритический режим
tdokr=[1.2;1.4;1.6];
C=y*sqrt((2*y*R*T00)/(y-1))*(Fa/V)*(pn/p00)^((y-1)/2*y);
for i=1:1:3
    J=(tdokr(i)-tkr)*C;
end
%по таблице определяем значения
p0d(1)=1.35*pn*10^5;
p0d(2)=1.15*pn*10^5;
p0d(3)=1.05*pn*10^5;
p0d=[135000;115000;105000];
for i=1:1:3
    rodokr(i)=(p0d(i)/(p00*10^5))^(1/y)*ro00;
    T0dokr(i)=(p0d(i)/(p00*10^5))^((y-1)/y)*T00;
    a0dokr(i)=(p0d(i)/(p00*10^5))^((y-1)/(2*y))*a00;
end
%c помощью Маthcad найдём значения для М
M1=0.669;
M2=0.451;
M3=0.265;
%Найдём значения для параметров в отверстии
roadokr=p00*(pn/p00)^(1/y);
Tadokr=(pn*10^5)/(R*roadokr);
aadokr=sqrt(y*R*Tadokr);
Vadokr(1)=M1*aadokr;
Vadokr(2)=M2*aadokr;
Vadokr(3)=M3*aadokr;
Qdokr(1)=roadokr*Vadokr(1)*Fa;
Qdokr(2)=roadokr*Vadokr(2)*Fa;
Qdokr(3)=roadokr*Vadokr(3)*Fa;
Qdokr=[0.0163;0.0110;0.0064];
Vadokr=[194;131;77];
Mdokr=[0.669;0.451;0.265];
rodokr=[2.4585;2.1925;2.0545];
T0dokr=[191;183;178];
a0dokr=[277;270;267];
T2=table(tdokr,p0d,rodokr,T0dokr,a0dokr,Vadokr,Qdokr)

%Построение графиков
t=[0;0.2;0.4;0.6;0.8;1.0;1.2;1.4;1.6];
p=[6;4.7;3.5;2.5;2.0;1.6;1.4;1.2;1];
tteor=[0;0.2;0.4;0.6;0.8;1.2;1.4;1.6];
ptteor=[600000;443860;332470;251900;192890;135000;115000;105000];
p00=6;
P=p/p00;
Pteor=ptteor/(p00*10^5);
h=0:0.01:1.6;
gr1=pchip(t,P,h);
gr2=pchip(tteor,Pteor,h);
figure(2)
plot(t,P,'x',h,gr1)
grid on
xlabel('t,c');
ylabel('p0/p00');
hold on
plot(tteor,Pteor,'o',h,gr2,'--')
legend('','Экспериментальное','','Теоретическое')
hold off 

ro=[7.1351;5.7531;4.6801;3.8386;3.1722;2.4585;2.1925;2.0545];
gr3=pchip(tteor,ro,h);
figure(3)
plot(tteor,ro,'x',h,gr3);
grid on
xlabel('t,c');
ylabel('ro, кг/м^3');

T=[293;268.8;247.5;228.7;211.9;191;183;178];gr3=pchip(tteor,ro,h);
gr4=pchip(tteor,T,h);
figure(4)
plot(tteor,T,'x',h,gr4);
grid on
xlabel('t,c');
ylabel('T, K');

a=[343;318;310;296;282;277;270;267];
gr5=pchip(tteor,a,h);
figure(5)
plot(tteor,a,'x',h,gr5);
grid on
xlabel('t,c');
ylabel('a, m/c');

Va=[313;290;283;270;258;194;131;77];
gr6=pchip(tteor,Va,h);
figure(6)
plot(tteor,Va,'x',h,gr6);
grid on
xlabel('t,c');
ylabel('V, m/c');

Q=[0.0720668;0.0459862;0.0391932;0.0292771;0.0213002;0.0163;0.0110;0.0064];
gr7=pchip(tteor,Q,h);
figure(7)
plot(tteor,Q,'x',h,gr7);
grid on
xlabel('t,c');
ylabel('Q, kg/c');

M=[1;1;1;1;1;0.669;0.451;0.265];
gr8=pchip(tteor,M,h);
figure(8)
plot(tteor,M,'x',h,gr8);
grid on
xlabel('t,c');
ylabel('M');