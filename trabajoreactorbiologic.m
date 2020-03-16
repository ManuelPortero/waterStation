function reactorbiologic
 
clear
clc
close all
%Datos
Ca1=0; Cb1=0  ;Cc1=0 ;Ca2=0; Cb2=0  ;Cc2=0 ;Ca3=0; Cb3=0  ;Cc3=0 ; %mg/l
Ca4=0; Cb4=0  ;Cc4=0 ;Ca5=0; Cb5=0  ;Cc5=0 ;%mg/l 
mumax1=0.50; mumax2= 0.50; %dia-1 ; m
ks1=1; ks2=1 ; %mgNH4/l
ko1=0.74; ko2= 1.75; %mgO2 /l; 
Yxs1=0.18; Yxs2= 0.08 ;%mgDQO/mgN
kd1=0.11; kd2=0.11 ; %d-1 
ki1= 480; ki2=34; % mg/l 
cO = 8;% mg/l
X1= 1030; X12= 1030; X13= 1030; X14= 1030; X15= 1030;% mg /l 
X2= 480; X22= 480 ;  X23= 480 ;X24= 480 ;  X25= 480 ;% mg/l
F=300 ;%l/d
Fr=1000;%l/d
V=500;%l
NH4e=3000; %mg/l
NO2e=0;%mg/l
NO3e=0;%mg/l
opciones = odeset('Abstol',1e-5);
cond0=[X1 X2 Ca1 Cb1 Cc1 X12 X22 Ca2 Cb2 Cc2 X13 X23 Ca3 Cb3 Cc3 X14 X24 Ca4 Cb4 Cc4 X15 X25 Ca5 Cb5 Cc5];
 
[tt yy]=ode45(@bio,[0 120], cond0,opciones);
% una vez resuelta la diferencial renombramos las variables para graficar
X1a=yy(:,1);
X2a=yy(:,2);
CAa=yy(:,23);
CBa=yy(:,24);
CCa=yy(:,25);
%graficamos
figure
hold on
grid on
plot(tt,X1a,'b')
plot(tt,X2a,'r')
xlabel (' dias ')
ylabel (' mg biomasa /l ') 
legend('X1','X2')

hold off
figure
grid on
hold on
plot(tt,CAa,'g')
plot(tt,CCa,'m')
xlabel (' dias ')
ylabel (' mg/l ') 
legend('NH4','NO3','Location','NorthEastOutside')
hold off
figure
grid on
hold on
plot(tt,CBa,'b')
plot(tt,CCa,'m')
xlabel (' dias ')
ylabel (' mg/l ') 
legend('NO2','NO3','Location','NorthEastOutside')
hold off

 
function dy=bio(t,y)
 
%Resolvemos los balances de la diferencial
X1=y(1);X12=y(6);X13=y(11);X14=y(16);X15=y(21);
X2=y(2);X22=y(7);X23=y(12);X24=y(17);X25=y(22);
NH4=y(3);NH42=y(8);NH43=y(13);NH44=y(18);NH45=y(23);
NO2=y(4);NO22=y(9);NO23=y(14);NO24=y(19);NO25=y(24);  
NO3=y(5);NO32=y(10);NO33=y(15);NO34=y(20);NO35=y(25); 
  
muns= mumax1*(NH4/(ks1+NH4+(NH4^2/ki1)))*(cO/(ko1+cO));
munb= mumax2*(NO2/(ks2+NO2+(NO2^2/ki2)))*(cO/(ko2+cO));
 
rs=(muns*X1/Yxs1);
rn=(munb*X2/Yxs2);
 
dy(1,1)=(muns-kd1)*X1;
dy(2,1)=(munb-kd2)*X2;
dy(3,1)=(F/V)*NH4e+Fr/V*NH45-(F/V)*NH4-rs;
dy(4,1)= (F/V)*NO2e+Fr/V*NO25-(F/V)*NO2+rs-rn;
dy(5,1)= (F/V)*NO3e+Fr/V*NO35-(F/V)*NO3+rn;
 
muns2= mumax1*(NH42/(ks1+NH42+(NH42^2/ki1)))*(cO/(ko1+cO));
munb2= mumax2*(NO22/(ks2+NO22+(NO22^2/ki2)))*(cO/(ko2+cO));
 
rs2=(muns2*X12/Yxs1);
rn2=(munb2*X22/Yxs2);
 
dy(6,1)=(muns2-kd1)*X12;
dy(7,1)=(munb2-kd2)*X22;
dy(8,1)=(F/V)*NH4-(F/V)*NH42-rs2;
dy(9,1)= (F/V)*NO2-(F/V)*NO22+rs2-rn2;
dy(10,1)= (F/V)*NO3-(F/V)*NO32+rn2;
 
muns3= mumax1*(NH43/(ks1+NH43+(NH43^2/ki1)))*(cO/(ko1+cO));
munb3= mumax2*(NO23/(ks2+NO23+(NO23^2/ki2)))*(cO/(ko2+cO));
 
rs3=(muns3*X13/Yxs1);
rn3=(munb3*X23/Yxs2);
 
dy(11,1)=(muns3-kd1)*X13;
dy(12,1)=(munb3-kd2)*X23;
dy(13,1)=(F/V)*NH42-(F/V)*NH43-rs3;
dy(14,1)= (F/V)*NO22-(F/V)*NO23+rs3-rn3;
dy(15,1)= (F/V)*NO32-(F/V)*NO33+rn3;
 
 
muns4= mumax1*(NH44/(ks1+NH44+(NH44^2/ki1)))*(cO/(ko1+cO));
munb4= mumax2*(NO24/(ks2+NO24+(NO24^2/ki2)))*(cO/(ko2+cO));
 
rs4=(muns4*X14/Yxs1);
rn4=(munb4*X24/Yxs2);
 
dy(16,1)=(muns4-kd1)*X14;
dy(17,1)=(munb4-kd2)*X24;
dy(18,1)=(F/V)*NH43-(F/V)*NH44-rs4;
dy(19,1)= (F/V)*NO23-(F/V)*NO24+rs4-rn4;
dy(20,1)= (F/V)*NO33-(F/V)*NO34+rn4;
 
 
muns5= mumax1*(NH45/(ks1+NH45+(NH45^2/ki1)))*(cO/(ko1+cO));
munb5= mumax2*(NO25/(ks2+NO25+(NO25^2/ki2)))*(cO/(ko2+cO));
 
rs5=(muns5*X15/Yxs1);
rn5=(munb5*X25/Yxs2);
 
dy(21,1)=(muns5-kd1)*X15;
dy(22,1)=(munb5-kd2)*X25;
dy(23,1)=(F/V)*NH44-(F/V)*NH45-Fr/V*NH45-rs5;
dy(24,1)= (F/V)*NO24-(F/V)*NO25-Fr/V*NO25+rs5-rn5;
dy(25,1)= (F/V)*NO34-(F/V)*NO35-Fr/V*NO35+rn5;
 
 
end
end

