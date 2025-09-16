% Runge-Kutta en Dimensión. Se hacen m  pasos 
%  para aproximar la solución de
% Y'=F(t,Y), t en [t0,T] con Y(t0)=Y0.
% El intervalo [t0,T] se particiona en m+1 puntos 
% igualmente espaciados contenidos en el vector t.
% La función F se especifica en este archivo
% La solución Y al tiempo t(i) es aproximado
% por el renglón i-esimo de YY, esto es:
% YY(i,:) aproxima a Y(t(i)).
%
% Para graficar la primera componente
% del vector solución, escribir
%  plot(t,YY(:,1))
% La segunda componente del vector solucion se
% grafica con 
%  plot(t,YY(:,2))
% 
%  
% es conveniente escribirlo en la forma
% Y'=F(t,Y), donde Y es un vector con 3 componentes
%  Y=(Y(1),Y(2),Y(3)) 
% Estas 3 líneas se escriben dentro de la función F
% de abajo.  


%Simular numéricamente la dinámica de contagios usando el método de
%Runge-Kutta

function [t,YY] = R_K_n(t0, T, Y0)
m=10000;  % Número de pasos de Runge-Kutta
h = (T-t0)/m;
t=linspace(t0,T,m+1);
n=length(Y0);
YY=zeros(m+1,n);
YY(1,:)=Y0;
 
for i=1:m
   k1=F(t(i),YY(i,:));
   k2=F(t(i)+h/2,YY(i,:)+h*k1/2);
   k3=F(t(i)+h/2,YY(i,:)+h*k2/2);
   k4=F(t(i)+h,YY(i,:)+h*k3);
   YY(i+1,:)=YY(i,:)+h/6*(k1+2*k2+2*k3+k4);
 end
 
 function W=F(t,Y)
 n=length(Y);
 W=zeros(1,n);

 % Modelo COVID-19 CDMX. Definimos las constantes:
 R0=2.83;
 Dinf=2.9;Dinc=5.2;Pgrave=0.138;Dr=14;Dhosp=4;
 Picu=0.05; Drh=12; Pm=0.03; Dicu=1; Dricu=7; Dm=8;

%Definimos: 
W(1)=-(R0/Dinf)*Y(3)*Y(1);
W(2)=(R0/Dinf)*Y(3)*Y(1)-Y(2)/Dinc;
W(3)=Y(2)/Dinc-Y(3)/Dinf;
W(4)=(1-Pgrave)/Dinf*Y(3)-Y(4)/Dr;
W(5)=(Pgrave)/Dinf*Y(3)-Y(5)/Dhosp;
W(6)=Y(5)/Dhosp-(1-Picu)/Drh*Y(6)-Picu/Dicu*Y(6);
W(7)=Picu/Dicu*Y(6)-(1-Pm)/Dricu*Y(7)-Pm/Dm*Y(7);
W(8)=Y(4)/Dr+(1-Picu)/Drh*Y(6)+(1-Pm)/Dricu*Y(7);

function [y]=R0(t)
  c1=2.8; c2=0.5; c3=2.8;
  I1=[0,15];I2=[16,29];
  if t<=I2(1)
    y=c1;
  end 
  if t>=I2(1) && t<=I2(2)
    y=c2;
  end  
  if t>=I2(2)+1
    y=c3;
  end  
  if t>=I1(2) && t<=I2(1)
    y=c1+(c2-c1)*(t-I1(2))/(I2(1)-I1(2));
  end
  if t>=I2(2) && t<=I2(2)+1
    y=c3-(c2-c3)*(t-I2(2)-1);
  end
  
  %definir un vector de valores iniciales cuando t=o
  %1 contagiado por cada millón de habitantes
  
  %ponemos en la  consola:
  %I0=1/1000000
  %S0= 1-I0
  %E0=0
  %L0=0
  %G0=0
  %H0=0
  %ICU0=0
  %REC0=0

  %Y0=[S0,E0,I0,L0,G0,H0,ICU0,REC0]

  %a) Hacer una gráfica de la proporción de contagios activos I(t), 
  % para t en los primeros 100 días. Al hacer la gráfica, usar los 
  % valores dados de las constantes del modelo: R0=2.83, etc.
  %[t,YY] = R_K_n(0, 100, Y0)
 % Graficar la proporción de contagios activos I(t)
figure;
plot(t, YY(:, 3), 'LineWidth', 2);
xlabel('Tiempo (días)');
ylabel('Proporción de contagios activos');
title('Dinámica de contagios activos de COVID-19 en CDMX');
grid on;

  %b) En la CDMX se reportaron los primeros contagios el 28 de febrero de 2020 
  % (año bisiesto). Según la gráfica de a), y tomando como t=0 el 28 de febrero
  % de 2020, ¿en qué fecha predijo el modelo que ocurriría el pico de la pandemia?
  
  %Encontrar el valor máximo de la proporción de contagios activos
[max_I, idx_max] = max(YY(:, 3));

% Encontrar la fecha correspondiente al valor máximo
fecha_pico = datetime('28-Feb-2020') + days(t(idx_max));

disp(['La fecha predicha para el pico de la pandemia es: ', datestr(fecha_pico)]);


%c) Basándose en un modelo matemático, la Secretaría de Salud de México anunció
%que el pico de la pandemia en la CDMX (primera ola) ocurriría en cierta fecha. 
%Averigua dicha fecha y compárala con la fecha que encontraste en b)

%Según el reporte epidemiológico de la CDMX, las proyecciones sin
%intervención apuntaban que el pcio de la pandemia en la CDMX seria entre
%mayo y junio, con nuestro modelo, que está basado en Runge-Kutta y datos
%del mismo reporte concordamos que si estuvo en este intervalo de
%mayo-junio 2020

