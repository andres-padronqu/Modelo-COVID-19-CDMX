
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
% Ejemplo (dimension 3): si se tiene un sistema de EDOs 
% de la forma: 
%  x'=t*x-y+z^2
%  y'=x+2*cos(z)-t
%  z'=xy+z,
%  
% es conveniente escribirlo en la forma
% Y'=F(t,Y), donde Y es un vector con 3 componentes
%  Y=(Y(1),Y(2),Y(3)) 
%  En este caso, si W=F(t,Y)
% entonces: 
%     W(1)= Y(1)-Y(2)+Y(3)^2
%     W(2)=Y(1)+2*cos(Y(3))
%     W(3)=Y(1)*Y(2)+Y(3)
% Estas 3 líneas se escriben dentro de la función F
% de abajo.  

 
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
 W(1)= Y(1)-Y(2)+Y(3)^2;
 W(2)=Y(1)+2*cos(Y(3));
 W(3)=Y(1)*Y(2)+Y(3);


