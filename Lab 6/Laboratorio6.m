%Primero se limpia y elimina la ventana de comandos y el área de trabajo
clear
clc
clf('reset')

format long;
syms x;

%Con el fin de realizar solo un programa, se da la opción de escoger el 
%método del disparo lineal y el método de las diferencias finitas 
M=input('Ingrese 1 para el método del disparo lineal y 2 para el método de las diferencias finitas: ');

%En caso de que el valor suministrado sea igual a 1 se procede con el
%método del disparo lineal de Runge-Kutta de orden N=4
%x(t)=1.25+0.4860896526*t-(2.25*t.^2)+(2*t.*atan(t))+(1/2).*(t.^2-1).*(log(1+t.^2))
%p(t)=(2*t)/(1+t^2)
%q(t)=(-2)/(1+t^2)
%r(t)=1
%h=[0.2 0.1]
%[0 4]
%alpha=1.25
%beta=-0.95

p=input('ingrese la función P(t): ');
q=input('ingrese la función q(t): ');
r=input('ingrese la función r(t): ');
sol=input('en caso de conocer la solución x(t) ingresela: ');
intervalo=input('ingrese el intervalo de trabajo: ');
hin=input('ingrese el tamaño de paso h: ');
alpha=input('ingrese el valor de alpha: ');
beta=input('ingrese el valor de beta: ');
hlong=length(hin);

for n=1:hlong
   
    syms t y u z v

    h=hin(n);    
    tinf=intervalo(1,1);
    tsup=intervalo(1,2);

    u1=alpha;
    u2=0;
    v1=0;
    v2=1;
    M=(tsup-tinf)/h;
    U=zeros(1,M+1);
    V=zeros(1,M+1);

    %Solucion primera ecuación diferencial ordinaria de forma U(t)

    T=zeros(1,M+1);

    for j=1:M+1
        T(1,j)=tinf+(j-1)*h;
    end

    U(1,1)=alpha;
    Y=zeros(1,M+1);
    Y(1,1)=u2;
    yf=y;
    uf=u;

    for i=2:M+1
        
        y=Y(1,i-1);
        u=U(1,i-1);
        t=T(1,i-1);
        
        g1=eval(p)*eval(yf)+eval(q)*eval(uf)+eval(r);
        f1=eval(yf);
        
        y=Y(1,i-1)+(h/2)*g1;
        u=U(1,i-1)+(h/2)*f1;
        t=T(1,i-1)+(h/2);
        
        g2=eval(p)*eval(yf)+eval(q)*eval(uf)+eval(r);
        f2=eval(yf);
        
        y=Y(1,i-1)+(h/2)*g2;
        u=U(1,i-1)+(h/2)*f2;
        t=T(1,i-1)+(h/2); 
        
        g3=eval(p)*eval(yf)+eval(q)*eval(uf)+eval(r);
        f3=eval(yf);
        
        y=Y(1,i-1)+(h/1)*g3;
        u=U(1,i-1)+(h/1)*f3;
        t=T(1,i-1)+(h/1); 
        
        g4=eval(p)*eval(yf)+eval(q)*eval(uf)+eval(r);
        f4=eval(yf);
        
        U(1,i)=U(1,i-1)+(h/6)*(f1+2*f2+2*f3+f4);
        Y(1,i)=Y(1,i-1)+(h/6)*(g1+2*g2+2*g3+g4);
        
    end    

    U;
    Y;

    %Solución segunda ecuación diferencial ordinaria de forma v(t)

    V(1,1)=v1;
    Z=zeros(1,M+1);
    Z(1,1)=v2;
    zf=z;
    vf=v;

    for j=2:M+1
        
        z=Z(1,j-1);
        v=V(1,j-1);
        t=T(1,j-1);
        
        g1=eval(p)*eval(zf)+eval(q)*eval(vf);
        f1=eval(zf);
        
        z=Z(1,j-1)+(h/2)*g1;
        v=V(1,j-1)+(h/2)*f1;
        t=T(1,j-1)+(h/2);
        
        g2=eval(p)*eval(zf)+eval(q)*eval(vf);
        f2=eval(zf);
        
        z=Z(1,j-1)+(h/2)*g2;
        v=V(1,j-1)+(h/2)*f2;
        t=T(1,j-1)+(h/2); 
        
        g3=eval(p)*eval(zf)+eval(q)*eval(vf);
        f3=eval(zf);
        
        z=Z(1,j-1)+(h/1)*g3;
        v=V(1,j-1)+(h/1)*f3;
        t=T(1,j-1)+(h/1); 
        
        g4=eval(p)*eval(zf)+eval(q)*eval(vf);
        f4=eval(zf);
        
        V(1,j)=V(1,j-1)+(h/6)*(f1+2*f2+2*f3+f4);
        Z(1,j)=Z(1,j-1)+(h/6)*(g1+2*g2+2*g3+g4);
        
    end
    V;
    Z;

    X=zeros(1,M+1);

    for k=1:M+1
        X(1,k)=U(1,k)+((beta-U(1,M+1))/V(1,M+1))*V(1,k);
    end

    X';
    W=X-U;

    if n==1
        Matriz=[T' X' W'];
    else
        Matriz2=[T' X' W'];
    end

end

%Creación de la matriz de resultados 
Matriz=[Matriz(:,1) Matriz(:,2)];
Matriz2=[Matriz2(:,1) Matriz2(:,2)];
t=Matriz(:,1);
Xexacto1=eval(sol);
t=Matriz2(:,1);
Xexacto2=eval(sol);
error1=Matriz(:,2)-Xexacto1;
error2=Matriz2(:,2)-Xexacto2;
R1=[Matriz Xexacto1 error1];
R2=[Matriz2 Xexacto2 error2];

%Impresión de los resultados 
fprintf('\nLos resultados con h=0.2 es\n:')
array2table(R1,'VariableNames',{'tj','xj','x exacto','error'})
fprintf('\nLos resultados con h=0.1 es\n:')
array2table(R2,'VariableNames',{'tj','xj','x exacto','error'})

T2long=length(Matriz2(:,1));
X2=Matriz2(:,2);
for d=3:2:T2long-2
    X2F((d-1)/2,1)=X2(d);
    errorF2((d-1)/2,1)=error2(d);
end
X2F=[alpha;X2F;beta];
errorF2=[0;errorF2;0];

Mfinal=[Matriz(:,1) Matriz(:,2) X2F Xexacto1 error1 errorF2];
array2table(Mfinal,'VariableNames',{'tj','xj(h=0.2)','xj(h=0.1)','x exacto','error(h=0.2)','error(h=0.1)'})

f1=figure
subplot(1,2,1)
plot(T,X,'r')
hold on
plot(T,W,'g')
hold on
plot(T,U,'b')
grid on
hold on
xlabel('t')
ylabel('y')
title('Disparo lineal aproximaciones')
legend('X(t)', 'error','U(t)','Location','southwest')

f2=figure
subplot(1,2,1)
plot(T,Xexacto2,'c')
grid on
hold on
xlabel('t')
ylabel('y')
title('Disparo lineal resultado real')
legend('X(t)', 'error','U(t)','Location','southwest')
%por hacer: poner nombres a los ejes de las graficas, hacer que la grafica
%de la solucion real este en otro recuadro puesto se superpone a la
%solucion calculada con matlab y ponerle nombre a cada una de la curvas 
