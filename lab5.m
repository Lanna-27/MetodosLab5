clear
clc
clf('reset')

%Pide la función incial
fun = input("Ingresa la función a evaluar en formato simbolico: ");
%Dar la opción de elegir derivacion o integracion
opc = input("Elige 1 si quieres hacer derivación ó 2 si quieres hacer integración: ");

%-------Punto 1, derivación polinomio interpolador de Newton-------
if(opc == 1)
    %punto 1
    %Derivada del polinomio interpolador de Newton de grado N
    %para aproximar numéricamente f'(t0)
    %Se piden los parametros necesarios
    t0=input('Ingrese el valor en el que desea evaluar la derivada: ');
    h=input('Ingrese el valor de h: ');

    t1 = t0+h;
    t2 = t0-h; 

    % x=[t2, t0, t1]; 

    % f = fun(x);

    x=[t0, t1, t2];

    n=length(x);

    f=fun(x);

    %Se crea la matriz de ceros para crear la tabla de diferencias finitas
    A=zeros(n,n);
    %Se ponen como primera columna el valor de las ordenadas
    A(1:n,1)=f;
    %Se calcula la diferencia entre dos puntos consecutivos de las abscisas
    %h=x1(2)-x1(1);

    %Se construye la tabla de diferencias divididas
    for J=2:n
        J1=J-1;
        for K=1:(n-J1)
            A(K,J)=(A(K+1,J1)-A(K,J1))/(x(K+J1)-x(K));
        end
    end
    A

    %Calculo del dato interpolado ( productoria)
    fprintf('\nLa matriz de diferencias dividas es : ')
    %Se enumeran las filas en la tabla de resultados
    num=[0:1:n-1];
    M=[num' x' A];
    array2table(M,'VariableNames',{'Punto','Xi','f(Xi)','Δƒ[xi]','Δƒ²[xi]'})

    coef = A(1,1:n);
    X=[x(1):0.1:x(n)];
    nd=length(X);
    %Construccion del polinomio interpolador con base en el pivote y
    %la tabla de diferencias
    syms s
    %Se inicializa el primer valor para el polinomio junto a la productoria
    polsig=coef(1);
    produ=1;
    polderv=coef(2);
    prodderv = 1;
    P = zeros(nd,0);
    %Se crea el polinomio de grado n-1
    for i=1:n-1
        %polinomio interpolado
        produ=produ*(s-x(i));
        polsig=polsig+produ*coef(i+1);
        
        if i>1
            prodderv=prodderv*(s-x(i));
            polderv=polderv+prodderv*coef(i+1);
            
        end    
    end

    %Se imprime el polinomio construido
    polsig
    polderv
    polexp=expand(polsig);
    disp('El polinomio de interpolacion de newton es: ')
    pretty(polexp)
    disp('Y su derivada es: ')
    polder=expand(polderv);
    pretty(polder)
    %Se convierte el polinomio construido en una función en linea para evaluar
    %los diferentes valores de s
    polfinal=inline(polexp);
    poldervfinal=inline(polder); 
    %Calculo del valor 'a' a interpolar
    Seval=t0/h;
    fprintf('\nEl valor del punto %.3f es : ',t0)
    poldervfinal(t0)

    %Se convierten los valores x1 correspondiente al método de
    %diferencias finitas
    ev=[-0.5:0.1:1.5];
    %Se convierten los valores x1 correspondiente al método de
    %diferencias finitas
    Sgrafica=(ev-x(1))/h;

    fun2 = @(x) -sin(x);
    %Grafica de los datos experimentales y el polinomio obtenido
    polevaluado=polfinal(ev);
    polevalderv=poldervfinal(ev);
    fplot(fun,[-0.5,1.5],'r')
    grid on
    axis([0 1.2 -1 1])
    hold on
    fplot(fun2,[-0.5,1.5],'m')
    hold on
    plot(ev,polevaluado,'b')
    hold on
    plot(ev,polevalderv, "g")
    hold on
    plot(t0,poldervfinal(t0),'k.','markersize',12)
    hold on
    xlabel('Rad')
    ylabel('F(x)')
    title('Derivación')
    legend('Datos reales', 'Derivada real','Interpolación','Derivada','Location','southwest')


%----PUNTOS 2 Y 3 (integración)--------
elseif(opc == 2)
    met = input("Ingrese 1 si quiere usar el método del trapecio o 2 si quiere usar el método de Simpson: ");
    
    %----------Punto 2, metodo del trapecio-----------
    if(met == 1)
        %Se solicitan los datos necesarios para resolver el problema

f=input('Ingrese la función f(x)='); % f=@(x) 
a=input('Ingrese el límite inferior de la integral:');
b=input('Ingrese el límite superior de la integral:');
n=input('Ingrese el número de nodos:');

m=n-1;
h=(b-a)/m; %Se calcula el parámetro h
sum=0; %Se separa el término de la sumatoria de la ecuación general y se inicializa

%Se evaluan los términos de la sumatoria desde el nodo x1 hasta x_(m-1)
for k=1:1:m-1
    x(k)=a+k*h; %Se calculan los xk
    y(k)=f(x(k));
    sum=sum+y(k);
end

M=zeros(n,2);


%Se aplica la fórmula del método del trapecio compuesto: (h/2)*((f0+fm)+2*(f1+f2+f3+...+fm-1))
T=double(h/2*(f(a)+f(b)+2*sum));
%Se calcula el valor real solucionando la integral definida
%R=double(int(f,a,b));
%Se imprime la solución aproximada y la solución analítica
fprintf('El área aproximada bajo la curva es: %10.15f',T);
%fprintf('El real bajo la curva es: %10.15f',R);

    %----------Punto 3, metodo de Simpson-----------
    elseif(met == 2)
        %Se solicitan los datos necesarios para resolver el problema

f=input('Ingrese la función f(x)='); % f=@(x) 
a=input('Ingrese el límite inferior de la integral:');
b=input('Ingrese el límite superior de la integral:');
n=input('Ingrese el número de nodos:');

m=(n-1);
h=(b-a)/m; %Se calcula el parámetro h

%Se separan los términos de la sumatorias de los términos pares e impares de la ecuación general y se
%inicializan
si=0;
sp=0; 
%Se evaluan los términos de las sumatorias 
for k=1:1:m-1
    x(k)=a+k*h; %Se calculan los xk
    y(k)=f(x(k));
    if rem(k,2)==1 %Se separan los términos impares y se realiza la sumatoria
        si=si+y(k);
    else
        sp=sp+y(k);%sumatoria términos pares
    end
end
%Se aplica la fórmula del método de simpson compuesto:
%(h/3)*((f0+fm)+4*(f1+f3+f5+...+términos impares)+2*(f2+f4+f6...términos pares))
S=double(h/3*(f(a)+f(b)+4*si+2*sp));
%Se calcula el valor real solucionando la integral definida
%R=double(int(f,a,b));
%Se imprime la solución aproximada y la solución analítica
fprintf('El área aproximada bajo la curva es: %10.15f',S);
%fprintf('El real bajo la curva es: %10.15f',R);

%Se grafica el área bajo la curva en el intervalo [a,b]
abcisas=a:0.01:b;
fun=subs(f,abcisas);
area(abcisas,fun)
xlabel ('x')
ylabel('y')
title('Área bajo la curva')
       

    end
end
