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
    x1=input('Ingrese el vector de  abscisas(x) equidistantes: ');
    f=input('Ingrese el vector de ordenadas(fx): ');
    while(length(vec_x) ~= length(vec_y))
        printf("Los vectores ingresados no son de la misma longitud\n")
        vec_x = input("Vector de X: ");
        vec_y = input("Vector de Y: ");
    end

    a=input('Ingrese un valor para evaluar la derivada: ');
    % pivote=input('Ingrese la casilla del dato que usara como pivote en la construcción del polinomio: ');
    % grado=input('¿De que grado quiere que sea el polinomio aproximado? : ');

    %Se calcula la magnitud del vector de abscisas
    n=length(x1);

    %Se rectifica que las abscisas sean equidistantes
    for z=1:n
        if x1(z+1)-x1(z)~=x1(z+2)-x1(z+1)
            disp('Los puntos de las abscisas no son equidistantes')
            x1=input('Ingrese el vector de  abscisas(x) equidistantes: ');
        else
            break
        end
    end

    %Se crea la matriz de ceros para crear la tabla de diferencias finitas
    b=zeros(n);
    %Se ponen como primera columna el valor de las ordenadas
    b(:,1)=f';
    %Se calcula la diferencia entre dos puntos consecutivos de las abscisas
    h=x1(2)-x1(1);

    %Se construye la tabla de diferencias finitas
    for j=2:n
    for i=j:n
    b(i,j)=(b(i,j-1)-b(i-1,j-1));
    end
    end

    %Calculo del dato interpolado ( productoria)
    fprintf('\nLa matriz de diferencias dividas es : ')
    %Se enumeran las filas en la tabla de resultados
    num=[0:1:n-1];
    M=[num' x1' b];
    array2table(M,'VariableNames',{'Punto','Xi','f(Xi)','Δƒ[xi]','Δƒ²[xi]','Δƒ³[xi]','Δƒ4[xi]','Δƒ5[xi]'})

    %Construccion del polinomio interpolador con base en el pivote y
    %la tabla de diferencias
    syms s
    %Se inicializa el primer valor para el polinomio junto a la productoria
    polsig=M(pivote,3);
    produ=1;
    %Se crea el polinomio del grado deseado por el usuario (maximo n-1)
    for i=1:grado
    produ=produ*(s-i+1)/factorial(i);
    polsig=polsig+produ*M(i+pivote,i+3);
    end

    %Se imprime el polinomio construido
    polsig;
    polsig=expand(polsig);
    disp('El polinomio de interpolacion de newton es: ')
    pretty(polsig)
    %Se convierte el polinomio construido en una función en linea para evaluar
    %los diferentes valores de s
    polfinal=inline(polsig);
    %Calculo del valor 'a' a interpolar
    Seval=(a-x1(pivote))/h;
    fprintf('\nEl valor del punto %.3f es : ',a)
    polfinal(Seval)

    %Se convierten los valores x1 correspondiente al método de
    %diferencias finitas
    Sgrafica=(x1-x1(pivote))/h;

    %Grafica de los datos experimentales y el polinomio obtenido
    polevaluado=polfinal(Sgrafica);
    plot(x1,f,'r')
    grid on
    axis([x1(1) x1(n) f(1) f(n)])
    hold on
    plot(x1,polevaluado,'b')
    hold on
    plot(a,polfinal(Seval),'r.','markersize',12)
    xlabel('T(°F)')
    ylabel('Pvapor(lb/in²)')
    title('Diferencias finitas')
    legend('Datos reales','Interpolación','Location','northwest')


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


%Se aplica la fórmula del método del trapecio: (h/2)*((f0+fm)+2*(f1+f2+f3+...+fm-1))
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

%Se separan los términos de la sumatorias de los térinos pares e impares de la ecuación general y se
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


#punto 2
#punto 3
