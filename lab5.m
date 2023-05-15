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

    xderv=[t0, t1, t2];

    f=fun(xderv);

    %Se crea la matriz de ceros para crear la tabla de diferencias finitas
    b=zeros(3);
    %Se ponen como primera columna el valor de las ordenadas
    b(:,1)=f';

    %Se construye la tabla de diferencias finitas
    for j=2:3
        for i=j:3
            b(i,j)=(b(i,j-1)-b(i-1,j-1));
        end
    end

    %Calculo del dato interpolado ( productoria)
    fprintf('\nLa matriz de diferencias dividas es : ')
    %Se enumeran las filas en la tabla de resultados
    num=[0:1:3-1];
    M=[num' x' b];
    array2table(M,'VariableNames',{'Punto','Xi','f(Xi)','Δƒ[xi]','Δƒ²[xi]'})

    %Construccion del polinomio interpolador con base en el pivote y
    %la tabla de diferencias
    syms s
    %Se inicializa el primer valor para el polinomio junto a la productoria
    polsig=M(1,3);
    produ=1;
    polderv=bderv(2,2);
    prodderv = 1;
    %Se crea el polinomio de grado n-1
    for i=1:(3-1)
        %polinomio interpolado
        produ=produ*(s-i+1)/factorial(i);
        polsig=polsig+produ*M(i+1,i+3);
        
        if i<n-1
        prodderv=prodderv*(s-i+2)/factorial(i);

        
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
    Seval=(t0-xderv(1))/h;
    fprintf('\nEl valor del punto %.3f es : ',t0)
    poldervfinal(t0)

    %Se convierten los valores x1 correspondiente al método de
    %diferencias finitas
    ev=[-0.5:0.1:1.5];
    %Se convierten los valores x1 correspondiente al método de
    %diferencias finitas
    Sgrafica=(ev-x(1))/h;
    Sgraficaderv=(ev-xderv(1))/h;

    fun2 = @(x) -sin(x);
    %Grafica de los datos experimentales y el polinomio obtenido
    polevaluado=polfinal(Sgrafica);
    polevalderv=poldervfinal(Sgraficaderv);
    fplot(fun,[-0.5,1.5],'r')
    grid on
    axis([0 1.2 -1 1])
    hold on
    fplot(fun2,[-0.5,1.5],'r')
    hold on
    plot(ev,polevaluado,'b')
    hold on
    plot(ev,polevalderv, "g")
    hold on
    plot(t0,polfinal(Seval),'r.','markersize',12)
    hold on
    xlabel('Rad')
    ylabel('F(x)')
    title('Derivación')
    legend('Datos reales','Interpolación','Derivada','Location','southwest')


%----PUNTOS 2 Y 3 (integración)--------
elseif(opc == 2)
    met = input("Ingrese 1 si quiere usar el método del trapecio o 2 si quiere usar el método de Simpson: ");
    
    %----------Punto 2, metodo del trapecio-----------
    if(met == 1)

    %----------Punto 3, metodo de Simpson-----------
    elseif(met == 2)

    end
end
