%Primero se limpia y elimina la ventana de comandos y el área de trabajo
clear
clc
clf('reset')

format long
syms t

%Con el fin de realizar solo un programa, se da la opción de escoger el 
%método del disparo lineal y el método de las diferencias finitas 
op=input('Ingrese 1 para el método del disparo lineal y 2 para el método de las diferencias finitas: ');

P=input('ingrese la función p(t) en formato anónimo: ');
Q=input('ingrese la función q(t) en formato anónimo: ');
R=input('ingrese la función r(t) en formato anónimo: ');
SOL=input('en caso de conocer la solución x(t) ingresela en formato anónimo: ');
intervalo=input('ingrese el intervalo de trabajo: ');
hin=input('ingrese el tamaño de paso h: ');
alpha=input('ingrese el valor de alpha: ');
beta=input('ingrese el valor de beta: ');
hlong=length(hin);

if(op == 1)
    %En caso de que el valor suministrado sea igual a 1 se procede con el
    %método del disparo lineal de Runge-Kutta de orden N=4
    %------Entradas-------
    %@(t) (2*t)/(1+t^2)
    %@(t) (-2)/(1+t^2)
    %1
    %@(t) 1.25+0.4860896526*t-(2.25*t.^2)+(2*t.*atan(t))+(1/2).*(t.^2-1).*(log(1+t.^2))
    %[0 4]
    %[0.2 0.1]
    %alpha=1.25
    %beta=-0.95

    disp('Método 1: Método del disparo lineal');
    disp('-------------------------------');

    for n=1:hlong

        Matriz = [];
        Matriz2 = [];

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

        T=intervalo(1):hin(j):intervalo(2);

        U(1,1)=alpha;
        Y=zeros(1,M+1);
        %y=u'
        Y(1,1)=u2;

        for i=2:M+1
            
            y=Y(1,i-1);
            u=U(1,i-1);
            t=T(i-1);
            
            g1=double(subs(P,t)*y+subs(Q,t)*u+subs(R));
            f1=double(y);
            
            y=Y(1,i-1)+(h/2)*g1;
            u=U(1,i-1)+(h/2)*f1;
            t=T(1,i-1)+(h/2);
            
            g2=double(subs(P,t)*y+subs(Q,t)*u+subs(R));
            f2=double(y);
            
            y=Y(1,i-1)+(h/2)*g2;
            u=U(1,i-1)+(h/2)*f2;
            t=T(1,i-1)+(h/2); 
            
            g3=double(subs(P,t)*y+subs(Q,t)*u+subs(R));
            f3=double(y);
            
            y=Y(1,i-1)+(h/1)*g3;
            u=U(1,i-1)+(h/1)*f3;
            t=T(1,i-1)+(h/1); 
            
            g4=double(subs(P,t)*y+subs(Q,t)*u+subs(R));
            f4=double(y);
            
            U(1,i)=U(1,i-1)+(h/6)*(f1+2*f2+2*f3+f4);
            Y(1,i)=Y(1,i-1)+(h/6)*(g1+2*g2+2*g3+g4);
            
        end    

        U;
        Y;

        %Solución segunda ecuación diferencial ordinaria de forma v(t)

        V(1,1)=v1;
        Z=zeros(1,M+1);
        Z(1,1)=v2;

        for j=2:M+1
            
            z=Z(1,j-1);
            v=V(1,j-1);
            t=T(1,j-1);
            
            g1=double(subs(P,t)*z+subs(Q,t)*v);
            f1=double(z);
            
            z=Z(1,j-1)+(h/2)*g1;
            v=V(1,j-1)+(h/2)*f1;
            t=T(1,j-1)+(h/2);
            
            g2=double(subs(P,t)*z+subs(Q,t)*v);
            f2=double(z);
            
            z=Z(1,j-1)+(h/2)*g2;
            v=V(1,j-1)+(h/2)*f2;
            t=T(1,j-1)+(h/2); 
            
            g3=double(subs(P,t)*z+subs(Q,t)*v);
            f3=double(z);
            
            z=Z(1,j-1)+(h/1)*g3;
            v=V(1,j-1)+(h/1)*f3;
            t=T(1,j-1)+(h/1); 
            
            g4=double(subs(P,t)*z+subs(Q,t)*v);
            f4=double(z);
            
            V(1,j)=V(1,j-1)+(h/6)*(f1+2*f2+2*f3+f4);
            Z(1,j)=Z(1,j-1)+(h/6)*(g1+2*g2+2*g3+g4);
            
        end
        V;
        Z;

        X=zeros(1,M+1);

        for k=1:M+1
            X(1,k)=double(U(1,k)+((beta-U(1,M+1))/V(1,M+1))*V(1,k));
        end
        
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
    Xexacto1=subs(SOL,t);
    t=Matriz2(:,1);
    Xexacto2=subs(SOL,t);
    error1=Matriz(:,2)-Xexacto1;
    error2=Matriz2(:,2)-Xexacto2;
    R1=[Matriz Xexacto1 error1];
    R2=[Matriz2 Xexacto2 error2];

    %Impresión de los resultados 
    if(hlong == 1)
        fprintf('\nLos resultados con h=0.2 son\n:')
        array2table(R1,'VariableNames',{'tj','xj','x exacto','error'})
    elseif(hlong == 2)
        fprintf('\nLos resultados con h=0.2 son\n:')
        array2table(R1,'VariableNames',{'tj','xj','x exacto','error'})

        fprintf('\nLos resultados con h=0.1 son\n:')
        array2table(R2,'VariableNames',{'tj','xj','x exacto','error'})
    end

    T2long=length(Matriz2(:,1));
    X2=Matriz2(:,2);
    for d=3:2:T2long-2
        X2F((d-1)/2,1)=X2(d);
        errorF2((d-1)/2,1)=error2(d);
    end
    X2F=[alpha;X2F;beta];
    errorF2=[0;errorF2;0];

    fprintf("Aproximaciones numericas de x''(t)\n")
    if(hlong == 1)
        Mfinal=[Matriz(:,1) Matriz(:,2) Xexacto1 error1];
        array2table(Mfinal,'VariableNames',{'tj','xj(h=0.2)','x exacto','error(h=0.2)'})
    elseif(hlong == 2)
        Mfinal=[Matriz(:,1) Matriz(:,2) X2F Xexacto1 error1 errorF2];
        array2table(Mfinal,'VariableNames',{'tj','xj(h=0.2)','xj(h=0.1)','x exacto','error(h=0.2)','error(h=0.1)'})
    end

    subplot(1,2,1)
    plot(T,X,'r')
    hold on
    plot(T,W,'g')
    hold on
    plot(T,U,'b')
    grid on
    hold on
    plot(T,Xexacto2,'co')
    hold on
    xlabel('t')
    ylabel('y''(t)')
    title('Disparo lineal aproximaciones')
    legend('X(t)', 'error','U(t)','X(t) exacto','Location','southwest')

    f2=figure;
    subplot(1,2,1)
    plot(T,Xexacto2,'c')
    grid on
    hold on
    xlabel('t')
    ylabel('y''(t)')
    title('Disparo lineal resultado real')
    legend('X(t) exacto','Location','southwest')

%Método por diferencias finitas
elseif(op == 2)
    %En caso de que el valor suministrado sea igual a 2 se procede con el
    %método de diferencias finitas de orden O(h^2)
    %x(t)=1.25+0.4860896526*t-(2.25*t.^2)+(2*t.*atan(t))+(1/2).*(t.^2-1).*(log(1+t.^2))
    %p(t)=(2*t)/(1+t^2)
    %q(t)=(-2)/(1+t^2)
    %r(t)=1
    %hin=[0.2 0.1 0.05 0.025]
    %intervalo = [0 4]
    %alpha=1.25   
    %beta=-0.95

    disp('Metodo 2: Metodo de diferencias finitas');
    disp('----------------------------------');

    %Inicializamos hmin y hmax con el h mas pequeño y mas grande respectvamente
    hmin=min(hin);
    hmax=max(hin);

    %Planteamos las ecuaciones
    for j=1:length(hin)
    T=intervalo(1):hin(j):intervalo(2);
    n=length(T)-1;
    for i=1:n-1
        if i==1
        B(1,1) = double( -(hin(j)^2)*subs(R,T(i+1))+ ((hin(j)/2)*subs(P,T(i+1))+1)*alpha);
        elseif i==n-1
        B(i,1) = double( -(hin(j)^2)*subs(R,T(i+1)) + ((-hin(j)/2)*subs(P,T(i+1))+1)*beta);
        else
        B(i,1) = double( -(hin(j)^2)*subs(R,T(i+1)) );
        end
    end

    %Crea matriz tridiagonal de las diapositivas para representar las ecuaciones para p y q
    diagP=[];
    diagI=[];
    diagS=[];
    diagP(1,1)=0;
    for i=1:n-1
        diagP(i,1) = double( 2+hin(j)^2*subs(Q,T(i+1)));
        if i~=1
        diagI(i-1,1) = double( (-hin(j)/2)*subs(P,T(i+1))-1);
        end
        if i~=n-1
        diagS(i+1,1) = double((hin(j)/2)*subs(P,T(i+1))-1);
        end
    end
    diagI(n-1,1)=0;
    A=spdiags([diagI, diagP, diagS],[-1;0;1],n-1,n-1);

    %Resolvemos el sistema de ecuaciones
    X = A\B;

    %Teniendo en cuenta la cantidad de h que haya ingresado el usuario se arma la tabla de aproximaciones
    for i=intervalo(1):hin(j):intervalo(2)
        if(i==intervalo(1))
        tab(1,j)=alpha;
        elseif(i==intervalo(2))
        tab((i-intervalo(1))/hmax+1,j)=beta;
        elseif mod((i-intervalo(1)),hmax)==0
        tab(int16((i-intervalo(1))/hmax)+1,j)=X(int16(i/hin(j)),1);
        end
    end
    end

    %Damos formato a la tabla para que quede igual a la del ejemplo
    tj= intervalo(1):hmax:intervalo(2);
    tj=tj';
    fun = sym(1.25+0.4860896526*t-(2.25*t.^2)+(2*t.*atan(t))+(1/2).*(t.^2-1).*(log(1+t.^2)));
    xj = [];
    for i=1:length(tj)
    xj(i,1)=double(subs(fun,tj(i,1)));
    end
    tabla = horzcat(tj,tab,xj);

    %Imprimos el encabezado de la tabla y despues mostramos la tabla de aproxmaciones
    fprintf("Aproximaciones numericas de x''(t)\n")
    if length(hin)==1
    fprintf("   tj\t\t       h=0.2\t\t   exacto\n")
    elseif length(hin)==2
    fprintf("   tj\t\t       h=0.2\t\t   h=0.1\t       exacto\n")
    elseif length(hin)==3
    fprintf("   tj\t\t       h=0.2\t\t   h=0.1\t       h=0.05\t\t   exacto\n")
    else
    fprintf("   tj\t\t       h=0.2\t\t   h=0.1\t       h=0.05\t\t   h=0.025\t       exacto\n")
    end
    disp(tabla)
    fprintf("\n")

    %Creamos la matriz de errores
    errores=[];
    for i=2:length(tabla(1,:))-1
    errores(:,i-1)=tabla(:,length(tabla(1,:)))-tabla(:,i);
    end
    errores = horzcat(tj,errores);

    %Imprimos el encabezado de la tabla y despues mostramos la tabla de errores
    fprintf("Errores de las aproximaciones numericas\n")
    if length(hin)==1
    fprintf("   tj\t\t       x(tj)-xj1\n")
    elseif length(hin)==2
    fprintf("   tj\t\t       x(tj)-xj1\t   x(tj)-xj1\n")
    elseif length(hin)==3
    fprintf("   tj\t\t       x(tj)-xj1\t   x(tj)-xj2\t       x(tj)-xj3\n")
    else
    fprintf("   tj\t\t       x(tj)-xj1\t   x(tj)-xj2\t       x(tj)-xj3\t   x(tj)-xj4\n")
    end
    disp(errores)
    fprintf("\n")

    %Creamos un arreglo con las leyendas de los graficos
    legends = {};
    for i=1:length(hin)
    legends(i) = cellstr(num2str(hin(i)));
    end
    legends(length(hin)+1)=cellstr("x(t)");

    %Imprimos todas las aproximaciones generadas en una misma grafica
    hold on
    for i=2:length(tabla(1,:))
    plot(tabla(:,1),tabla(:,i));
    end
    hold off
    legend(legends);
    title("Aproximación numérica a la solución de x''(t)")
    xlabel("tj")
    ylabel("x''(t)")
    legends = legends(1:end-1);

    fig2=figure;
    %Limpiamos el grafico anterior e Imprimos todas lo errores generados en una misma grafica
    hold on
    for i=2:length(errores(1,:))
    plot(errores(:,1),errores(:,i));
    end
    hold off
    title("Errores de las aproximaciones numéricas");
    xlabel("tj");
    ylabel("x(tj)");
    legend(legends);
    
end