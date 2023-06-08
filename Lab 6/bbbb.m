clear;
clc;
clf('reset');
format long;
syms t u v;

% Method 1: Linear Shooting Method
disp('Method 1: Linear Shooting Method');
disp('-------------------------------');

op = input('Enter 1 for the Linear Shooting Method or 2 for the Finite Differences Method: ');
p = input('Enter the function P(t): ');
q = input('Enter the function q(t): ');
r = input('Enter the function r(t): ');
sol = input('If you know the exact solution x(t), enter it here, or press enter: ');
intervalo = input('Enter the working interval [t0, tf]: ');
hin = input('Enter the step size h: ');
alpha = input('Enter the value of alpha: ');
beta = input('Enter the value of beta: ');

if op == 1
    % Linear Shooting Method
   
    disp('Running Linear Shooting Method...');
   
    hlong = length(hin);
   
    for n = 1:hlong
        h = hin(n);
        tinf = intervalo(1);
        tsup = intervalo(2);
        u1 = alpha;
        u2 = 0;
        v1 = 0;
        v2 = 1;
        M = (tsup - tinf) / h;
        U = zeros(1, M + 1);
        V = zeros(1, M + 1);

        % Solve the first ordinary differential equation for U(t)
        % Solve the second ordinary differential equation for V(t)
        fprintf('\n');
        fprintf('For h = %f\n', h);
        fprintf('   t             U            V\n');
        fprintf('--------------------------------\n');

        for i = 1:M + 1
            fprintf('%f   %f   %f\n', T(i), Y(i), Z(i));
            fprintf('\n');
        end

        fprintf('\n');
        fprintf('Error of the Method\n');
        fprintf('-------------------\n');

        for i = 1:M + 1
            errorU = abs(Y(i) - subs(sol));
            errorV = abs(Z(i) - subs(sol));
            fprintf('%f     %f\n', T(i), errorU);
            fprintf('%f     %f\n', T(i), errorV);
            fprintf('\n');
        end
    end

elseif op == 2
    % Method 2: Finite Differences Method
   
    disp('Method 2: Finite Differences Method');
    disp('----------------------------------');
   
    syms x;
    x = sym('x(t)');
    p = sym('p(t)');
    q = sym('q(t)');
    r = sym('r(t)');
   
    eqn = diff(p * diff(x)) + q * diff(x) + r == 0;
    sol = dsolve(eqn);
   
    disp('Differential Equation:');
    disp(eqn);
    disp('Solution:');
    disp(sol);
   
    P = matlabFunction(p);
    Q = matlabFunction(q);
    R = matlabFunction(r);
    SOL = matlabFunction(sol);
   
    a = intervalo(1);
    b = intervalo(2);
    N = (b - a) / hin;
    x = zeros(1, N + 1);
    y = zeros(1, N + 1);
    t = zeros(1, N + 1);
    t(1) = a;
    x(1) = alpha;
    y(1) = beta;

    fprintf('\n');
    fprintf('  t             x             y\n');
    fprintf('--------------------------------\n');

    for i = 1:N
        t(i + 1) = a + i * hin;
        x(i + 1) = (2 * hin^2 * R(t(i)) + (1 - hin * P(t(i))) * x(i) - hin^2 * Q(t(i)) * y(i)) / (1 + hin * P(t(i)));
        y(i + 1) = y(i) + hin * x(i + 1);
        fprintf('%f   %f   %f\n', t(i + 1), x(i + 1), y(i + 1));
        fprintf('\n');
    end

    figure(1);
    plot(t, x);
    hold on;
    plot(t, y);
   
else
    disp('Invalid input. Please enter either 1 or 2 for the method selection.');
end