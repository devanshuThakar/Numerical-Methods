%This is progrom solve a Boundary value problem using Finite-difference method
E = 3e04*6.89e06;
w = 4751.4;
L = 3.04799;
I = 800*(0.0254)^4;

dx = 2; %In feet
xl = 0;
xu = 10;

x = zeros((xu-xl)/dx + 1, 1);
for i=2:length(x)
    x(i) = x(i-1) + dx;
end
y = zeros((xu-xl)/dx + 1, 1);
f = zeros((xu-xl)/dx - 1, 1);
g = zeros((xu-xl)/dx - 2, 1);
e = zeros((xu-xl)/dx - 2, 1);

f = f -2;
e=e+1;
g=g+1;
R = zeros((xu-xl)/dx - 1, 1);
for i = 1:(xu-xl)/dx - 1
    R(i) = (0.304799^3)*w*(5*x(i+1) - 0.5*x(i+1)*x(i+1))/(E*I);
end
X_finite=ThomasAlgorithm(e,f,g,R);

%f is the digonal, array. e is the lower band digonal and g is upper one
function Thomas = ThomasAlgorithm(e,f,g,R)
    n = length(f);
    X = zeros(n,1);
    %Decomposition
    for i=2:n
        %index of e,will be i-1
        e(i-1) = e(i-1)/f(i-1);
        %index of g,will be i
        f(i) = f(i) - e(i-1)*g(i-1);
    end
    %Forward substitution
    for i=2:n
        R(i) = R(i) - R(i-1)*e(i-1);
    end
    % Backward Substitution to get X
    X(n) = R(n)/f(n);
    for i=n-1:-1:1
        X(i) = (R(i) - g(i)*X(i+1))/f(i);
    end
    Thomas = X;
end
