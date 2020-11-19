%Problem Parameters
rho = 2.7;
C = 0.2174;
k = 0.49;
alpha = k/(rho*C);

Ll = 0;
Lu = 10;
tl = 0;
tu = 10;
dt = 0.2;
dx = 1;

lambda = (alpha*dt)/(dx^2);

%Arrays for explicit method
x = zeros((Lu-Ll)/dx + 1, 1);
for j=2:length(x)
    x(j) = x(j-1) + dx;
end
t = zeros((tu-tl)/dt + 1, 1);
%Here spatial variation is in rows, and time variation in columns
T = zeros((tu-tl)/dt + 1, (Lu-Ll)/dx + 1);
%Setting the BCs
T(:, 1) = 100;
T(:, (Lu-Ll)/dx + 1) = 50;
%Loop for explicit method
%Temperature is known at t=0
for i=2:length(t)
    t(i) = t(i-1) + dt;
    %For every time, a system of linear equations are need to be solved using Thomas algorithm
    %Loop to set R
    R = zeros((Lu-Ll)/dx-1,1);
    R(1) = 2*lambda*T(i-1,1) + 2*(1-lambda)*T(i-1,2) + lambda*T(i-1,3);
    for j=2:length(R) - 1
        R(j) = lambda*T(i-1,j-1+1) + 2*(1-lambda)*T(i-1,j+1) + lambda*T(i-1,j+1+1);
    end
    R(length(R)) = lambda*T(i-1,(Lu-Ll)/dx-1) + (2*(1-lambda)*T(i-1,(Lu-Ll)/dx)) + 2*lambda*T(i-1,(Lu-Ll)/dx + 1);
    
    f = zeros((Lu-Ll)/dx-1,1) + 2*(1+lambda);
    e = zeros((Lu-Ll)/dx-2,1) + (-1*lambda);
    g = zeros((Lu-Ll)/dx-2,1) + (-1*lambda);
    
    T(i,2:(Lu-Ll)/dx) = ThomasAlgorithm(e,f,g,R);
end

%Loop to plot temperature
hold on
Legend = string(zeros(floor(((tu-tl)/dt + 1)/10)));
k = 5; %Variable to control the number of curves on the plot
for i=1:k:(tu-tl)/dt + 1
    plot(x,T(i,:))
    Legend(floor(i/k)+1) = ("t=" + t(i));
end
title("Plot of Temperature(T) vs. distance (x), by Crank-Nicolsen method");
xlabel("x(cm)");
ylabel("Temperature (C)");
legend(Legend);
hold off


%f is the digonal, array. e is the lower band digonal and g is upper one
function Thomas = ThomasAlgorithm(e,f,g,R)
    n = length(f);
    X = zeros(1, n);
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