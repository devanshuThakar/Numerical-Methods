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
x(1) = Ll;
x((Lu-Ll)/dx + 1)=Lu;
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
    %Temperature is to evaluated at all locations
    for j=2:length(x)-1
        x(j) = x(j-1) + dx;
        T(i,j) = T(i-1,j) + lambda*(T(i-1, j+1) - 2*T(i-1, j) + T(i-1, j-1));
    end
end

%Loop to plot temperature
hold on
Legend = string(zeros(floor(((tu-tl)/dt + 1)/5)));
for i=1:5:(tu-tl)/dt + 1
    plot(x,T(i,:))
    Legend(floor(i/5)+1) = ("t=" + t(i));
end
title("Plot of Temperature(T) vs. distance (x)");
xlabel("x(cm)");
ylabel("Temperature (C)");
legend(Legend);
hold off