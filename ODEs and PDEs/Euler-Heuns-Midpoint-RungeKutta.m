clear all;
clc;
h = 1; %time step
t = zeros((50-0)/h + 1, 1);

%Euler's method
x_euler = zeros((50-0)/h + 1, 1);
v_euler = zeros((50-0)/h + 1, 1);
% Loop for Euler method
for i=2:(50-0)/h+1
    f_eval = 0;
    t(i) = t(i-1)  + h;
    v_euler(i) = v_euler(i-1) + (func(t(i-1), v_euler(i-1), x_euler(i-1)))*h;
    f_eval = f_eval + 1;
    x_euler(i) = x_euler(i-1) + v_euler(i-1)*h;
end

%Heun's method
x_heun = zeros((50-0)/h + 1, 1);
v_heun = zeros((50-0)/h + 1, 1);
t = zeros((50-0)/h + 1, 1);
%Loop for Heun's method
for i=2:(50-0)/h + 1
    t(i) = t(i-1)  + h;
    %Predictor
    vp = v_heun(i-1) + func(t(i-1), v_heun(i-1), x_heun(i-1))*h;
    xp = x_heun(i-1) + v_heun(i-1)*h;
    %Corrector
    v_heun(i) = v_heun(i-1) +(func(t(i-1), v_heun(i-1), x_heun(i-1)) + func(t(i-1), vp, xp))*h/2;
    x_heun(i) = x_heun(i-1) + (v_heun(i-1) + vp)*h/2; 
end

%Midpoint-method
v_midpoint = zeros((50-0)/h + 1, 1);
x_midpoint = zeros((50-0)/h + 1, 1);
%Loop for mid-point method
for i=2:(50-0)/h+1
    t(i) = t(i-1)  + h;
    %Mid-point
    vmp = v_midpoint(i-1) + func(t(i-1), v_midpoint(i-1), x_midpoint(i-1))*h/2;
    xmp = x_midpoint(i-1) + v_midpoint(i-1)*h/2;
    %Next-Value
    v_midpoint(i) = v_midpoint(i-1) + func(t(i-1) + h/2, vmp, xmp)*h;
    x_midpoint(i) = x_midpoint(i-1) + v_midpoint(i)*h;
end

%Fourth order Range-Kutta method
x_RK4 = zeros((50-0)/h + 1, 1);
v_RK4 = zeros((50-0)/h + 1, 1);
%Loop for fourth order Range-Kutta method
for i=2:(50-0)/h+1
    k1v = func(t(i-1), v_RK4(i-1), x_RK4(i-1));
    k1x = v_RK4(i-1);
    k2v =  func(t(i-1)+ h/2, v_RK4(i-1) + (k1v*h)/2, x_RK4(i-1) + (k1x*h)/2);
    % To obatin V at  (t(i-1)+ h/2, v_RK4(i-1) + (k1v*h)/2, x_RK4(i-1) +  (k1x*h)/2), Euler's linear interpolation has been used
    %V_k2x = v_RK4(i-1) + k2v*h/2;
    % % Check this again #############
    k2x = v_RK4(i-1) + k2v*h/2;
    
    k3v = func(t(i-1) + h/2, v_RK4(i-1) + (k2v*h)/2, x_RK4(i-1) + (k2x*h)/2);
    k3x = v_RK4(i-1) + k3v*h/2;
    k4v = func(t(i-1) + h, v_RK4(i-1) + k3v*h, x_RK4(i-1) + k3x*h);
    k4x = v_RK4(i-1)+ k4v*h;
    
    v_RK4(i) = v_RK4(i-1) + (k1v + 2*k2v + 2*k3v + k4v)*h/6;
    x_RK4(i) = x_RK4(i-1) + (k2x + 2*k2x + 2*k3x + k4x)*h/6;
end

figure(1)
a = plot(t,v_euler, t,v_heun, t, v_midpoint, t, v_RK4);
title('Plot of velocity(m/s) vs. time(d)');
xlabel('Time t(s)');
ylabel('Velocity (m/s)');
legend('Euler', 'Heuns', 'Midpoint', 'Fourth order RK');
grid on

figure(2)
b = plot(t,x_euler, t,x_heun, t, x_midpoint, t, x_RK4);
title('Plot of position(m) vs. time(d)');
xlabel('Time t(s)');
ylabel('Position (m/s)');
legend('Euler', 'Heuns', 'Midpoint', 'Fourth order RK');
grid on

% figure(1)
% a = plot(t,v_euler);
% title('Plot of velocity(m/s) vs. time(d)');
% xlabel('Time t(s)');
% ylabel('Velocity (m/s)');
% legend('Euler with h=0.1');
% grid on
% 
% figure(2)
% b = plot(t,x_euler);
% title('Plot of position(m) vs. time(d)');
% xlabel('Time t(s)');
% ylabel('Position (m/s)');
% legend('Euler with h=0.1');
% grid on


%To evaluate the function
function f = func(t, v, x)
    g = 9.81;
    L = 30;
    m = 68.1;
    cd = 0.25;
    k = 40;
    Y = 8;  %Gamma, cord damping coefficient
    if(x<= L)
        f = g - (sigmoid(v)*cd*(v^2))/m;
    end
    if(x>L)
        f= g - (sigmoid(v)*cd*(v^2))/m - (k*(x-L))/m - (Y*v)/m;
    end
end

function sign = sigmoid(x)
    if(x==0)
        sign = 0;
    end
    if (x~=0)
        sign = abs(x)/x;
    end
end