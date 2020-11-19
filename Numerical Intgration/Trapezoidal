mu = 0.005;                         %SI units
Q = 10e-06;                         %m3/s
rho = 1000;                         %kg/m3
L = 10;                             %cm
x = [0,2,4,5,6,7,10];               %cm
r = [2,1.35,1.34,1.6,1.58,1.42,2];  %mm
n = length(x);
fxr = zeros(length(r), 1);
for i=1:n
    fxr(i) = 1/r(i)^4;
end

%Assuming the pipe to be of constant radius of average radius
ravg = Trapezoidal(x,r)/L;

Integral = Trapezoidal(x, fxr);

Delta_P = (-8*mu*Q*Integral*10^(-2))/(pi*10^(-12));
Delta_P_avg = (-8*mu*Q*10*10^(-2))/(pi*10^(-12)*ravg^4);
Re = (2*Q*rho)/(pi*ravg*10^(-3)*mu);



%Pass arrays of x, and fx of same size
%The function will return the numerical integration of f(x)from x=x1 to x=xn
function Trapezoid = Trapezoidal(x, fx)
    n = length(x);
    SUM = 0;
    for i=1:n-1
        SUM = SUM + (x(i+1)-x(i))*(fx(i+1) + fx(i))/2;
    end
    Trapezoid = SUM;
end
