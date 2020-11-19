clear all;
clc;
t = [0,10,20,30,35,40,45,50];
c = [10,35,55,52,40,37,32,34];
Q = 4;   %m3/min


%Calculation of M manually
Integration = 0;
i=1;
h1 = 10;
Integration = Integration + (3*h1*(c(i) + (3*c(i+1)) + (3*c(i+2)) + c(i+3)))/8;
i=i+3;
h1 = 5;
Integration = Integration + (h1*(c(i) + (4*c(i+1)) + c(i+2)))/3;
i=i+2;
Integration = Integration + (h1*(c(i) + (4*c(i+1)) + c(i+2)))/3;

M_manual = Q*Integration;
M_function = Q*NumericalIntegration(c,t);
%This is the function to do numerical integration, on passing two arrayd
% Though this is a generic function to do numerical integration using
% Trapezoidal and Simpson (1/3) and (3/8) rule, the drawback of this
% function is that for an interval where there are two possibilities as
% (1) Simpson's 3/8 rule + Trapezoidal rule OR (2) Simpson's 1/3 rule +Simpson's 1/3 rule
% This function selecets the former one i.e. Simpson's 3/8 rule + Trapezoidal rule
% rather than the later one. Thus reducing the accuracy. 
function Integration = NumericalIntegration(fx, x)
    n = length(fx);
    i = 1;
    SUM = 0;
    while(i<n)
        h1 = x(i+1)-x(i);
        h2 = 0;
        h3 = 0;
        if(i<n-1)
            h2 = x(i+2)-x(i+1);
        end
        if (i<n-2)
            h3 = x(i+3) - x(i+2);
        end
        if(h1==h2 && h2==h3 && h3==h1)
            %Using Simpson's 3/8 rule
            SUM = SUM + (3*h1*(fx(i) + (3*fx(i+1)) + (3*fx(i+2)) + fx(i+3)))/8;
            i = i+3;
            continue;
        end
        if (h1 == h2)
            %Using Simpson's 1/3 rule
            SUM = SUM + (h1*(fx(i) + (4*fx(i+1)) + fx(i+2)))/3;
            i = i+2;
            continue;
        end
        %Using trapezoidal rule
        SUM = SUM + (h1*(fx(i)+fx(i+1)))/2;
        i = i+1;
    end
    Integration = SUM;
end