clear all;
clc;
format long;

xprev = 0.2; %Initial Guess
xtrue = 0.144932848274478; %This is the value of x-obtained from the last iteration of Newton-Raphson method
xnext=0;
ea = 100;
Table = zeros(1, 5);
iter = 0;

while (ea>0.1)
    iter = iter + 1;
    Table(iter, 1) = iter;
    Table(iter, 2) = xprev;  
    xnext = xprev - (func(xprev)/fdash(xnext));
    Table(iter, 3) = xnext;
    ea = abs((xnext-xprev)/xnext)*100;
    Table(iter, 4) = ea;
    et = abs((xnext-xtrue)/xnext)*100;
    Table(iter, 5) = et;
    xprev = xnext;
end
display("The caluculated root is " + xnext + " ;with iterations = " + iter);