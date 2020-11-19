clear all;
clc;
format long;
%Initial Guess
xl = 0.1;
xu = 0.2;
xr_old = 0;
ea = 100;
fxu = func(xu);
fxl = func(xl);

iter = 0;
Table = zeros(7, 5);

while(ea>1)
    iter=iter+1;
    Table(iter, 1) = iter;
    Table(iter, 2) = xl; 
    Table(iter, 3) = xu; 
    xr = (xu+xl)/2;
    Table(iter, 4) = xr; 
    fxr = func(xr);
    if (fxr*fxl < 0)
        fxu = fxr;
        ea = abs((xr-xr_old)/xr)*100;
        xu = xr;
        %ea = abs((xu-xl)/xu)*100;     % Approximate error
        Table(iter, 5) = ea; 
        xr_old= xr;
        continue;
    end 
    if (fxr*fxu < 0)
        fxl = fxr;
        ea = abs((xr-xr_old)/xr)*100;
        xl = xr;
        %ea = abs((xu-xl)/xu)*100;
        Table(iter, 5) = ea;
        xr_old= xr;
        continue;
    end    
    if(fxr == 0)
        xr_old= xr;
        break;
    end
end

display("The caluculated root is " + xr + " ;with iterations = " + iter);
