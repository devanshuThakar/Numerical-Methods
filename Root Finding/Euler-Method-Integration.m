%Parameters of model
Ta = 21;  %Ambient Temperature C
k = 0.017;             

%Dependent variables
Tnum = zeros(11, 1);    %Array to store temperautre calculated numerically
Texcat = zeros(11, 1);  %Array to store temperature calculated anlytically

%Independent variables
t = zeros(11, 1);    

%Setting the intial values 
To = 68;
Tnum(1) = To;
Texcat(1) = To;

for i=2:11
    t(i) = 1+t(i-1);
    Texcat(i) = Ta + 47*exp(-k*t(i));
    Tnum(i) = Tnum(i-1) - (k*(Tnum(i-1)-Ta)*(t(i)-t(i-1)));
end

a = plot(t, Tnum, t, Texcat);
title('Plot of Tnum and Texcat');
xlabel('Time t(min)');
ylabel('Temperature (C)');
legend('Tnum','Texcat');