clear all;
clc;
m1 = 40;
m2 = 10;
m3 = 50;

u1 = 0.2;
u2 = 0.5;
u3 = 0.3;
g=9.8;

B = [+m1*g*sin(pi/6) + m1*u1*g*cos(pi/6); +m2*g*sin(pi/6) + m2*u2*g*cos(pi/6); m3*g*sin(pi/3) - m3*u3*g*cos(pi/3)];
A = [1,0,-m1; -1,1,-m2; 0,1,m3];

X = GaussEli(A,B);

function Gauss = GaussEli(A,B)
    %Gauss-Elimination Algorithm
    %Takes two input A(n X n) and B(n X 1)
    n=length(B);
    X = zeros(n,1);
    %Loop to Do Forward Elimination
    for i=1:n                       % This loop is for the iteration over all rows. 
        for j=i+1:n                 % This loop will run on all the rows below the ithe row. This loop is used for pivoting. 
            ratio = A(j,i)/A(i,i);  % This is the ratio, thats need to be multiplied to ith row, to make the below row zero. It is a21/a11 for the first time
            for k=i:n               %Iterating over all columns. This loop is for the subration from all the elements, after mutiplying a21/a11
                A(j,k) = A(j,k) - A(i,k)*ratio;
            end
            B(j) = B(j) - B(i)*ratio;        
        end
    end

    %Loop to do Back Substitution
    for i=n:-1:1        %This is the main loop, ti iterate over all xi's
        SUM = 0;        % This keeps track of Sigma a12x22 + a13x33 and so, on
        for j=i+1:n     %This will update the SUM
            SUM = SUM + (A(i, j)*X(j,1));
        end    
        X(i) = (B(i) - SUM)/A(i,i);
    end
    Gauss = X;
end

