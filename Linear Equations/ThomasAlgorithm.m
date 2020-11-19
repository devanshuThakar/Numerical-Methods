%This is progrom to solve a tri-digonal (banded Matrix) system of equation by Thomas Algorithm 
n =4;
% F = zeros(n,1);
% G = zeros(n-1,1);
% E = zeros(n-1,1);
X = zeros(n,1);
A = [2.04,-1,0,0; -1,2.04,-1,0; 0,-1,2.04,-1; 0,0,-1,2.04 ];
R = [40.8; 0.8; 0.8; 200.8];
% F = [2.04; 2.04; 2.04; 2.04];
% G = [-1;-1;-1];
% E = [-1;-1;-1];

%Decomposition
for i=1:n-1
    ratio = A(i+1,i)/A(i,i);
    A(i+1,i+1) = A(i+1,i+1) - A(i,i+1)*ratio;
    A(i+1,i) = ratio;
end
%Forward substitution
for i=2:n
    R(i) = R(i) - R(i-1)*A(i,i-1);
end
% Backward Substitution to get X
for i=n:-1:1
    if(i==n)
        X(i) = R(i)/A(i,i);
        continue;
    end
    X(i) = (R(i) - X(i+1)*A(i,i+1))/A(i,i);    
end