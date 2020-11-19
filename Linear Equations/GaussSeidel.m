A = [3,-0.1,-0.2; 0.1,7,-0.3; 0.3,-0.2,10];
B = [7.85;-19.3;71.4];
Xo = [0;0;0];               %Initial Guess
n = length(B);
Xnew = zeros(n,1);
Ea = zeros(n,1) + 100;      %To hold the approximate errors

converge = false;
while ~converge
    converge = true;
    for i=1:n
        if(Ea(i)>1)        
            %Loop to get SUM
            SUM = 0;
            for j=1:n
                SUM = SUM + A(i,j)*Xo(j);
            end
            SUM = SUM - A(i,i)*Xo(i);
            Xnew(i) = (B(i) - SUM)/A(i,i);
            Ea(i) = abs((Xnew(i) - Xo(i))/Xnew(i))*100;
            Xo(i) = Xnew(i);
            converge = false;
        end
    end
end