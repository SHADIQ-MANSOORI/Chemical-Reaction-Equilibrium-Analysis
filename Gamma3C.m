%This is a function file for Gamma for Wilson 3 componenet system for
%wilson

function g = Gamma3C(Del,X1,X2,X3)
    for i = 1:3
        term1 = log((X1*Del(i,1))+(X2*Del(i,2))+(X3*Del(i,3)));
        term2 = (X1*Del(1,i))/(X1+(X2*Del(1,2))+(X3*Del(1,3)));
        term3 = (X2*Del(2,i))/((X1*Del(2,1))+X2+(X3*Del(2,3)));
        term4 = (X3*Del(3,i))/((X1*Del(3,1))+(X2*Del(3,2))+X3);
        lng = 1 - term1 - term2 - term3 - term4;
        g(i) = exp(lng);
    end
end