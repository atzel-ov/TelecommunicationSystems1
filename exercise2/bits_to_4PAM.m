function X = bits_to_4PAM(a,b)

N = 100;
X = zeros(N/2,1);

for i=1:(N/2)
    for j=1:(N/2)
        if(a(i)==0 && b(j) ==0)
            X(i) = 3;
        elseif(a(i)==0&&b(j)==1)
            X(i) = 1;
        elseif(a(i)==1 &&b(j)==1)
            X(i)=-1;
        else
            X(i)=-3;
        end
    end
end