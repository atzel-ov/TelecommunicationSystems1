function X = bits_to_2PAM(b)

N = 100;
X = zeros(N,1);

for i=1:(N)

    if(b(i) == 0)
        X(i) = 1;
    else
        X(i) = -1;
    end

end
