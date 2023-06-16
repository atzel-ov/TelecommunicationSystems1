function [Xi_n, Xq_n]  = bits_to_PSK_16(b)
    
    L = length(b);

    X = zeros(L/4,1);
    k = 1;

    for i = 1 : 4 : L-1 
        
        % b = 0000
        if(b(i)==0 && b(i+1)==0 && b(i+2)==0 && b(i+3)==0)
            X(k) = 1;
            k = k + 1;
        end

        % b = 0001
        if(b(i)==0 && b(i+1)==0 && b(i+2)==0 && b(i+3)==1)
            X(k) = exp(1*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 0011
        if(b(i)==0 && b(i+1)==0 && b(i+2)==1 && b(i+3)==1)
            X(k) = exp(2*(2*pi)*j/16);
            k = k + 1;
        end
        
        % b = 0010
        if(b(i)==0 && b(i+1)==0 && b(i+2)==1 && b(i+3)==0)
            X(k) = exp(3*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 0110
        if(b(i)==0 && b(i+1)==1 && b(i+2)==1 && b(i+3)==0)
            X(k) = exp(4*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 0111
        if(b(i)==0 && b(i+1)==1 && b(i+2)==1 && b(i+3)==1)
            X(k) = exp(5*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 0101
        if(b(i)==0 && b(i+1)==1 && b(i+2)==0 && b(i+3)==1)
            X(k) = exp(6*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 0100
        if(b(i)==0 && b(i+1)==1 && b(i+2)==0 && b(i+3)==0)
            X(k) = exp(7*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1100
        if(b(i)==1 && b(i+1)==1 && b(i+2)==0 && b(i+3)==0)
            X(k) = exp(8*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1101
        if(b(i)==1 && b(i+1)==1 && b(i+2)==0 && b(i+3)==1)
            X(k) = exp(9*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1111
        if(b(i)==1 && b(i+1)==1 && b(i+2)==1 && b(i+3)==1)
            X(k) = exp(10*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1110
        if(b(i)==1 && b(i+1)==1 && b(i+2)==1 && b(i+3)==0)
            X(k) = exp(11*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1010
        if(b(i)==1 && b(i+1)==0 && b(i+2)==1 && b(i+3)==0)
            X(k) = exp(12*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1011
        if(b(i)==1 && b(i+1)==0 && b(i+2)==1 && b(i+3)==1)
            X(k) = exp(13*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1001
        if(b(i)==1 && b(i+1)==0 && b(i+2)==0 && b(i+3)==1)
            X(k) = exp(14*(2*pi)*j/16);
            k = k + 1;
        end

        % b = 1000
        if(b(i)==1 && b(i+1)==0 && b(i+2)==0 && b(i+3)==0)
            X(k) = exp(15*(2*pi)*j/16);
            k = k + 1;
        end

        Xi_n = real(X);
        Xq_n = imag(X);

    end
end