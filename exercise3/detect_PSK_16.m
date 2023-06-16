function [est_X, est_bit_seq] = detect_PSK_16(Y)

    L = length(Y);
    M = 16;
    neighbor = zeros(M,1);
    est_X = zeros(L,2);
    est_bit_seq = zeros(4*L,1);

    for i = 1:L
        for m = 1:M
            % Finding every neighbor
            neighbor(m) = sqrt((Y(i,1) - cos((m-1)*2*pi/M)).^2 + (Y(i,2) - sin((m-1)*2*pi/M)).^2);
        end
        
        % Finding the Nearest 
        nearest_neighbor = min(neighbor);
        
        % Estimating 
        for n = 1:M
            if(neighbor(n) == nearest_neighbor)
                est_X(i,1) = cos((n-1)*2*pi/M);
                est_X(i,2) = sin((n-1)*2*pi/M);
            end
        end

    end
    
    k = 1;
    for j = 1:4:4*L-1

        % 0000
        if(est_X(k,1) == cos(0*2*pi/M) && est_X(k,2) == sin(0*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 0;
        end
        
        % 0001
        if(est_X(k,1) == cos(1*2*pi/M) && est_X(k,2) == sin(1*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 1;
        end

        % 0011
        if(est_X(k,1) == cos(2*2*pi/M) && est_X(k,2) == sin(2*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 1;
        end
        
        % 0010
        if(est_X(k,1) == cos(3*2*pi/M) && est_X(k,2) == sin(3*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 0;
        end

        % 0110
        if(est_X(k,1) == cos(4*2*pi/M) && est_X(k,2) == sin(4*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 0;
        end

        % 0111
        if(est_X(k,1) == cos(5*2*pi/M) && est_X(k,2) == sin(5*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 1;
        end

        % 0101
        if(est_X(k,1) == cos(6*2*pi/M) && est_X(k,2) == sin(6*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 1;
        end

        % 0100
        if(est_X(k,1) == cos(7*2*pi/M) && est_X(k,2) == sin(7*2*pi/M))
            est_bit_seq(j) = 0;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 0;
        end

        % 1100
        if(est_X(k,1) == cos(8*2*pi/M) && est_X(k,2) == sin(8*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 0;
        end

        % 1101
        if(est_X(k,1) == cos(9*2*pi/M) && est_X(k,2) == sin(9*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 1;
        end

        % 1111
        if(est_X(k,1) == cos(10*2*pi/M) && est_X(k,2) == sin(10*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 1;
        end

        % 1110
        if(est_X(k,1) == cos(11*2*pi/M) && est_X(k,2) == sin(11*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 1;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 0;
        end
        
        % 1010
        if(est_X(k,1) == cos(12*2*pi/M) && est_X(k,2) == sin(12*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 0;
        end

        % 1011
        if(est_X(k,1) == cos(13*2*pi/M) && est_X(k,2) == sin(13*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 1;
            est_bit_seq(j+3) = 1;
        end

        % 1001
        if(est_X(k,1) == cos(14*2*pi/M) && est_X(k,2) == sin(14*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 1;
        end

        % 1000
        if(est_X(k,1) == cos(15*2*pi/M) && est_X(k,2) == sin(15*2*pi/M))
            est_bit_seq(j) = 1;
            est_bit_seq(j+1) = 0;
            est_bit_seq(j+2) = 0;
            est_bit_seq(j+3) = 0;
        end

        k = k + 1;
    end

end