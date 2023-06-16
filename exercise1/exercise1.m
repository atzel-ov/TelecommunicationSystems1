clear;
close all;
clc;



%% %%%%% A1 %%%%%%%%

T = .001;
over = 10;
A = 4;

Ts = T/over;

figure
for a = 0:0.5:1
    [phi,t] = srrc_pulse(T, over, A, a);
    plot(t,phi)
    hold on
end
title("SRRC functions for different roll-offs")
legend('a = 0','a = 0.5','a = 1')
grid on


%%%%%%% A2 %%%%%%%%

%a

Fs = 1/Ts;
Nf = 2048;
f = -Fs/2:Fs/Nf:Fs/2-Fs/Nf;

figure
for a = 0:0.5:1
    [phi,t] = srrc_pulse(T, over, A, a);
    Phi = fftshift(abs(Ts.*fft(phi,Nf)).^2);
    plot(f,Phi)
    hold on
end
title("SRRC energy spectral density for different roll-offs")
legend('a = 0','a = 0.5','a = 1')
grid on


%b  -Same with (a) for semilogy instead of plot

figure
for a = 0:0.5:1
    [phi,t] = srrc_pulse(T, over, A, a);
    Phi = fftshift(abs(Ts.*fft(phi,Nf)).^2);
    semilogy(f,Phi)
    hold on
end
title("SRRC energy spectral density in a logarithmic scale")
legend('a = 0','a = 0.5','a = 1')
grid on


%%%%%%% A3 %%%%%%%%

%a

for a = 0:0.5:1
    BW = 0.5*(1+a)/T
end

figure
for a = 0:0.5:1
    [phi,t] = srrc_pulse(T, over, A, a);
    Phi = fftshift(abs(Ts.*fft(phi,Nf)).^2);
    semilogy(f,Phi)
    hold on
end
c = (T/(10^5)).*ones(length(f));
semilogy(f,c)

title("SRRC energy spectral density in a logarithmic scale")
legend('a = 0','a = 0.5','a = 1')
grid on



%% %%%%% B1 %%%%%%%%


%1

T = .001;
over = 10;
A = 4;

Ts = T/over;

figure
for a = 0:0.5:1

    subplot(3,1,2*a+1)
   
    for k = 0:1:2
        
        
        
        [phi,t] = srrc_pulse(T, over, A, a);
        delayed = [zeros(1,k*over) phi(1:(end-k*over))];
        % k*over->end we lose some information on the delayed phi
        % That doesnt impact the multiplication with phi since phi is 0
        % after k*over
        hold on
        plot(t, delayed)
    end
    hold on
    title("SRRC functions for different roll-offs and k")
    legend('k = 0','k = 1','k = 2')
    grid on
end

figure
for a = 0:0.5:1

    subplot(3,1,2*a+1)

    for k = 0:1:2
        
        [phi,t] = srrc_pulse(T, over, A, a);
        delayed = [zeros(1,k*over) phi(1:(end-k*over))];
        fpf = phi.*delayed;     
        hold on
        plot(t,fpf)
        legend('k = 0','k = 1','k = 2')
        grid on

        
    end
end

figure
for a = 0:0.5:1
    subplot(3,1,2*a+1)
    for k = 0:1:3
        [phi,t] = srrc_pulse(T, over, A, a);
        delayed = [zeros(1,k*over) phi(1:(end-k*over))];
        fpf = phi.*delayed; % product phi*phi(t-kT)
        rff = Ts*sum(fpf).*ones(size(t));
        plot(t,rff)
        legend('k = 0','k = 1','k = 2','k = 4')
        grid on
        hold on

    end
    axis([-4*T 4*T -0.2 1.2])
end



%% %%%%% C %%%%%%%%

T = .001;
over = 10;
A = 4;
a = 0.5;

Ts = T/over;


N = 100;

b = (sign(randn(N, 1)) + 1)/2;

X = bits_to_2PAM(b);

t_delta = [0:Ts:N*T - Ts];

X_delta = 1/Ts*upsample(X, over);

figure
stem(t_delta, X_delta)
title('X_d_e_l_t_a')



[phi, t] = srrc_pulse(T, over, A, a);
    
tconv = [min(t)+min(t_delta):Ts:max(t)+max(t_delta)];

x_conv = Ts*conv(phi,X_delta);

figure
plot(tconv, x_conv)
grid on
title('Transmitter output')






[phi, t] = srrc_pulse(T, over, A, a);
trec = [min(tconv)+min(t):Ts:max(tconv)+max(t)];
z = Ts*conv(x_conv,phi);

figure
plot(trec,z)
hold on
stem([0 : N-1]*T, X);
grid on 
title('Receiver')






