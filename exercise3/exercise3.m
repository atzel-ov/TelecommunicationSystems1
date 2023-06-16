clear
close all
clc

%% %%%%%%%%%%%%%%
%---- Part 1----%
%%%%%%%%%%%%%%%%%

% 1
N = 100;

bit_seq = (sign(randn(4*N, 1)) + 1)/2;

% 2
[Xi_n, Xq_n]  = bits_to_PSK_16(bit_seq);

T = .01;
over = 10;
Ts = T/over;
Fs = 1/Ts;
A = 2;
a = 0.5;

[phi, t] = srrc_pulse(T, over, A, a);

% 3
t_delta = 0:Ts:N*T-Ts;

Xi_delta = (1/Ts)*upsample(Xi_n, over);
Xq_delta = (1/Ts)*upsample(Xq_n, over);

tconv = t(1)+t_delta(1):Ts:t(end)+t_delta(end);

X_i = Ts*conv(Xi_delta, phi);
X_q = Ts*conv(Xq_delta, phi);

figure
plot(tconv, X_i, 'b')
xlabel('Time(s)')
ylabel('X_I')
grid on

figure
plot(tconv, X_q, 'b')
xlabel('Time(s)')
ylabel('X_Q')
grid on


T_total = (N+2*A)*T;

Nf = 2048;
f = linspace(-Fs/2,Fs/2-Fs/Nf, Nf);

X_i_en = fftshift(abs(Ts*fft(X_i,Nf)).^2);
X_q_en = fftshift(abs(Ts*fft(X_q,Nf)).^2);

P_Xi = X_i_en/T_total;
P_Xq = X_q_en/T_total;

figure
plot(f,P_Xi,'b');
xlabel('Frequency(Hz)')
ylabel('P_X_I')
grid on

figure
plot(f,P_Xq,'b');
xlabel('Frequency(Hz)')
ylabel('P_X_Q')
grid on


% 4
F0 = 200;

cos_carrier = (2*cos(2*pi*F0*tconv))';
sin_carrier = (-2*sin(2*pi*F0*tconv))';

X_i_c = X_i.*cos_carrier;
X_q_s = X_q.*sin_carrier;

figure
plot(tconv, X_i_c, 'r')
xlabel('Time(s)')
ylabel('X_I * 2cos(2πFot)')
grid on

figure
plot(tconv, X_q_s, 'r')
xlabel('Time(s)')
ylabel('X_Q *(-2sin(2πFot))')
grid on


X_i_c_en = fftshift(abs(Ts*fft(X_i_c,Nf)).^2);
X_q_s_en = fftshift(abs(Ts*fft(X_q_s,Nf)).^2);

P_Xi_c = X_i_c_en/T_total;
P_Xq_s = X_q_s_en/T_total;

figure
plot(f,P_Xi_c,'r');
xlabel('Frequency(Hz)')
ylabel('P_X_I')
grid on

figure
plot(f,P_Xq_s,'r');
xlabel('Frequency(Hz)')
ylabel('P_X_Q')
grid on


% 5
X_t = X_i_c + X_q_s;

figure
plot(tconv, X_t, 'm')
xlabel('Time(s)')
ylabel('X(t)')
grid on



X_en = fftshift(abs(Ts*fft(X_t,Nf)).^2);

P_x = X_en/T_total;

figure
plot(f,P_x,'m');
xlabel('Frequency(Hz)')
ylabel('P_X')
grid on

% 7
SNRdb = 20;

var_w = (1/Ts)*10^(-SNRdb/10);

W = sqrt(var_w)*randn(1, length(tconv));

Y = X_t + W';

%% 8
Y_i = Y.*cos_carrier/2;
Y_q = Y.*sin_carrier/2;

Y_i_en = fftshift(abs(Ts*fft(Y_i,Nf)).^2);
Y_q_en = fftshift(abs(Ts*fft(Y_q,Nf)).^2);

P_Yi = Y_i_en/T_total;
P_Yq = Y_q_en/T_total;

figure
plot(tconv, Y_i, 'r')
xlabel('Time(s)')
ylabel('Y_I')
grid on

figure
plot(tconv, Y_q, 'r')
xlabel('Time(s)')
ylabel('Y_Q')
grid on

figure
plot(f,P_Yi,'r');
xlabel('Frequency(Hz)')
ylabel('P_Y_I')
grid on

figure
plot(f,P_Yq,'r');
xlabel('Frequency(Hz)')
ylabel('P_Y_Q')
grid on


% 9
t_out = tconv(1)+t(1):Ts:tconv(end)+t(end);

T_total = (N+4*A)*T;

Y_i_delta = Ts*conv(Y_i, phi);
Y_q_delta = Ts*conv(Y_q, phi);

P_Yi_delta = fftshift(abs(Ts*fft(Y_i_delta,Nf)).^2)/T_total;
P_Yq_delta = fftshift(abs(Ts*fft(Y_q_delta,Nf)).^2)/T_total;

figure
plot(t_out, Y_i_delta, 'b')
xlabel('Time(s)')
ylabel('Y_I_,_δ')
grid on

figure
plot(t_out, Y_q_delta, 'b')
xlabel('Time(s)')
ylabel('Y_Q_,_δ')
grid on

figure
plot(f,P_Yi_delta,'b');
xlabel('Frequency(Hz)')
ylabel('P_Y_I_,_δ')
grid on

figure
plot(f,P_Yq_delta,'b');
xlabel('Frequency(Hz)')
ylabel('P_Y_Q_,_δ')
grid on

% 10

Y_i_k = zeros(100,1);
Y_q_k = zeros(100,1);

k = 1;

for i = 2*A*over+1:over:length(t_out)-2*A*over
    Y_i_k(k) = Y_i_delta(i);
    Y_q_k(k) = Y_q_delta(i);
    k = k+1;
end

Y = [Y_i_k,Y_q_k];


scatterplot(Y)
grid on


% 11
[est_X, est_bit_seq] = detect_PSK_16(Y);


% 12
X = [Xi_n, Xq_n];
num_of_symbol_errors = symbol_errors(est_X, X);
disp(num_of_symbol_errors)


% 13
num_of_bit_errors = bit_errors(est_bit_seq, bit_seq);
disp(num_of_bit_errors)



%% %%%%%%%%%%%%%%
%---- Part 2----%
%%%%%%%%%%%%%%%%%

p_error_symbol = zeros(14,1);
p_error_s = zeros(14,1);
p_error_bit = zeros(14,1);
p_error_b = zeros(14,1);

j = 1;

for SNRdb = -2:2:24

    for K = 1:1000
        var_w = (1/Ts)*10^(-SNRdb/10);
        
        W = sqrt(var_w)*randn(1, length(tconv));
        
        Y = X_t + W';
        
        Y_i = Y.*cos_carrier/2;
        Y_q = Y.*sin_carrier/2;
       
        t_out = tconv(1)+t(1):Ts:tconv(end)+t(end);
        
        T_total = (N+4*A)*T;
        
        Y_i_delta = Ts*conv(Y_i, phi);
        Y_q_delta = Ts*conv(Y_q, phi);
        
        
        Y_i_k = zeros(100,1);
        Y_q_k = zeros(100,1);
        
        k = 1;
        
        for i = 2*A*over+1:over:length(t_out)-2*A*over
            Y_i_k(k) = Y_i_delta(i);
            Y_q_k(k) = Y_q_delta(i);
            k = k+1;
        end
        
        Y = [Y_i_k,Y_q_k];
        
        
        [est_X, est_bit_seq] = detect_PSK_16(Y);
        
        
        X = [Xi_n, Xq_n];
        num_of_symbol_errors(K) = symbol_errors(est_X, X);
        
        num_of_bit_errors(K) = bit_errors(est_bit_seq, bit_seq);
    end
    
    % Experimental error
    p_error_symbol(j) = sum(num_of_symbol_errors)/(K*length(X));

    SNR = 10^(SNRdb/10);
    % Theoretical error
    p_error_s(j) = 2*Q(sqrt(2*(SNR))*sin(pi/16));
    
    % Experimental error
    p_error_bit(j) = sum(num_of_bit_errors)/(K*length(bit_seq));

    % Theoretical error
    p_error_b(j) = p_error_s(j)/4;

    j = j+1;

end

SNRdb = -2:2:24;
figure
semilogy(SNRdb, p_error_symbol, 'b')    
hold on
semilogy(SNRdb, p_error_s, 'g') 
xlabel('SNR_d_b')
ylabel('P^H^A^T(E_s_y_m_b_o_l)')
grid on

figure
semilogy(SNRdb,p_error_bit, 'b')
hold on
semilogy(SNRdb,p_error_b, 'g')
xlabel('SNR_d_b')
ylabel('P^H^A^T(E_b_i_t)')
grid on