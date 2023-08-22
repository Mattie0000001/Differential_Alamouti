function Y = get_received_symbols(N, P, H, S)
j = sqrt(-1);
n1 = sqrt(N/2) * randn + sqrt(N/2) * randn * j;
n2 = sqrt(N/2) * randn + sqrt(N/2) * randn * j;
Noise = [n1; n2];
Y = sqrt(P) * S * H + Noise;

end