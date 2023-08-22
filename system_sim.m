clear all
clc

%% basic system constant
M = 8; % M-ary Modulation

Lamda = 0.05; % wavelength

Location_AP = 0.5;
IRS = [0 0 0];
AP = [Location_AP 0 0];
User = [30, 50, 0];

Dis_h1 = norm(AP - User); % distance between AP and User
Dis_g1 = norm(IRS - User); % distance between AP and User

N_dbm = -85;  % noise, in dBm
N = db2pow(N_dbm); % noise, in mW

Alpha = 3; % path loss component
Beta_dB = -30; % reference path gain, in dB
Beta = db2pow(Beta_dB); % reference path gain, in unit 1

% P_dbm = 15 : 1 : 35; % transmit power, in dBm
P_dbm = 30; % transmit power, in dBm
P = db2pow(P_dbm); % transmit power, in mW

Gamma = pow2db(P / N); % transmit SNR

Path_loss_g1 = Beta * Dis_g1 ^ (- Alpha);
Path_loss_h1 = Beta * Dis_h1 ^ (- Alpha);

L_pack = 100; % Packedge length, even number
Right_cnt = 0;
Error_cnt = 0;
Num_pack = 10; % Numbers of packedge
L_block = L_pack / 2;

%% IRS coefficient constant
Num_y = 105; 
Num_z = Num_y; 
Num_IRS = Num_y * Num_z; % number of elements
Delta = Lamda / 2; % the width of the elements
Epsilon = Delta / Location_AP;
Area_single = Delta ^ 2; % area of single elements
xi = 1; % array occupation ratio

IRS_UPA_y = repmat(-(Num_y - 1) / 2 * Delta : Delta : (Num_y -  1) / 2 * Delta, [Num_z, 1]);
IRS_UPA_y = IRS_UPA_y(:);
IRS_UPA_z = repmat((-(Num_z - 1) / 2 * Delta : Delta : (Num_z -  1) / 2 * Delta)', [1, Num_y]);
IRS_UPA_z = IRS_UPA_z(:);
IRS_UPA = [zeros(Num_IRS, 1), IRS_UPA_y, IRS_UPA_z];
IRS_UPA = IRS + IRS_UPA;

Vec_AP_IRS = AP - IRS_UPA;
Dis_g0 = sqrt(sum(abs(Vec_AP_IRS) .^ 2, 2)); % distance between AP and User, is an array
Path_loss_g0 = Area_single * Location_AP ./ (4 * pi .* Dis_g0 .^ 3);

IRS_vec_reflection = ones(Num_IRS, 1); % the reflection vector of IRS

%% channel model constant
% IRS reflection link
Gain_g0 = sqrt(Area_single * Location_AP ./ (4 * pi .* Dis_g0 .^ 3)) .* exp(- 1i * 2 * pi .* Dis_g0 ./ Lamda);
Array_reflection = diag(IRS_vec_reflection);

Mu_g1 = zeros(1, Num_IRS);
Sigma_g1 = 2 * Path_loss_g1 * eye(Num_IRS);

% Direct link between AP and User
Sigma_h1 = Path_loss_h1;

%% variables initialization
y_old = [0, 0]; % y_k
y_new = [0, 0]; % y_k+1
S_old = [0, 0; 0, 0]; % X_k
S_new = [1, 1; -1, 1]*sqrt(1/2); % X_k+1
H = [0; 0]; 

%% get the set of U
U_set = get_U_set(M);

for np = 1 : Num_pack
    %% Data initial
    [Data, Index] = initial_data(M, L_block);    
    %% channel gain
    H = get_channels(Gain_g0, Array_reflection, Mu_g1, Sigma_g1, Sigma_h1);
    % first transmission
    y_new = get_received_symbols(N, P, H, S_new);
    %% transmission
    for k = 1 : L_block
      S_old = S_new;
      y_old = y_new;
      %% symbol modulation
      S_new = get_modulate_symbols(k, Data);
      U_theo = S_new * S_old';
      %% new transmission
      y_new = get_received_symbols(N, P, H, S_new);
      %% demodulation
            U_k = demodulate(U_set, y_old, y_new);
      %% count the error
      if isequal(round(U_theo*10^6)/10^6, round(U_k*10^6)/10^6)
        Right_cnt = Right_cnt + 1;
      else
        Error_cnt = Error_cnt + 1;
      end
    end
end

SER = 2 * Error_cnt / (Num_pack * L_pack);
  



