function H = get_channels(Gain_g0, Array_reflection, Mu_g1, Sigma_g1, Sigma_h1)
% IRS reflection link
g1 = mvnrnd(Mu_g1, Sigma_g1, 1);
g1 = g1(:);
g1_bar = Gain_g0' * Array_reflection * g1;

% Direct link between AP and User
h1 = sqrt(Sigma_h1) * randn;

H = [h1; g1_bar];
end

