function U_set = get_U_set(M)
i = 1 : M;
Theta = 2 * pi .* (i - 1) ./ M;
Symbol = sqrt(1/2) * exp(1i .* Theta);
X = zeros(2, 2, M^2);
for i = 1 : M
    for j = 1 : M
        index = M * (i - 1) + j;
        X(:, :, index) = [Symbol(i), Symbol(j); -conj(Symbol(j)), conj(Symbol(i))];
    end
end

p = (M^2+M^4)/2;
U = zeros(2, 2, p);
for i = 1 : M^2
    for j = i : M^2
        index = 1/2*(2*M^2-i)*(i-1)+j;
        X_k1 = X(:, :, i);
        X_k2 = X(:, :, j);
        U(:, :, index) = X_k1 * conj(X_k2.') ;
    end
end 

U_set = U;

end