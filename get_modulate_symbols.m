function X = get_modulate_symbols(k, Data)
X1 = Data(k,:);
X2 = [-conj(X1(2)), conj(X1(1))];
X = [X1; X2];
end

