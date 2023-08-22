function U_k = demodulate(U_set, y_old, y_new)
L = size(U_set, 3);
U_temp = U_set(:,:,1);
Temp = real((y_new') * U_set(:,:,1) * y_old);
for i = 2 : L
   Temp_new = real((y_new')*U_set(:,:,i)*y_old);
   if Temp_new > Temp
       U_temp = U_set(:,:,i);
       Temp = Temp_new;
   end
end
U_k = U_temp;
end