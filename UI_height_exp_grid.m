%% check height
N_layer = 20;
h = 0;
dh = 1;
alpha = 0.1;
h_line = zeros(N_layer, 1);
for num = 0 : N_layer
    h_line (num +1) = h;
    h = h + ((1 + alpha)^num)*dh;
end
disp(h);
disp(((1+alpha)^(N_layer+1)-1)/(alpha));
h_line2=zeros(N_layer, 1);
for num = 0 : N_layer-1
    h_line2(num + 2) = ((1+alpha)^(num+1)-1)/(alpha);
end
plot(h_line);
hold on;
plot(h_line2);
hold on;