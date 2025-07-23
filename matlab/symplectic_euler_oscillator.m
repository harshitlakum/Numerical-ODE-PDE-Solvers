% symplectic_euler_oscillator.m
% Semi-implicit (symplectic) Euler for y''(t) + y(t) = sin(t)

T = 10;
h = 0.01;
N = round(T/h) + 1;
t = linspace(0, T, N);

x = zeros(1, N);   % position
v = zeros(1, N);   % velocity

x(1) = 0;
v(1) = 1;

for n = 2:N
    x(n) = x(n-1) + h * v(n-1);
    v(n) = v(n-1) + h * (-x(n-1) + sin(t(n-1)));
end

figure;
subplot(1,2,1)
plot(t, x, 'b')
xlabel('Time t')
ylabel('Position x')
title('Position vs Time')
grid on

subplot(1,2,2)
plot(t, v, 'r')
xlabel('Time t')
ylabel('Velocity v')
title('Velocity vs Time')
grid on