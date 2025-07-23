% euler_logistic_compare.m
% Explicit Euler for logistic growth IVP, compared to the exact solution

T = 10;
y0 = 1;
hs = [0.2, 0.1, 0.05, 0.01];
colors = lines(length(hs));
figure; hold on;

for k = 1:length(hs)
    h = hs(k);
    t = 0:h:T;
    y = zeros(1,length(t));
    y(1) = y0;
    for n = 2:length(t)
        y(n) = y(n-1) + h*(1 - y(n-1)/100)*y(n-1);
    end

    plot(t, y, '--', 'Color', colors(k,:), 'DisplayName', sprintf('Euler (h=%.2f)',h));

    if k == 1 % Plot exact solution only once
        t_dense = linspace(0, T, 1000);
        y_exact = 100 ./ (1 + 99*exp(-t_dense));
        plot(t_dense, y_exact, 'k-', 'LineWidth',1.5, 'DisplayName','Exact Solution');
    end
end

xlabel('t')
ylabel('y(t)')
title('Logistic Growth: Explicit Euler vs. Exact')
legend show
grid on
hold off