


scatter([-2.5, 0, 2], [1.8E11, 2.9E11, 3.7E11], 100, "r", "filled")
hold on
xlabel("Back Gate Voltage (V)")
ylabel("n (cm^-2)")
scatter([-2.5, 0, 2], [2.3E11, 2.9E11, 4.3E11], 75, "b", "v", "filled")
xlim([-3, 3])
ylim([1E11, 5E11])
legend(["ETH", "Chris"])
grid on