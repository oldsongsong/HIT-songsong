function pb_plot_dynamic_results(dynamic, par)
% PB_PLOT_DYNAMIC_RESULTS Plot quasi-static drainage dynamics outputs.
figure('Color','w','Position',[120 120 1100 800]);

subplot(2,2,1);
plot(dynamic.time*1e3, dynamic.gap*1e6, 'b-', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('minimum gap [\mum]');
title('Gap evolution during drainage'); grid on;

subplot(2,2,2);
plot(dynamic.time*1e3, dynamic.urel*1e3, 'r-', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('approach speed [mm/s]');
title('Relative approach speed'); grid on;

subplot(2,2,3);
plot(dynamic.time*1e3, dynamic.Fhydro*1e6, 'k-', 'LineWidth', 1.5); hold on;
plot(dynamic.time*1e3, dynamic.Ftheory*1e6, 'k--', 'LineWidth', 1.2);
plot(dynamic.time*1e3, dynamic.Fdlvo*1e6, 'm-.', 'LineWidth', 1.2);
yline(par.Fdrive*1e6, 'b:', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('force [\muN]');
legend('numerical hydro','Taylor hydro','disjoining','driving force','Location','best');
title('Force balance'); grid on;

subplot(2,2,4);
plot(dynamic.time*1e3, dynamic.nCells, 'g-', 'LineWidth', 1.5); hold on;
plot(dynamic.time*1e3, dynamic.minDr*1e6, 'c--', 'LineWidth', 1.5);
xlabel('time [ms]'); ylabel('cell count / min \Deltar [\mum]');
legend('cell count','min \Deltar','Location','best');
title('Dynamic adaptive-grid response'); grid on;

sgtitle('Adaptive particle-bubble drainage simulation');
end
