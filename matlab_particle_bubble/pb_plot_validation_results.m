function pb_plot_validation_results(validation)
% PB_PLOT_VALIDATION_RESULTS Plot validation outputs for the adaptive framework.
figure('Color','w','Position',[100 100 1100 800]);

subplot(2,2,1);
loglog(validation.gap*1e6, validation.forceNumerical, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
loglog(validation.gap*1e6, validation.forceTheory, 's--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('gap [\mum]'); ylabel('drainage force [N]');
legend('numerical','Taylor theory','Location','southwest');
title('Hydrodynamic drainage force'); grid on;

subplot(2,2,2);
semilogx(validation.gap*1e6, 100.0*validation.relativeError, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('gap [\mum]'); ylabel('relative error [%]');
title('Validation error'); grid on;

subplot(2,2,3);
semilogx(validation.gap*1e6, validation.nCells, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogx(validation.gap*1e6, validation.minDr*1e6, 's--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('gap [\mum]'); ylabel('cell count / min \Deltar [\mum]');
legend('cell count','min \Deltar','Location','best');
title('Adaptive-grid response'); grid on;

subplot(2,2,4);
if isfield(validation.sampleGrid, 'centers')
    stairs(validation.sampleGrid.faces*1e6, [validation.sampleGrid.pressure(1); validation.sampleGrid.pressure; 0], 'LineWidth', 1.2); hold on;
    yyaxis right;
    plot(validation.sampleGrid.centers*1e6, validation.sampleGrid.indicator, 'k-', 'LineWidth', 1.5);
    ylabel('refinement indicator [-]');
    yyaxis left;
    ylabel('pressure [Pa]');
    xlabel('r [\mum]');
    title(sprintf('Sample adaptive grid at gap = %.2f \\mum', validation.sampleGrid.gap*1e6));
    grid on;
else
    text(0.1, 0.5, 'No sample grid stored.', 'FontSize', 12);
    axis off;
end

sgtitle('Particle-bubble lubrication validation with modular adaptive framework');
end
