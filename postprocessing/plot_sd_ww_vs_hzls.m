function [] = plot_sd_ww_vs_hzls( info_WWLS, info_HZLS, Fh_ref, pars )

% function [] = plot_sd_ww_vs_hzls( info_NMLS_w_noise, info_WWLS, info_HZLS, pars )
% Purpose: Plot the results given by 'Driver_SD_WW_vs_HZLS'.

% Created:     2020.05.23
% Last change: 2022.11.02

%   Jan 26, 2021:
%       Added plot for the line search with safeguarded quadratic
%       interpolation.
%   Aug 23, 2020:
%       Added if statement for for diff_f_rel and for the legend.
%   Aug 19, 2020:
%       Added Rosenbrock's function.

options_plot;

cost_array_WWLS = [info_WWLS.cost];
gradnorm_array_WWLS = [info_WWLS.gradnorm];
opt_variable_array_WWLS = [info_WWLS.opt_variable];
normalized_gradnorm_WWLS = gradnorm_array_WWLS/gradnorm_array_WWLS(1);
diff_f_rel_WWLS = abs(cost_array_WWLS - Fh_ref)/abs(Fh_ref);

cost_array_HZLS = [info_HZLS.cost];
gradnorm_array_HZLS = [info_HZLS.gradnorm];
opt_variable_array_HZLS = [info_HZLS.opt_variable];
normalized_gradnorm_HZLS = gradnorm_array_HZLS/gradnorm_array_HZLS(1);
diff_f_rel_HZLS = abs(cost_array_HZLS - Fh_ref)/abs(Fh_ref);

%--------------------------------------------------------------------------
% Convergence plot
%--------------------------------------------------------------------------
% Set to epsilon machine the values that are exactly zero in order to
% avoid "holes" in the plot
% diff_f_rel_NMLS_w_noise(diff_f_rel_NMLS_w_noise==0) = eps;
diff_f_rel_WWLS(diff_f_rel_WWLS==0) = eps;
diff_f_rel_HZLS(diff_f_rel_HZLS==0) = eps;

figure(1);

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 20;

% Reference lines:
semilogy( sqrt(eps)*ones(length(diff_f_rel_WWLS), 1 ), '--', 'Color', [ 0.75 0.75 0.75 ], 'LineWidth', 2 );
hold on
grid on
semilogy( eps*ones(length(diff_f_rel_WWLS),1), '--', 'Color', [ 0.75 0.75 0.75 ], 'LineWidth', 2 );

% Cost function and gradient norm:
% For standard line search (weak Wolfe, WWLS)
handle_array(1) = semilogy( diff_f_rel_WWLS, '-o', 'Color', green, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    green, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(diff_f_rel_WWLS) );

handle_array(2) = semilogy( normalized_gradnorm_WWLS, '--o', 'Color', green, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    green, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(normalized_gradnorm_WWLS) );

handle_array(3) = semilogy( opt_variable_array_WWLS, ':o', 'Color', green, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    green, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(opt_variable_array_WWLS) );

% For Hager-Zhang line search (HZLS)
handle_array(4) = semilogy( diff_f_rel_HZLS, '-^', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    red, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(diff_f_rel_HZLS) );

handle_array(5) = semilogy( normalized_gradnorm_HZLS, '--^', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    red, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(normalized_gradnorm_HZLS) );

handle_array(6) = semilogy( opt_variable_array_HZLS, ':^', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    red, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(opt_variable_array_HZLS) );

drawnow;
for i=1:6
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend( handle_array, ...
    {'$|f_{k} - f_{\ast}|/|f_{\ast}|$, WW', ...
    '$\|g_{k}\|/\|g_{0}\|$, WW', ...
    '$\|X_{k} - X_{\ast}\|/\|X_{\ast}\|$, WW', ...
    '$|f_{k} - f_{\ast}|/|f_{\ast}|$, HZ', ...
    '$\|g_{k}\|/\|g_{0}\|$, HZ', ...
    '$\|X_{k} - X_{\ast}\|/\|X_{\ast}\|$, HZ'}, ...
    'FontSize', 11, 'Location', 'NE' );

drawnow;
for i=1:6
    lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end

xlabel('iteration $k$ of steepest descent')
xticks( 1:stride:length(diff_f_rel_WWLS) );
xticklabels( 0:stride:length(diff_f_rel_WWLS) );
xlim( [ 0, length(diff_f_rel_WWLS) ] )
ylim( [ 5e-17, 10 ] )

% Make the window of the figure looking a bit nicer
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 1 0.5 0.66];

%--------------------------------------------------------------------------
% Save the plot to file
if strcmp(pars.fgname,'easy_quadratic')
    fileName = [ 'plots/', pars.fgname, '_', pars.var_type, '_hz_vs_ww' ];
%     export_fig Plots/easy_quadratic_matrix_hz_vs_ww.pdf -pdf -cmyk -transparent;
elseif strcmp(pars.fgname,'rosenbrock_function')
    fileName = [ 'plots/', pars.fgname, '_hz_vs_ww' ];
%     export_fig Plots/rosenbrock_function_hz_vs_ww.pdf -pdf -cmyk -transparent;
elseif strcmp(pars.fgname,'goldstein_price')
%     export_fig Plots/goldstein_price_hz_vs_ww.pdf -pdf -cmyk -transparent;
end
pause(0.5)
saveas( gcf, fileName, 'epsc' );
fprintf('+-----------------------------------------------------------------+\n');
fprintf('Saved graph to file %s.eps.\n', fileName);
%--------------------------------------------------------------------------

end