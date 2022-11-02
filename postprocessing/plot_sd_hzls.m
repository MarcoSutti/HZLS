function [] = plot_sd_hzls( info_HZLS, Fh_ref )

% function [] = plot_sd_hzls( info_HZLS, Fh_ref )
% Purpose: Plot the results given by ...

% Created:     2021.02.17
% Last change: 2021.11.02

%   Jan 17, 2021:
%       Created.

options_plot;

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
diff_f_rel_HZLS(diff_f_rel_HZLS==0) = eps;

figure(1);

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 1;

% Reference lines:
semilogy( sqrt(eps)*ones(length(diff_f_rel_HZLS), 1 ), '--', 'Color', [ 0.75 0.75 0.75 ], 'LineWidth', 2 );
hold on
grid on
semilogy( eps*ones(length(diff_f_rel_HZLS),1), '--', 'Color', [ 0.75 0.75 0.75 ], 'LineWidth', 2 );

% Cost function and gradient norm:
% For Hager-Zhang line search (HZLS)
handle_array(1) = semilogy( diff_f_rel_HZLS, '-^', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    red, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(diff_f_rel_HZLS) );

handle_array(2) = semilogy( normalized_gradnorm_HZLS, '--^', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    red, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(normalized_gradnorm_HZLS) );

handle_array(3) = semilogy( opt_variable_array_HZLS, ':^', 'Color', red, ...
    'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...,
    red, 'MarkerSize', 7, 'MarkerIndices', 1:stride:length(opt_variable_array_HZLS) );

drawnow;
for i=1:3
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend( handle_array, ...
    {'$|f_{k} - f_{\ast}|/|f_{\ast}|$, HZ', ...
    '$\|g_{k}\|/\|g_{0}\|$, HZ', ...
    '$\|X_{k} - X_{\ast}\|/\|X_{\ast}\|$, HZ'}, ...
    'FontSize', 11, 'Location', 'NE' );

drawnow;
for i=1:3
    lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end

xlabel('iteration $k$ of steepest descent')
xticks( 1:stride:length(diff_f_rel_HZLS) );
xticklabels( 0:stride:length(diff_f_rel_HZLS) );
% xlim( [ 0, length(diff_f_rel_HZLS) ] )
ylim( [ 5e-17, 100 ] )

% Make the window of the figure looking a bit nicer
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 1 0.5 0.66];

end