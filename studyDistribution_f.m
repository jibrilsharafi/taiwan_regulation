function studyDistribution_f(f)
    % To plot the residuals, we reduce the number of points by a reduction
    % factor
    reduction_factor = 1e3;
    f_reduced = f(1:reduction_factor:end);
    % To plot an historam of the distribution of the frequency, we choose the
    % number of bins according to Sturge’s Rule
    n_bins = floor(1 + 3.322*log10(length(f)));
    
    % Create a subplot
    figure()
    subplot(1,2,1)
    hold on
    grid on
    xlabel('Frequency [Hz]')
    ylabel('Count [-]')
    legend(Location="northwest")
    histogram(f,NumBins=n_bins,DisplayName=sprintf('Mean = %2.4f, stf = %0.4f', [mean(f), std(f)]))
    
    subplot(1,2,2)
    hold on
    grid on
    xlabel('Time [s]')
    ylabel('Residuals [Hz]')
    legend
    plot(1:reduction_factor:length(f), f_reduced-mean(f_reduced), 'o', DisplayName='Residuals')
    plot([0 length(f)], [0 0], 'r--', LineWidth=1, DisplayName='Mean value')
end