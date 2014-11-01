exectuable_path = fullfile('build','simsimhash');

datadir = '/home/njclimer/data/p4';

datasets = {'enron','as-733','enron_noempty','p2p-Gnutella','reality_mining_voices'};
%close all;
fclose all;
for i=1:numel(datasets)
    system([exectuable_path ' ' fullfile(datadir,datasets{i}) ' '...
        datasets{i} '.coeff ' datasets{i} '.anom']);
    
    coeffs = dlmread([datasets{i} '.coeff']);
    cutoff = coeffs(1);
    coeffs(1)=[];
    mr = mean(abs(diff(coeffs)));
    md = median(coeffs);
    disp(cutoff-md+3*mr);
    cutoff = md-3*mr;
    figure;
    scatter(find(coeffs>cutoff),coeffs(coeffs>cutoff),'blue')
    hold on;
    scatter(find(coeffs<cutoff),coeffs(coeffs<cutoff),'red','x');
    plot([0,numel(coeffs)],[cutoff,cutoff],'--');
    if (any(coeffs<cutoff))
        legend('passing','anomaly','cutoff');
    else
        legend('passing','cutoff');
    end
    title(['Simsimhash for ',datasets{i}]);
    xlabel('Iteration (0 indexed)');
    ylabel('Simularity Score');    
end

