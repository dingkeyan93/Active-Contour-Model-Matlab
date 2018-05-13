function printStats(ops, total_stats)
disp('  ================ time statistics (MEX) ==================');
for i=1:length(ops),
    disp(sprintf('%26s\t %.6g s', ops{i}.name, total_stats(i)));
end