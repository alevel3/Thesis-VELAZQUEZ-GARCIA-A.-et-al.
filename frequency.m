function [range_final,freq,data_mean,data_std,ind_test1] =frequency(range,data);


for i=1:max(size(range))-1
    ind_test1{i} = find(data >range(i) & data <range(i+1));
    ind_test(i) = max(size(find(data >range(i) & data <range(i+1))));
end

for i=1:max(size(range))-1
data_mean(i) = nanmean(data(ind_test1{i}));
data_std(i) = nanstd(data(ind_test1{i}));
end

freq = ind_test./sum(ind_test);
range_final = range(1:end-1);
end
