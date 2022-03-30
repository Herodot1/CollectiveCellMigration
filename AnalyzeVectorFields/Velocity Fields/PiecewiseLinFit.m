N = 50*12;
x = 1:N;
y = randn(1, N);
p = cell(1, N-60);


In = 19; 
x = log10(VelFieldData(In).tau);
y = log10(nanmean(VelFieldData(In).MSDAll(1:length(VelFieldData(In).tau),:),2));
NumSamples = 10;

p = [];
for i = 1:400%length(x)-NumSamples
p(i,:) = polyfit(x(i:i+NumSamples-1), y(i:i+NumSamples-1), 1);
end
mean(p(1:400,1))
figure;
plot(p(:,1))
