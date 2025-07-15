clear;
clc;

num_sim = 1;
tti = 10^-3;

% load .mat file from scripts in 1) DL Throughput Simulation with Perfect 
% ACK NACK subfolder
% ex. load("filename.mat", "SNRIn", "repeat_counter_arr_sim")
% load for HARQ Type-III with IR

SNRIdx = length(SNRIn);

latency_HARQ = [];

for i=1:num_sim
    latency_HARQ_temp = [];

    for j=1:SNRIdx
        latency_arr = [];

        for l=1:length(repeat_counter_arr_sim{i, j})
            repeat_data = repeat_counter_arr_sim{i, j};

            repeat_data(repeat_data > 3) = 3;
            
            latency_temp = 2 * tti + tti + 2 * (repeat_data(1, l)+1) * 3 * tti + 2*(repeat_data(1, l)+1) * tti;
            latency_arr = [latency_arr, latency_temp];
        end
        latency_HARQ_temp = [latency_HARQ_temp, mean(latency_arr)];
    end
    latency_HARQ = [latency_HARQ; latency_HARQ_temp];
end

% load .mat file from scripts in 1) DL Throughput Simulation with Perfect 
% ACK NACK subfolder
% ex. load("filename.mat", "SNRIn", "repeat_counter_arr_sim")
% load for HARQ Type-I without CC

latency_ARQ = [];

for i=1:num_sim
    latency_ARQ_temp = [];

    for j=1:SNRIdx
        latency_arr = [];

        for l=1:length(repeat_counter_arr_sim{i, j})
            repeat_data = repeat_counter_arr_sim{i, j};

            repeat_data(repeat_data > 3) = 3;
            
            latency_temp = 2 * tti + tti + 2 * (repeat_data(1, l)+1) * 3 * tti + 2*(repeat_data(1, l)+1) * tti;
            latency_arr = [latency_arr, latency_temp];
        end
        latency_ARQ_temp = [latency_ARQ_temp, mean(latency_arr)];
    end
    latency_ARQ = [latency_ARQ; latency_ARQ_temp];
end

% load .mat file from scripts in 1) DL Throughput Simulation with Perfect 
% ACK NACK subfolder
% ex. load("filename.mat", "SNRIn", "repeat_counter_arr_sim")
% load for HARQ Type-I with CC

latency_ARQ_CC = [];

for i=1:num_sim
    latency_ARQ_CC_temp = [];

    for j=1:SNRIdx
        latency_arr = [];

        for l=1:length(repeat_counter_arr_sim{i, j})
            repeat_data = repeat_counter_arr_sim{i, j};

            repeat_data(repeat_data > 3) = 3;
            
            latency_temp = 2 * tti + tti + 2 * (repeat_data(1, l)+1) * 3 * tti + 2*(repeat_data(1, l)+1) * tti;
            latency_arr = [latency_arr, latency_temp];
        end
        latency_ARQ_CC_temp = [latency_ARQ_CC_temp, mean(latency_arr)];
    end
    latency_ARQ_CC = [latency_ARQ_CC; latency_ARQ_CC_temp];
end

% Burst transmission

latency_CC_temp = 2 * tti + tti + 3 * tti + 4 * tti; % Always 3 times repeatition
latency_CC = ones(1, SNRIdx) * latency_CC_temp;

plot(SNRIn, latency_CC, "LineWidth", 2.5, "Color", "blue", 'LineStyle', '-')
hold on
plot(SNRIn, mean(latency_HARQ, 1), 'LineWidth', 1.5, 'Color', 'black', 'LineStyle', '--', 'Marker', '*', 'MarkerSize', 10, 'MarkerIndices',[1:3:length(SNRIn)])
plot(SNRIn, mean(latency_ARQ, 1), "LineWidth", 2, "Color", "magenta", "LineStyle", "--")
plot(SNRIn, mean(latency_ARQ_CC, 1), 'LineWidth', 1.5, 'Color', 'green', 'LineStyle', '--', 'Marker', 'diamond', 'MarkerSize', 10, 'MarkerIndices',[1:3:length(SNRIn)])
grid on
xlabel("SINR [dB]")
ylabel("Average Latency [s]")
legend(["Burst Transmission with CC","HARQ Type III with IR", ...
    "HARQ Type I without CC", "HARQ Type I with CC"], "Location", "best", 'NumColumns',2)
ylim([0 0.04])

ax = gca; 
ax.FontSize = 12;
