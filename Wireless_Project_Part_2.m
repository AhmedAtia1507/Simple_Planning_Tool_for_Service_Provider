function Wireless_Project_Part_2
close all;
clear;clc;
%% Cluster Size vs. SIR plot
SIR = 1:0.01:30;
io = [6 2 1];
n = 4;
N = zeros(length(io), length(SIR));
figure;
%Here, we calculate the cluster size for every value of SIR and for
%different sectorization method
%io = 6 for omni-directional
%io = 2 for 120-degrees sectorization
%io = 1 for 60-degrees sectorization
for index_1 = 1 : length(io)
    for index_2 = 1 : length(SIR)
        N(index_1, index_2) = cluster_size_calc(SIR(index_2), n, io(index_1));
    end
    plot(SIR, N(index_1,:))
    hold on;
end
title("Cluster Size vs. Minimum Signal-to-Interference Ratio");
xlabel("SIR (dB)");
ylabel("Cluster Size");
grid on;
legend('Omni-Directional', '120 degrees sectorization', '60 degrees sectorization');
%% 
%Here, we plot the number of cells and traffic intensity per cell at 
%SIR = 19 dB and user density = 1400 users/km^2
plot_part_2_3(19);
%%
%Here, we plot the number of cells and traffic intensity per cell at 
%SIR = 14 dB and user density = 1400 users/km^2
plot_part_2_3(14);
%%
%Here, we plot the number of cells and cell radius at SIR = 14 dB
plot_part_4_5(14);
%%
%Here, we plot the number of cells and cell radius at SIR = 19 dB
plot_part_4_5(19);
end

function N = cluster_size_calc(SIR_min,n,io)
N_calculated = (1 / 3) .* ((((io .* (10 ^ (SIR_min / 10))) .^ (1 / n)) + 1) .^ 2);
%N_calculated = (1 / 3) .* ((((io .* (10 ^ (SIR_min / 10))) .^ (1 / n)) ).^ 2); 
N_values = zeros(1, (floor(N_calculated) * floor(N_calculated)));
N_index = 1;
%Here, we get the allowed values of the cluster size that follows the
%equation N = i^2 + ik + k^2
for i = 0 : floor(N_calculated)
    for k = 0 : floor(N_calculated)
        N_values(N_index) = (i ^ 2) + (i * k) + (k ^ 2);
        N_index = N_index + 1;
    end
end
%Now, we remove the values that are less than N_calculated
N_values(N_values < N_calculated) = [];
%Then, the minimum value in the rest is the cluster size we are looking for
N = min(N_values);
end

function plot_part_2_3(SIR)
n = 4; City_Area = 100; S = 340; Au = 0.025; 
User_density = 1400;
GoS = 0.01 : 0.01 : 0.3;
traffic_intensity_per_cell = zeros(3, length(GoS));
num_cells = zeros(3, length(GoS));
io = [6 2 1]; n_sectors = [1 3 6];

%In order to plot the number of cells and traffic intensity per cell, we
%calculate them at each value of GoS and for each sectorization method
for index_1 = 1 : length(io)
    N = cluster_size_calc(SIR, n, io(index_1));
    num_channels_per_sector = floor(S / (N * n_sectors(index_1)));
    
    for index_2 = 1 : length(GoS)
        syms A;
        Pr_blocking = (sum((A .^ (0:num_channels_per_sector)) ./ factorial((0:num_channels_per_sector)))).*GoS(index_2) == ...
                ((A.^num_channels_per_sector) ./ factorial(num_channels_per_sector)) ;
        traffic_temp = solve(Pr_blocking,A);
        traffic_temp = double(traffic_temp);

        if ~isscalar(traffic_temp)
            for index = 1 : length(traffic_temp)
                if isreal(traffic_temp(index))
                    traffic_intensity_per_sector = traffic_temp(index);
                end
            end
        end
        traffic_intensity_per_cell(index_1,index_2) = traffic_intensity_per_sector * n_sectors(index_1);

        total_traffic_intensity = User_density * City_Area * Au;

        num_cells(index_1,index_2) = ceil((total_traffic_intensity) / traffic_intensity_per_cell(index_1,index_2));
    end
end
    figure
    for index_1 = 1 : length(io)
        subplot(1,2,1)
        plot(GoS,num_cells(index_1,:));
        hold on;
        subplot(1,2,2)
        plot(GoS, traffic_intensity_per_cell(index_1,:))
        hold on;
    end
    subplot(1,2,1);
    title(['Number of Cells vs. Grade of Service at SIR = ', num2str(SIR),' dB']);
    xlabel("GoS");
    ylabel("Number of Cells");
    legend('Omni-Directional', '120 degrees sectorization', '60 degrees sectorization');
    grid on;
    subplot(1,2,2);
    title(['Traffic Intensity per Cell vs. Grade of Service at SIR = ', num2str(SIR),' dB']);
    xlabel("GoS");
    ylabel("Traffic Intensity per Cell (Erlang)");
    legend('Omni-Directional', '120 degrees sectorization', '60 degrees sectorization');
    grid on;
end

function plot_part_4_5(SIR)
n = 4; City_Area = 100; S = 340; Au = 0.025; 
User_density = 100 : 2000;
GoS = 0.02;
num_cells = zeros(3, length(User_density));
Cell_Radius = zeros(3,length(User_density));
io = [6 2 1]; n_sectors = [1 3 6];

for index_1 = 1 : length(io)
    N = cluster_size_calc(SIR, n, io(index_1));
    num_channels_per_sector = floor(S / (N * n_sectors(index_1)));
    syms A;
    Pr_blocking = (sum((A .^ (0:num_channels_per_sector)) ./ factorial((0:num_channels_per_sector)))).*GoS == ...
                ((A.^num_channels_per_sector) ./ factorial(num_channels_per_sector)) ;
    traffic_temp = solve(Pr_blocking,A);
    traffic_temp = double(traffic_temp);

    if ~isscalar(traffic_temp)
        for index = 1 : length(traffic_temp)
            if isreal(traffic_temp(index))
                traffic_intensity_per_sector = traffic_temp(index);
            end
        end
    end
    traffic_intensity_per_cell = traffic_intensity_per_sector * n_sectors(index_1);
    for index_2 = 1 : length(User_density)

        total_traffic_intensity = User_density(index_2) * City_Area * Au;

        num_cells(index_1,index_2) = ceil((total_traffic_intensity) / traffic_intensity_per_cell);
        
        Cell_Area = City_Area / num_cells(index_1,index_2);
        Cell_Radius(index_1,index_2) = sqrt(Cell_Area / (1.5 * sqrt(3)));
    end
end
    figure
    for index_1 = 1 : length(io)
        subplot(1,2,1)
        plot(User_density,num_cells(index_1,:));
        hold on;
        subplot(1,2,2)
        plot(User_density, Cell_Radius(index_1,:));
        hold on;
    end
    subplot(1,2,1);
    title(['Number of Cells vs. User Density at SIR = ', num2str(SIR),' dB']);
    xlabel("User Density");
    ylabel("Number of Cells");
    legend('Omni-Directional', '120 degrees sectorization', '60 degrees sectorization');
    grid on;
    subplot(1,2,2);
    title(['Cell Radius vs. User Density at SIR = ', num2str(SIR),' dB']);
    xlabel("User Density");
    ylabel("Cell Radius (Km)");
    legend('Omni-Directional', '120 degrees sectorization', '60 degrees sectorization');
    grid on;
end