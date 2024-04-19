function Wireless_Project_Part_1
close all;
clear;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 340;                    %Total number of channels in a cluster
fc = 900;                   %Carrier Frequency (MHz)
h_Tx = 20;                  %Effective Height of Base Station (m)
h_Rx = 1.5;                 %Effective Height of Mobile Station (m)
P_Rx_min = -95;             %Mobile Station Sensitivity (dBm)
Au = 0.025;                 %Traffic intensity per user (Erlang)
n = 4;                      %Path loss exponent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)
    GoS = input("Please enter the grade of service in percentage:\n");
    if (GoS < 0 || GoS > 100)
        fprintf("Please enter a valid value for grade of service\n");
    else
        break;
    end
end
GoS = GoS / 100;            %Grade of Service

City_Area = input("Please enter the city area (in Km^2):\n");
User_density = input("Please enter the user density:\n");
SIR_min = input("Please enter the minimum signal-to-interference ratio (in dB):\n");

while(1)
    Sect_choice = input(['Please choose the sectorization method:\n',...
                    '1. Omni-directional\n', '2. 120 degrees sectorization\n',...
                    '3. 60 degrees sectorization\n']);
    if (Sect_choice ~= 1) && (Sect_choice ~= 2) && (Sect_choice ~= 3)
        fprintf("Please enter a valid choice\n");
    else
        break;
    end
end
%io: number of interfering co-channels
%n_sectors: number of sectors in a cluster
if (Sect_choice == 1) %Omni-directional
    io = 6;
    n_sectors = 1;
elseif (Sect_choice == 2) %120 degrees sectorization
    io = 2;
    n_sectors = 3;
elseif (Sect_choice == 3) %60 degrees sectorization
    io = 1;
    n_sectors = 6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cluster Size Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = cluster_size_calc(SIR_min, n, io); %We use this function to calculate cluster size
fprintf(['Cluster Size = ',num2str(N),'\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Number of cells Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_channels_per_sector = floor(S / (N * n_sectors));
%K = 0:num_channels_per_sector;
%Here, we solve the equation of probability of blocking (Erlang-B) to get 
%the traffic intensity per sector
syms A;
Pr_blocking = (sum((A .^ (0:num_channels_per_sector)) ./ factorial((0:num_channels_per_sector)))).*GoS == ...
                ((A.^num_channels_per_sector) ./ factorial(num_channels_per_sector)) ;
traffic_temp = solve(Pr_blocking,A);
%As the solve function produces symbollic solutions, we should transform
%them into numeric values in order to use them
%The solutions then will be in the form of a column vector whose values are
%sorted from top to bottom and the values could be real of complex
traffic_temp = double(traffic_temp);

%We are looking for the biggest real solution in this vector
if ~isscalar(traffic_temp)
    for index = 1 : length(traffic_temp)
        if isreal(traffic_temp(index))
            traffic_intensity_per_sector = traffic_temp(index);
        end
    end
end

traffic_intensity_per_cell = traffic_intensity_per_sector * n_sectors;

%Total traffic intensity = total number of users * intensity per user
%total number of users = user density * total area
total_traffic_intensity = User_density * City_Area * Au;

num_cells = ceil((total_traffic_intensity) / traffic_intensity_per_cell);

fprintf(['Total Number of Cells = ', num2str(num_cells), ' cells\n']);
fprintf(['Traffic Intensity per Cell = ', num2str(traffic_intensity_per_cell), ' Erlang\n']);
fprintf(['Traffic Intensity per Sector = ', num2str(traffic_intensity_per_sector), ' Erlang\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cell Radius Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cell_Area = City_Area / num_cells;
%As the cell is in hexagonal shape, its area = 1.5 * R^2 * sqrt(3)
Cell_Radius = sqrt(Cell_Area / (1.5 * sqrt(3)));
fprintf(['Cell Radius = ', num2str(Cell_Radius), ' Km\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tx Power Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hata Model for urban medium-sized city
CH = 0.8 + (1.1 * log10(fc) - 0.7) * h_Rx - 1.56 * log10(fc);
%Path loss calculation
Lu = 69.55 + 26.16 * log10(fc) - 13.82 * log10(h_Tx) - ...
    CH + (44.9 - 6.55 * log10(h_Tx)) * log10(Cell_Radius);

%Transmitted power = received power + path loss
P_Tx = P_Rx_min + Lu;
fprintf(['Base Station Transmitted Power = ', num2str(P_Tx), ' dBm\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tx Power Calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance = [(0 : 0.01 : (Cell_Radius)), Cell_Radius, ...
           (Cell_Radius :0.01 : (Cell_Radius + 1))];
Path_loss = zeros(1, length(distance));
P_Received = zeros(1, length(distance));

%Here, we calculate the path loss at every distance in the vector above in
%order to get the received power at these distances

for index = 1 : length(distance)
    Path_loss(index) = 69.55 + 26.16 * log10(fc) - 13.82 * log10(h_Tx) - ...
        CH + (44.9 - 6.55 * log10(h_Tx)) * log10(distance(index));

    P_Received(index) = P_Tx - Path_loss(index);
    if distance(index) == Cell_Radius
        %We just here specify the received power at the cell radius in
        %order to mark it in the plot
        P_Rx_Cell_Radius = P_Received(index);
    end
end
figure(1);
plot(Cell_Radius, P_Rx_Cell_Radius,'*r');
hold on;
plot(distance,P_Received);
title("Recevied Power vs. Distance from Base Station");
xlabel("Distance (Km)");
ylabel("Recevied Power (dBm)");
legend("Received Power @ Cell Radius");
grid on;
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
