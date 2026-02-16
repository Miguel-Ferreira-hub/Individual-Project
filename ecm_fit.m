% This script fits parameters to the equivalent circuit model (ECM) of the
% battery by solving the optimisation problem, and estimates the state of
% charge (SoC) through the predicted open circuit voltage (OCV).
close all;

% Plot ECM fit graphs if set to True
display=true; 

% Data filenames
data_list = [100, 90, 80, 70, 60, 50, 40, 30, 20, 00];

% Labels - true SoCs
labels = [100.0, 90.7, 81.4, 72.0, 62.7, 53.4, 44.1, 34.8, 25.4, 0.00];

for i = 1:length(data_list)
    % Construct filename
    filename = sprintf('%02d.txt', data_list(i));

    % Read data
    data = readtable(filename);
    freq = data.freq__Hz;
    ZI = data.Z1_ohm;
    ZII = data.Z2_ohm;
    Z = sgolayfilt(ZI + 1i*ZII, 3, 11);
    label = sprintf('%.1f%%', labels(i));

    % Fit ECM for each SoC
    [x_opt, RMSEI, RMSEII]=fit(Z, freq, label, display); 
end

% Discharging data
discharge=readtable("Full Discharge.txt"); % Import discharging data
current_discharge=discharge.I_mA; % Current used during discharge over time
time_discharge=discharge.time_s./3600; % Time during discharge
V_discharge=discharge.E_V; % Voltage readings during discharge
capacity_discharge=-trapz(time_discharge, current_discharge); % Calculation of battery capacity
SoC_discharge=((capacity_discharge-(time_discharge.*abs(current_discharge)))./capacity_discharge).*100; % SoC over time for the discharge case

% Charging data
charge=readtable("Full Charge.txt"); % Import charging data
current_charge=charge.I_mA; % Current used during charge over time
time_charge=charge.time_s./3600; % Time during charge
V_charge=charge.E_V; % Voltrage readings during charge
capacity_charge=trapz(time_charge, current_charge); % Calculation of battery capacity
SoC_charge=((0+(time_charge.*abs(current_charge)))./capacity_charge).*100; % SoC over time for the charge case

% Construct OCV array
SoC=[0.00, 16.1, 25.4, 34.8, 44.1, 53.4, 62.7, 72.0, 81.4, 90.7, 100]; % True SoC based on OCV
filenames=linspace(10, 90, 9); % Names of data files
V_array=zeros(1,length(s)+2); % OCV array for each SoC
V_array(1)=V_charge(1); % Adding initial OCV value for 0.00% SoC
V_array(11)=4.008; % Adding final OCV value for 100% SoC

for i = 1:length(filenames)
    % Construct filename
    filename = sprintf('Discharge to %02d.txt', filenames(i));

    % Read data, skipping the first row (header)
    data = readmatrix(filename, 'NumHeaderLines', 1);
    freq = data(:, 1);
    I = data(:, 2);
    V = data(:, 3);
    V_array(i+1)=V(end);
end

% Plot voltage over time for both charging and discharging cases
figure('Name', 'Voltage over time for both cases');
plot(time_discharge, V_discharge, 'DisplayName', 'Discharge' ,'LineWidth', 1);
hold on
plot(time_charge, V_charge, 'DisplayName', 'Charge', 'LineWidth', 1)
xlabel('$Time (hr)$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$Voltage (V)$', 'Interpreter', 'latex', 'FontSize', 15);
legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
box on

% Plot SoC over time for both charging and discharging cases
figure('Name', 'SoC over time for both cases');
plot(SoC_discharge, V_discharge, 'DisplayName', 'Discharge', 'LineWidth', 1);
hold on
plot(SoC_charge, V_charge, 'DisplayName', 'Charge', 'LineWidth', 1);
xlabel('$SoC (\%)$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$Voltage (V)$', 'Interpreter', 'latex', 'FontSize', 15);
legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
xlim([0 100]);
box on

% Terminal voltage array
Vt=[2.75, NaN, 3.423, 3.512, 3.555, 3.59, 3.64, 3.723, 3.82, 3.904, 4.008];

% OCV estimation from fitted ECM
OCV=[Vt(1)+(400e-3)*(sum(x_opt10(1:6))), NaN, Vt(3)+(400e-3)*(sum(x_opt9(1:6))), ... 
    Vt(4)+(400e-3)*(sum(x_opt8(1:6))), Vt(5)+(400e-3)*(sum(x_opt7(1:6))), ...
    Vt(6)+(400e-3)*(sum(x_opt6(1:6))), Vt(7)+(400e-3)*(sum(x_opt5(1:6))), ...
    Vt(8)+(400e-3)*(sum(x_opt4(1:6))), Vt(9)+(400e-3)*(sum(x_opt3(1:6))), ...
    Vt(10)+(400e-3)*(sum(x_opt2(1:6))), Vt(11)+(400e-3)*(sum(x_opt1(1:6)))];

% Estimated vs True SoC-OCV Curves
figure('Name', 'OCV-SoC');
valid_idx = ~isnan(OCV); % Ignore NaN
plot(SoC(valid_idx), OCV(valid_idx), 'DisplayName', 'Estimated OCV', 'LineWidth', 1);
hold on
plot(SoC, V_array, 'DisplayName', 'True OCV', 'LineWidth', 1);
hold off
ylabel('$OCV (V)$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('$SoC (\%)$', 'Interpreter', 'latex', 'FontSize', 15);
legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
box on

% Interpolate SoC-OCV curve with estimated OCV
SoCinterp = interp1(V_array, SoC, OCV, 'spline'); 

SoCRes=((SoC-SoCinterp)./SoC).*100; % Percentage difference between true SoC and estimated SoC
x=linspace(0, 100, 200); % Linear space between 0.00% and 100% SoC
y=spline(SoC, V_array, xq); % Spline interpolation of SoC curve

% Interpolated vs True SoC-OCV Curves
figure; 
plot(SoC, V_array, 'DisplayName', 'Original OCV Curve', 'LineWidth', 1);
hold on
plot(x, y, 'DisplayName', 'Interpolated OCV Curve', 'LineWidth', 1);
ylabel('$OCV (V)$', 'Interpreter', 'latex', 'FontSize', 15);
xlabel('$SoC (\%)$', 'Interpreter', 'latex', 'FontSize', 15);
hold off
legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
box on

% Function to solve optimisation problem and fit ECM
function [x_opt, RMSEI, RMSEII]=fit(Zdata, frequencydata, SoC, display)
    %Frequency
    frequency = logspace(3, -2, 51); % Frequency domain
    frequency=frequency(:); % Convert to column matrix
    omega=frequency*2*pi; % Angular frequency

    %Model Function
    modelFun = @(x, time) fun(x, time); % Function to optimise (ECM impedance)

    %Initial Parameter Guesses
    x0 = [0.01, 0.001, 0.001, 0.02, 0.002, 0.002, 0.01, 0.01, 0.01, 1, 7, 0.002]; 
    lb = [0.001, 0.001, 0.001, 0.001, 0.001, 0.0001, 1e-8, 1e-8, 1e-7, 1, 10, 1e-3]; 
    ub = [0.2, 0.03, 0.03, 0.03, 0.03, 0.01, 1, 1, 1, 1, 2000, 200];

    %Objective Function
    objective = @(x) [(real(modelFun(x, omega)) - real(Zdata));
                  (imag(modelFun(x, omega)) - imag(Zdata))];

    %lsqnonlin algorithm
    options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12);
    x_opt = lsqnonlin(objective, x0, lb, ub, options);

    %Fitted Model
    Z = modelFun(x_opt, omega);

    %RMSE
    RMSEI=[];
    RMSEII=[];
    RMSEIvalue=sqrt(sum((real(Zdata)-real(Z)).^2)./length(frequencydata));
    RMSEIIvalue=sqrt(sum((imag(Zdata)-imag(Z)).^2)./length(frequencydata));
    disp(RMSEIvalue);
    disp(RMSEIIvalue);

    if display==true % Plots
        figure('Name', SoC);
        plot(omega./(2*pi), M(real(Z), imag(Z)), 'DisplayName', 'ECM Fit', 'LineWidth', 1);
        hold on
        plot(frequencydata, M(real(Zdata), imag(Zdata)),'DisplayName', 'Data', 'LineWidth', 1);
        hold off
        ylabel('$Magnitude (dB)$', 'Interpreter', 'latex', 'FontSize', 15);
        xscale('log');
        xlabel('$Frequency (Hz)$', 'Interpreter', 'latex', 'FontSize', 15);
        legend('location', 'northwest', 'FontSize', 10, 'Box', 'off');
        box on

        figure('Name', SoC);
        plot(omega./(2*pi), P(real(Z), imag(Z)), 'DisplayName', 'ECM Fit', 'LineWidth', 1);
        hold on
        plot(frequencydata, P(real(Zdata), imag(Zdata)), 'DisplayName', 'Data', 'LineWidth', 1);
        hold off
        ylabel('$Phase (rad)$', 'Interpreter', 'latex', 'FontSize', 15);
        xscale('log');
        xlabel('$Frequency (Hz)$', 'Interpreter', 'latex', 'FontSize', 15);
        box on

        figure('Name', SoC);
        plot(real(Z), -imag(Z), 'DisplayName', 'ECM Fit', 'LineWidth', 1);
        hold on
        plot(real(Zdata), -imag(Zdata), 'DisplayName', 'Data', 'LineWidth', 1);
        hold off
        ylabel('$-ZII (\Omega)$', 'Interpreter', 'latex', 'FontSize', 15);
        xlabel('$ZI (\Omega)$', 'Interpreter', 'latex', 'FontSize', 15);
        box on
    end
end

% Function for impedance model to be optimised
function Z = fun(x, omega)
    % Matrix with parameters
    R0 = x(1);  
    R1 = x(2);  
    R2 = x(3);
    R3 = x(4);  
    R4 = x(5);  
    R5 = x(6);
    C1 = x(7);  
    C2 = x(8);
    C3 = x(9);  
    C4 = x(10);  
    C5 = x(11);
    sigma = x(12);

    % RC elements in parallel
    Z1 = R1 ./ (1 + 1i * omega * R1 * C1);
    Z2 = R2 ./ (1 + 1i * omega * R2 * C2);
    Z3 = R3 ./ (1 + 1i * omega * R3 * C3);
    Z4 = R4 ./ (1 + 1i * omega * R4 * C4);

    % Warburg + C5 in parallel, in series with R5
    ZC5 = 1 ./ (1i * omega * C5);
    ZW = sigma * (1 - 1i) ./ sqrt(omega);
    Z_parallel = 1 ./ (1 ./ ZC5 + 1 ./ ZW);
    Z5 = R5 + Z_parallel;

    % Total Impedance
    Z = R0 + Z1 + Z2 + Z3 + Z4 + Z5;
end

%Magnitude of impedance
function magnitude=M(ZI, ZII) 
    magnitude=20*log10(((ZI.^2+ZII.^2).^(1/2)));
end

%Phase of impedance
function phase=P(ZI, ZII) 
    phase=atand(ZII./ZI);
end