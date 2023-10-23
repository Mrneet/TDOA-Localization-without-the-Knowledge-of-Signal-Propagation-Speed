% Generate Fig. 1-2

% Y. Sun, K. C. Ho, Y. Yang, and L. Chen, "An Asymptotically Optimal
% Estimator for Source Location and Propagation Speed by TDOA," IEEE Signal
% Processing Letters, vol. 30, pp. 1037-1041, Aug. 2023.

% Created by Y. Sun, K. C. Ho, Oct. 2023

clear all; 
close all; clc;

rng('default');
warning off

% ----- simulation setting -----
noiseOrRange = 'nse';
% noiseOrRange = 'err';

% sensor position
senPos = round([
               0       0       0
           46.64  -87.12   62.94
          124.02   -7.98   81.16
           44.37   63.23  -74.60
           57.98  -23.62   82.68
          105.02  -51.72   26.47
          -81.56  104.48  -80.49
          -19.51  -18.06   28.50
          -41.64  -42.63   26.96
           24.02  -17.98   61.81
           ]'*2);

srcLoc = [1346;547;373];

c = 1500;

switch noiseOrRange
    case 'nse'
        % ******* vs. noise power config *******
        nsePwrDB = -90:5:-30;
        errLvlDB = -50;
        NumEnsembles = 1000;
    case 'err'
        % ******* vs. sensor error config *******
        nsePwrDB = -60;  % 10log(rad^2)
        errLvlDB = -40:5:20;
        NumEnsembles = 1000;
end

[N,M] = size(senPos);
K = length(nsePwrDB);       % number of noise levels
R = length(errLvlDB);       % number of ranges

nse = randn(M-1,NumEnsembles);
err = randn(N,M,NumEnsembles);
nse = nse - mean(nse,2);
err = err - mean(err,3);

disp('Simulation is running ...');
totalTime = zeros(6,1); boxdata = cell(0);

r = sqrt(sum((repmat(srcLoc,1,M)-senPos).^2,1))';
rd = r(2:end) - r(1);
td = rd/c;

aa = [1,3,7,10,4,1,9,7,2,1,3];
SS = kron(diag(aa(1:M)),eye(N));
nAlg = 6;

for is = 1:R   % loop through ranges
    disp(['Sensor position error: ',num2str(errLvlDB(is)),', ',num2str(is),'/',num2str(R),' ...']);
    
    Qs = 10^(errLvlDB(is)/10) * SS;
    
    for in = 1:K    % loop through noise powers
        disp(['Noise power (10log(\sigma^2): ',num2str(nsePwrDB(in)),', ',num2str(in),'/',num2str(K),' ...']);
        Qt = 10^(nsePwrDB(in)/10) * (ones(M-1, M-1)+eye(M-1))/2;
        
        % Calculate CRLB
        CRBcart = CRLB_TDOA_UPS(senPos, srcLoc, c, Qt,Qs);
        crlb_u(is,in) = trace(CRBcart(1:N,1:N));
        crlb_c(is,in) = CRBcart(end,end);
        
        COV = COV_Proj_UPS(senPos,srcLoc,c,Qt,Qs);
        cov_u(is,in) = trace(COV(1:N,1:N));
        cov_c(is,in) = COV(end,end);
        
        % -- Obtaining source location estimate --
        Qs1 = Qs(1:N:end,1:N:end);
        pos = zeros(N,NumEnsembles);
        vel = zeros(NumEnsembles,1);

        for i = 1:NumEnsembles
            td_m = td + sqrtm(Qt)*nse(:,i);
            senPos_m = senPos + err(:,:,i)*sqrtm(Qs1);

            [sol_u,sol_v] = TDOALoc_Proj_UPS(senPos_m,td_m,Qt,Qs);
            pos(:,i) = sol_u;
            vel(i) = sol_v;
        end

        mse_u(in,is) = mean(sum((pos - srcLoc).^2,1));
        mse_c(in,is) = mean((vel - c).^2);

        avBia_u(in,is) = norm(mean(pos - srcLoc,2));
        avBia_c(in,is) = abs(mean(vel - c));
    end
end
totalTime

switch noiseOrRange
    case 'nse'
        xlabtext = '$10log(\sigma_t^2(s^2))$';
        xdata = nsePwrDB;
        yl_mse = [5, 85;
            -20,80];
        yl_bias = [-40,63;
            -100,62];
    case 'err'
        xlabtext = '$10log(\sigma_s^2(m^2))$';
        xdata = errLvlDB;
        yl_mse = [35,75;
            -10,60];
        yl_bias = [20,70;
            -20,50];
end
    
% MSE of position estimate
figure;
subplot(2,1,1)
plot(xdata, 10*log10(reshape(mse_u,[],1)), 'o', 'LineWidth', 1.5, 'DisplayName', 'DNSP');hold on;grid on;
plot(xdata, 10*log10(crlb_u), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
plot(xdata, 10*log10(cov_u), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-Cov');
xlabel(xlabtext,'Interpreter','latex', 'FontSize', 13);
ylabel('$10log(MSE(\hat{\bf u})(m^2))$','Interpreter','latex', 'FontSize', 13);
ylim(yl_mse(1,:));
lgd11 = legend('Show');
set(lgd11, 'FontSize',11, 'Location', 'Northwest');

% MSE of speed estimate
subplot(2,1,2)
plot(xdata, 10*log10(reshape(mse_c,[],1)), 'o', 'LineWidth', 1.5, 'DisplayName', 'DNSP');hold on;grid on;
plot(xdata, 10*log10(crlb_c), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
plot(xdata, 10*log10(cov_c), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-Cov');
xlabel(xlabtext,'Interpreter','latex', 'FontSize', 13);
ylabel('$10log(MSE(\hat{c})(m^2/s^2))$','Interpreter','latex', 'FontSize', 13);
ylim(yl_mse(2,:));
lgd2 = legend('Show');
set(lgd2, 'FontSize',11, 'Location', 'Northwest');


% Bias of angle estimate
figure;
subplot(2,1,1)
plot(xdata, 20*log10(reshape(avBia_u,[],1)), 'o', 'LineWidth', 1.5, 'DisplayName', 'DNSP');hold on;grid on;
xlabel(xlabtext,'Interpreter','latex', 'FontSize', 13);
ylabel('$20log(Bias(\hat{\bf u})(m))$','Interpreter','latex', 'FontSize', 13);
ylim(yl_bias(1,:));
h3 = legend('Show');
set(h3, 'FontSize',11, 'Location', 'Northwest');

% Bias of c estimate
subplot(2,1,2)
plot(xdata, 20*log10(reshape(avBia_c,[],1)), 'o', 'LineWidth', 1.5, 'DisplayName', 'DNSP');hold on;grid on;
xlabel(xlabtext,'Interpreter','latex', 'FontSize', 13);
ylabel('$20log(Bias(\hat{c})(m/s))$','Interpreter','latex', 'FontSize', 13);
ylim(yl_bias(2,:));
h3 = legend('Show');
set(h3, 'FontSize',11, 'Location', 'Northwest');



