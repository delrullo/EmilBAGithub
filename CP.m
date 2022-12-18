clc, clear, close all

%% Data import
demand = importdata('C:\Users\emilb\Desktop\Vedvarende Energinetværk\CP\electricity_demand.csv');
solarPV_CF = importdata('C:\Users\emilb\Desktop\Vedvarende Energinetværk\CP\pv_optimal.csv');
wind_onshore_CF = importdata('C:\Users\emilb\Desktop\Vedvarende Energinetværk\CP\onshore_wind_1979-2017.csv');

    % Removing countries
% solarPV - Malta(:,23), Siberia{SERBIA}(:,29), Cyprus(:,6), and Bosnia and Herzegovina(:,4)

all_solar_CF = solarPV_CF.data(:,[1:3,5,7:22,24:28,30:end]);

% demand -  Siberia{SERBIA}(:,28), Cyprus(:,6), and Bosnia and Herzegovina(:,4)

load = demand.data(:,[1:3,5,7:27,29:end]); % year 2015

% onshore wind - No removals

all_wind_CF = wind_onshore_CF.data;


    % Isolating year 2015 for wind and solar
% 2015 is index 315577 to 324336

solar_CF = all_solar_CF(315577:324336,:);
wind_CF = all_wind_CF(315577:324336,:);

%% Power Generation

solar_gen = solar_CF .* load;
wind_gen = wind_CF .* load;


%% Normalizing data
% norm_wind = zeros(8760,28);
% norm_solar = zeros(8760,28);

% for i = 1:width(load) % country
%     for j = 1:height(load) % hour
%         
%         norm_wind(j,i) = (wind_gen(j,i) .* mean(load(:,i))) ./ mean(wind_gen(:,i));
%         norm_solar(j,i) = (solar_gen(j,i) .* mean(load(:,i))) ./ mean(solar_gen(:,i));
%     
%     end
% end

norm_wind = (wind_gen .* mean(load)) ./ mean(wind_gen);
norm_solar = (solar_gen .* mean(load)) ./ mean(solar_gen);
%% Mismatch calculation

a = 0:0.2:1; 
alpha = ["a1",'a2','a3','a4','a5','a6'];
delta = struct(alpha(1),zeros(8760,28),alpha(2),zeros(8760,28),alpha(3),zeros(8760,28),alpha(4),zeros(8760,28),alpha(5),zeros(8760,28),alpha(6),zeros(8760,28));


for i = 1:width(load) % country
    for j = 1:length(a) % alpha value

        delta.(alpha{j})(:,i) = a(j) .* norm_wind(:,i) + (1 - a(j)) .* norm_solar(:,i) - load(:,i);

    end
end

%% Mismatch for Germany(:,6) and Denmark(:,7)

sub_title = ["$\alpha = 0.0 $" '$\alpha = 0.2$' '$\alpha = 0.4 $' '$\alpha = 0.6 $' '$\alpha = 0.8 $' '$\alpha = 1.0 $'];
    
    % Denmark(:,7)
figure(1)
t1_DK = tiledlayout(2,3);
title(t1_DK,'Mismatch for Denmark','Interpreter','latex')

for i = 1:length(a)
    nexttile
    plot(delta.(alpha{i})(:,7))
    title(sub_title{i},'Interpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    xlabel('Time [h]')
    ylabel('$\Delta$',Rotation=360)
    xlim([0 8760]);
    grid on
end

print('Mismatch_DNK', '-depsc');  

    % Germany(:,6)
figure(2)
t1_DE = tiledlayout(2,3);
title(t1_DE,'Mismatch for Germany','Interpreter','latex')

for i = 1:length(a)
    nexttile
    plot(delta.(alpha{i})(:,6))
    title(sub_title{i},'Interpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    xlabel('Time [h]')
    ylabel('$\Delta$',Rotation=360)
    xlim([0 8760]);
    grid on
end

print('Mismatch_DEU', '-depsc');  

%% Aggregated European mismatch


deltaEU = struct(alpha(1),zeros(8760,1),alpha(2),zeros(8760,1),alpha(3),zeros(8760,1),alpha(4),zeros(8760,1),alpha(5),zeros(8760,1),alpha(6),zeros(8760,1));

for j = 1:length(a) % alpha value
    for i = 1:width(load) % country
    
        deltaEU.(alpha{j}) = delta.(alpha{j})(:,i) + deltaEU.(alpha{j});
  
    end
end

%% Plot of aggregated European mismatch and comparison to individual countries

figure(3)
t1_DK = tiledlayout(2,3);
title(t1_DK,'Mismatch for Europe','Interpreter','latex')

for i = 1:length(a)
    nexttile
    plot(deltaEU.(alpha{i}),'b')
    hold on
    plot(delta.(alpha{i}),'r')
    title(sub_title{i},'Interpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    xlabel('Time [h]')
    ylabel('$\Delta$',Rotation=360)
    xlim([0 8760]);
%     legend('European-wide mismatch', 'Individual country mismatch')
    grid on
end

print('Mismatch_EU', '-depsc'); 

%% Mismatch variance

var_delta_DK = zeros(length(a),1);
var_delta_DE = zeros(length(a),1);
var_delta_EU = zeros(length(a),1);


    % Determine the variance of Denmark and Germany
for i = 1:length(a)

    var_delta_DK(i) = var(delta.(alpha{i})(:,7));
    var_delta_DE(i) = var(delta.(alpha{i})(:,6));
    var_delta_EU(i) = var(deltaEU.(alpha{i}));

end

% divide by mean(load)^2

var_load_delta_DK = var_delta_DK / mean(load(:,7))^2;
var_load_delta_DE = var_delta_DE / mean(load(:,6))^2;
var_load_delta_EU = var_delta_EU / mean(mean(load))^2;

%% Plot variance of mismatches DK, DE, EU
    % DK
figure(4)
plot(a,var_load_delta_DK,'.-',MarkerSize=20,LineWidth=2)

title('Mismatch Variance for Denmark')
xlabel('$\alpha$')
ylabel('Var$(\vec{\Delta}) / \langle L_{\mbox{DNK}}\rangle ^2$')
ylim([0 max(var_load_delta_DK)+0.5])
grid
print('Mismatch_variance_DNK', '-depsc'); 

    % DE
figure(5)
plot(a,var_load_delta_DE,'.-',MarkerSize=20,LineWidth=2)

title('Mismatch Variance for Germany')
xlabel('$\alpha$')
ylabel('Var$(\vec{\Delta}) / \langle L_{\mbox{DEU}}\rangle ^2$')
ylim([0 max(var_load_delta_DE)+0.5])
grid
print('Mismatch_variance_DEU', '-depsc');  

    % EU
    figure(6)
plot(a,var_load_delta_EU,'.-',MarkerSize=20,LineWidth=2)

title('Mismatch Variance for Europe')
xlabel('$\alpha$')
ylabel('Var$(\vec{\Delta}) / \langle L_{\mbox{EU}}\rangle ^2$')
ylim([0 max(var_load_delta_EU)+0.5])
grid
print('Mismatch_variance_EU', '-depsc'); 

%% Covariance matrix and eigen values

delta_cov = struct(alpha(1),zeros(28,28),alpha(2),zeros(28,28),alpha(3),zeros(28,28),alpha(4),zeros(28,28),alpha(5),zeros(28,28),alpha(6),zeros(28,28));

    % Determine covarince matrices for each alpha 0:02:1

for i = 1:length(a)

    delta_cov.(alpha{i}) = cov(delta.(alpha{i}));

end

    % Determine eigenvalues and eigenvectors

for i = 1:length(a)

    eigen.vec.(alpha{i}) = 0;
    eigen.val_mat.(alpha{i}) = 0;
    [eigen.vec.(alpha{i}),eigen.val_mat.(alpha{i})] = eig(delta_cov.(alpha{i}));

    % Extract eigenvalues
    eigen.val.(alpha{i}) = diag(eigen.val_mat.(alpha{i}));

    % Sort in decending order
    eigen.val.(alpha{i}) = flip(eigen.val.(alpha{i}));
    eigen.vec.(alpha{i}) = flip(eigen.vec.(alpha{i}),2);

    % Cumulative sum of eigenvalues
    eigen.cum_val.(alpha{i}) = cumsum(eigen.val.(alpha{i}));

    % Set sum to 1
    eigen.cum_val_1.(alpha{i}) = eigen.cum_val.(alpha{i}) / max(eigen.cum_val.(alpha{i}));
    eigen.val_1.(alpha{i}) = eigen.val.(alpha{i}) / max(eigen.val.(alpha{1}));

end


eigen_val_mat = zeros(length(a),length(a));
eigen_val_mat_cum = zeros(length(a),length(a));
eigen_val_mat_cum_tot = zeros(width(load),length(a));
eigen_val_mat_cum_unscale = zeros(length(a),length(a));

for i = 1:length(a)
   
    % Create matrix with the six largest eigenvalues for each alpha increment
    eigen_val_mat(:,i) = eigen.val.(alpha{i})(1:6);

    % Cumulative eigenvalues
    eigen_val_mat_cum_unscale(:,i) = eigen.cum_val.(alpha{i})(1:6);
    eigen_val_mat_cum_tot(:,i) = eigen.cum_val.(alpha{i});

    % Cumulative eigenvalues scaled
    eigen_val_mat_cum(:,i) = eigen_val_mat_cum_unscale(:,i) ./ max(eigen_val_mat_cum_tot(:,i));

end

    % Determine the change between principal components
mat(:,1) = eigen_val_mat_cum(1,:);
mat(:,2) = eigen_val_mat_cum(2,:) -  eigen_val_mat_cum(1,:);
mat(:,3) = eigen_val_mat_cum(3,:) - eigen_val_mat_cum(2,:);
mat(:,4) = eigen_val_mat_cum(4,:) - eigen_val_mat_cum(3,:);
mat(:,5) = eigen_val_mat_cum(5,:) - eigen_val_mat_cum(4,:);
mat(:,6) = eigen_val_mat_cum(6,:) - eigen_val_mat_cum(5,:);
mat(:,7) = 1 - eigen_val_mat_cum(6,:);

     % impact of eigenvalues on system
figure(7)
area(a,mat)

title('PCA eigenvalues')
xlabel('$\alpha$')
ylabel('$\lambda_k$')
legend('$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$','$7 \leq k \leq 32$','Interpreter','latex',Location='best')
print('PCA_eigenvalue', '-depsc');

%%

% country names
cou_nam = ({'Belgium', 'Bulgaria', 'Czech Republic', 'Denmark', 'Germany','Estonia','Ireland',...
    'Greece','Spain','France','Croatia','Italy', 'Latvia', 'Lithuania', 'Luxembourg',...
     'Hungary', 'Malta', 'Netherlands','Austria', 'Poland', 'Portugal','Romania','Slovenia', 'Slovakia',...
        'Finland', 'Sweden', 'Malta'});
numRegions = length(cou_nam);

PC_nam = ["$\lambda_1$" "$\lambda_2$" "$\lambda_3$" "$\lambda_4$" "$\lambda_5$" "$\lambda_6$"];
set(0,'defaultTextInterpreter','latex');
%% Alpha = 0.0

figure(8)
tl_a1 = tiledlayout(1,3);
cmap = jet(30); % Number of color increments
colormap(cmap)

for j = 1:3
    % Plots the first three principal components    
    nexttile
    
    % Set color value for range [-1 1]
    color_val_a1 = round(interp1([-1 1],[1 30],eigen.vec.a1(:,j)));
    
    worldmap("Europe");
    geoshow('landareas.shp');
    bordersm('countries','facecolor','white');


    for i = 1:numRegions
        bordersm(cou_nam{i},'facecolor',cmap(color_val_a1(i),:)); 
    end

title(PC_nam(j));

end

colorbar
set(gcf, 'Position', get(0, 'Screensize'));
print('Principal_components_alpha_00', '-depsc');

%% Alpha = 1.0

figure(9)
tl_a6 = tiledlayout(2,3);
colormap(cmap)

for j = 1:6
% Plots the first three principal components
nexttile

% Set color value for range [-1 1]
color_val_a6 = round(interp1([-1 1],[1 30],eigen.vec.a6(:,j)));

worldmap("Europe");
geoshow('landareas.shp');
bordersm('countries','facecolor','white'); 

for i = 1:numRegions
        bordersm(cou_nam{i},'facecolor',cmap(color_val_a6(i),:)); 
end

title(PC_nam(j));

end

colorbar
set(gcf, 'Position', get(0, 'Screensize'));
print('Principal_components_alpha_1', '-depsc');  
