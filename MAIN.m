clear all

%Please manually insert desired parameters


Dopant = -4.4E11   % Single value in cm^-2. This dopant conc is applied to both sides
%Field = [0] 
Field = linspace(-2.5, 2.5, 11)
%Pot = linspace(-0.3, 0.7, 11)
Pot = [0]

%Pot = [ 0, 0.02, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2] %Potential array at 0A
%Field = [ 7000, 5000, 2000, 0, -10000, -50000] %E field array at bottom gate in V/cm



%GLOBAL SETTINGS
%Field = [20000, 15000, 10000, 5000, 1000, 0, -100000, -200000] 
%Pot = [1, 0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8, -1]


%2E11 ZOOM SETTINGS
%Pot = [ 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2]
%Field = [ 7000, 6500, 6000, 5500, 5000, 4500, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 500, 0, -1000, -10000, -30000, -50000] %2E11 zoom


%2.5E11 ZOOM SETTINGS
%Pot = [-0.25, -0.20, -0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15, 0.20, 0.25] % array in eV
%Field = [15000, 13500, 11500, 10000, 8500, 6500, 5000] % array in V/cm


%3E11 ZOOM SETTINGS
%Pot = [-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0] % array in eV
%Field = [20000, 18500, 16500, 15000, 13500, 11500, 10000] % array in V/cm


%Please insert desired filenames for output graphs to be saved under
filename1 = sprintf("%5.3GOccupation.png", abs(Dopant));  %Subband Occupation Graph filename, remove "-" from front
filename2 = sprintf("%5.3GEqui_Density.png", abs(Dopant)); %Equi-Electron Density Graph filename, remove "-" from front

%Convert field from V/cm to V/A
%Field = Field*(1E-8);


%%%%%%%%%%%%%%%%  PERFORM CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%
Data = calcbands(Dopant, Field, Pot); %Run simulation


%%%Optional: save the final iteration graph which shows conduction band,
%carrier density, and wavefunctions of lowest 2 subbands
%filename = ''
%saveas(gcf, '')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Calculate which subbands are occupied
subband_occ = zeros([length(Field), length(Pot)]);
for i = 1:length(Field)
    for j = 1:length(Pot)
        subband_occ(i,j) = (Data.Sub1.Occ(i,j) + Data.Sub2.Occ(i,j) + Data.Sub3.Occ(i,j));
    end
end



figure %Create new graph window

colormap cool; %Make the plot look nice:
%3 subbands occupied = ?
%2 subbands occupied = Pink 
%1 subbands occupied = Purple
%0 subbands occupied = Light Blue


%{
%%%%  Plot subbands occupancy  %%%%
%Spread the matrix data into vectors so they may be plotted easily
scatter(Data.Vbot(:), Data.Vtop(:), [], subband_occ(:), 'filled')
xlabel("V_b_o_t (V)")
ylabel("V_t_o_p (V)")
xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
D = sprintf('%5.3G', Dopant);  %Express doping conc. in scientific notation
                      
title("Delta Doping: " + D)


%Create some invisable points to easily make 
%a constant legend for the scatter plot
hold on
h=gobjects(3,1);
h(1)=scatter(nan,nan, [], 0,'filled');
h(2)=scatter(nan,nan,[], 1,'filled');
h(3)=scatter(nan,nan,[], 2,'filled');
legend(h, ["0" "1" "2"])
lgd = legend;
lgd.Title.String = "# of Occupied Subbands";
hold off
%}

%saveas(gcf, filename1)




%{
%%%% Plot a color map of electron density in QW as a function of gate voltage  %%%%
figure

contourf(Data.Vbot, Data.Vtop, Data.Well.Conc, "ShowText", true, "LabelFormat", "%0.3G")
xlabel("V_b_o_t (V)")
ylabel("V_t_o_p (V)")
xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
title("Carrier Concentration (cm^-^2)")
%}

%{
scatter(Data.Vtop, Data.Well.Conc, "filled")
xlabel("V_t_o_p (V)")
ylabel("n (cm^-2)")
%xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
%ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
title("Carrier Concentration VS Top Gate Bias")
%}


scatter(Data.Vbot, Data.Well.Conc, "filled")
xlabel("V_b_o_t (V)")
ylabel("n (cm^-2)")
%xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
%ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
title("Carrier Concentration VS Bottom Gate Bias")


%{
%Plot line corresponding to subband degeneracy crossover
%See below for reference line of best fits

hold on
plot([-0.83, 0.05], [-0.083, -0.107], "r", "LineWidth", 2)
hold off
%}

%LINE OF SUBBAND CROSSOVER FOR 2E11 Delta Doping
% X = [-0.83, -0.78, -0.71, -0.6, -0.49, -0.39, -0.28, -0.17, -0.06, 0.05]
% Y = [-0.087, -0.087, -0.087, -0.087, -0.087, -0.097, -0.097, -0.097,
% -0.107, -0.107]
% LINE OF BEST FIT: plot([-0.83, 0.05], [-0.083, -0.107], "r", "LineWidth", 2)

%saveas(gcf, filename2)


%TO DO: Overlay contour over subband occupancy
%Fix x and y axis of all plots to easily overlap each other
%Look into "Hold on" function