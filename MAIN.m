clear all

%{
Calcbands and MAIN written by Christopher Cravey
christopheracravey@gmail.com
Northwestern University, Prof. Matthew Grayson lab 2023
All other code original to Aquila Subband Plotter
written by Dr. Martin Rother,  martin.rother@web.de
%}


%%%Please manually define desired parameters%%%
%Dopant: Delta doping amount in cm^-2 (double)
%FrontGate: Front gate values to be iterated over (array)
%BackGate: Back gate values to be iterated over (array)

Dopant = -4.4E11; 
FrontGate = 0;
BackGate = [-2.5, -2, -1.5, -1, -0.5, -0, 0.5, 1, 1.5, 2];



%GLOBAL SETTINGS
%BackGate = [2, 1.75, 1.5, 1.25, 1.0, 0.75, 0.50, 0.25, 0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -1.75, -2] 
%FrontGate = [1, 0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8, -1]


%2E11 ZOOM SETTINGS
%Dopant = -2E11;
%FrontGate = [ 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2]
%BackGate = [ -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

%2.5E11 ZOOM SETTINGS
%Dopant = -2.5E11;
%FrontGate = [-0.25, -0.20, -0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15, 0.20, 0.25] % array in eV
%BackGate = [] 

%3E11 ZOOM SETTINGS
%Dopant = -3E11;
%FrontGate = [-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0] % array in eV
%BackGate = []


%ETH Top Gate Compare
%Dopant = -1.07E12;
%FrontGate = [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0, -0.1, -0.2, -0.3];
%BackGate = 0;

%ETH Bottom Gate Compare
%Dopant = -4.4E11; 
%FrontGate = 0;
%BackGate = [-2.5, -2, -1.5, -1, -0.5, -0, 0.5, 1, 1.5, 2];



%%%%Please insert desired filenames for output graphs to be saved under%%%%
filename1 = sprintf("%5.3GOccGlobal.png", abs(Dopant));  %Subband Occupation Graph filename, remove "-" from front
filename2 = sprintf("%5.3GConcGlobal.png", abs(Dopant)); %Equi-Electron Density Graph filename, remove "-" from front



%%%%%%%%%%%%%%%%  PERFORM CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%
Data = calcbands(Dopant, BackGate, FrontGate); %Run simulation

%{
%Optional: save the final iteration graph which shows conduction band,
%carrier density, and wavefunctions of lowest 2 subbands
filename = ''
saveas(gcf, '')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Calculate which subbands are occupied
subband_occ = zeros([length(BackGate), length(FrontGate)]);
for i = 1:length(BackGate)
    for j = 1:length(FrontGate)
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

%saveas(gcf, filename1)





%%%% Plot a color map of electron density in QW as a function of gate voltage  %%%%
figure

contourf(Data.Vbot, Data.Vtop, Data.WellConc, "ShowText", true, "LabelFormat", "%0.3G")
xlabel("V_b_o_t (V)")
ylabel("V_t_o_p (V)")
xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
title("Carrier Concentration (cm^-^2)")


%saveas(gcf, filename2)

%}



%{
%%Plot carrier conc VS top gate potential

figure

scatter(Data.Vtop, -1 * Data.WellConc, "filled", "blue")
xlabel("V_t_o_p (V)")
ylabel("n (cm^-2)")
%xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
%ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
title("Carrier Concentration VS Top Gate Bias")


%Plot ETH top gate results from sample D170301B for comparison
hold on
scatter(linspace(-0.3, 0.7, 101), 1E11*[0.3599
0.39858
0.44695
0.4957
0.54377
0.59188
0.6412
0.69061
0.73893
0.78821
0.83792
0.88561
0.93511
0.98418
1.03222
1.08464
1.13147
1.18068
1.22791
1.27912
1.32871
1.37051
1.41852
1.47975
1.51588
1.56743
1.62104
1.67656
1.72253
1.76949
1.82914
1.86717
1.9199
1.97216
2.01719
2.06647
2.12664
2.17181
2.22231
2.27569
2.32374
2.37563
2.43069
2.47017
2.53898
2.59108
2.63232
2.66041
2.71166
2.76536
2.82922
2.87381
2.92487
2.96422
3.01205
3.07302
3.11471
3.21302
3.24237
3.26248
3.30529
3.37254
3.4253
3.45067
3.51188
3.52298
3.58562
3.64234
3.66282
3.7006
3.72429
3.74907
3.71879
3.71319
3.74134
3.74156
3.75228
3.76663
3.79862
3.78339
3.77799
3.78407
3.84253
3.82767
3.8197
3.84945
3.84626
3.83997
3.87822
3.85581
3.87663
3.89598
3.92347
3.90118
3.88417
3.8847
3.88476
3.93672
3.89823
3.92152
3.95151
], "filled", "red");
legend(["NU Simulation", "ETH"]);
%}



%%Plot carrier conc VS bottom gate potential

figure

scatter(Data.Vbot, Data.WellConc, "filled", "blue");
xlabel("V_b_o_t (V)");
ylabel("n (cm^-2)");
%xlim([min(Data.Vbot, [], "all") max(Data.Vbot, [], "all")])
%ylim([min(Data.Vtop, [], "all") max(Data.Vtop, [], "all")])
title("Carrier Concentration VS Bottom Gate Bias");

%Plot back gate eth experimental results
hold on;
scatter([-2.5, 0, 2], [-1.8E11, -2.9E11, -3.7E11], "filled", "red");
legend(["NU Simulation", "ETH"]);




