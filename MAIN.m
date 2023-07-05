clear all

%Please manually insert desired parameters


Dopant = -2.5E11   % Single value in cm^-2. This dopant conc is applied to both sides
Pot = [0.25, 0.20, 0.15, 0.10, 0.05, 0, -0.05, -0.10, -0.15, -0.20, -0.25] % array in eV
Field = [15000, 13500, 11500, 10000, 8500, 6500, 5000] % array in V/cm

%GLOBAL SETTINGS
%Field = [20000, 15000, 10000, 5000, 1000, 0, -100000, -200000] 
%Pot = [1, 0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8, -1]


%2E11 ZOOM SETTINGS
%Field = [ 7000, 5000, 2000, 0, -10000, -50000] %2E11 zoom
%Pot = [ 7000, 5000, 3500, 2000, 0, -50000, -100000]


%2.5E11 ZOOM SETTINGS
%Pot = [0.25, 0.20, 0.15, 0.10, 0.05, 0, -0.05, -0.10, -0.15, -0.20, -0.25] % array in eV
%Field = [15000, 13500, 11500, 10000, 8500, 6500, 5000] % array in V/cm


%3E11 ZOOM SETTINGS
%Pot = [0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0] % array in eV
%Field = [20000, 18500, 16500, 15000, 13500, 11500, 10000] % array in V/cm


%Please insert desired filenames for output graphs to be saved under
filename1 = sprintf("%5.3GOccupationZoom.png", abs(Dopant));  %Subband Occupation Graph filename, remove "-" from front
filename2 = sprintf("%5.3GEqui_DensityZoom.png", abs(Dopant)); %Equi-Electron Density Graph filename, remove "-" from front

%Convert field from V/cm to V/A
Field = Field*(1E-8);

Data = calcbands(Dopant, Field, Pot); %Run simulation
%filename = ''
%saveas(gcf, '')



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



%Plot subbands occupancy
%Spread the matrix data into vectors so they may be plotted easily
scatter(Data.Vbot(:), Data.Vtop(:), [], subband_occ(:), 'filled')
xlabel("V_b_o_t")
ylabel("V_t_o_p")
D = sprintf('%5.3G', Dopant);  %Express doping conc. in scientific notation
                      
title("Delta Doping: " + D)

saveas(gcf, filename1)




%Calculate total carrier concentration inside the QW
tot_carrier_conc = Data.Sub1.Conc + Data.Sub2.Conc + Data.Sub3.Conc;


%Plot a color map of electron density in QW as a function of gate voltage
%TO DO: Add equi-density lines to color plot (Look up 2d contor plot)
figure

contourf(Data.Vbot, Data.Vtop, tot_carrier_conc, "ShowText", true, "LabelFormat", "%0.3G")
xlabel("V_b_o_t")
ylabel("V_t_o_p")
title("Delta Doping: " + D)

saveas(gcf, filename2)


%TO DO: Overlay contour over subband occupancy
%Fix x and y axis of all plots to easily overlap each other
%Look into "Hold on" function