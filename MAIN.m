clear all

%Please manually insert desired parameters
Pot = [0, -0.02, -0.05, -0.07, -0.1, -0.12, -0.15, -0.18, -0.2] % array in eV
%Field = [20000, 15000, 10000, 5000, 1000, 0, -100000, -200000] % array in V/cm
Field = [15000, 13000, 11000, 9000, 7000, 5000]
Dopant = -2E11   % Single value in cm^-2. This dopant conc is applied to both sides


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


%Overlay contour over subband occupancy
%Look into "Hold on" function