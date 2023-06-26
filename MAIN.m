clear all



Pot = [1, 0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8, -1] % in eV
Field = [0.0002, 0.0001, 0.00005, 0.00001, 0, -0.00001, -0.00005, -0.0001] % in 
Dopant = -2E11   % in cm^-2. This dopant conc is applied to both sides


Data = calcbands(Dopant, Field, Pot)

%Calculate which subbands are occupied
subband_occ = zeros([length(Field), length(Pot)])
for i = 1:length(Field)
    for j = 1:length(Pot)
        subband_occ(i,j) = (Data.Sub1.Occ(i,j) + Data.Sub2.Occ(i,j) + Data.Sub3.Occ(i,j))
    end
end

%Spread the matrix data into vectors so they may be plotted easily
Vbot = Data.Vbot(:);
Vtop = Data.Vtop(:);
Occupancy = subband_occ(:);

colormap cool; %Make the plot look nice:
%3 subbands occupied = ?
%2 subbands occupied = Pink 
%1 subbands occupied = Purple
%0 subbands occupied = Light Blue

scatter(Vbot, Vtop, [], occupancy, 'filled')
xlabel("V_b_o_t")
ylabel("V_t_o_p")
title("Delta Doping: " + Dopant)
