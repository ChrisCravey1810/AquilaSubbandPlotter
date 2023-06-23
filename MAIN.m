clear all



Pot = [0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5] % in eV
Field = [0.0000003, 0.0000002, 0.0000001, 0, -0.0000001, -0.0000002, -0.0000003] % in 
Dopant = -2E11   % in cm^-2. This dopant conc is applied to both sides


Data = calcbands(Dopant, Field, Pot)

%Calculate which subbands are occupied
subband_occ = zeros([length(Field), length(Pot)])
for i = 1:length(Field)
    for j = 1:length(Pot)
        subband_occ(i,j) = (Data.Sub1.Occ(i,j) + Data.Sub2.Occ(i,j) + Data.Sub3.Occ(i,j))
    end
end

y 

scatter(Data.Vbot, Pot, [], subband_occ)

