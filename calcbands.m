function Data = calcbands(DeltaDopF, DeltaDopB, BackGate, FrontGate)

global aquila_control
%{
Calcbands and MAIN written by Christopher Cravey
christopheracravey@gmail.com
Northwestern University, Prof. Matthew Grayson lab 2023
All other code original to Aquila Subband Plotter
written by Dr. Martin Rother,  martin.rother@web.de
%}

% Calculate which subbands are occupied for variety of boundary conditions
% Quantum well structure is hard coded into this function, must manually
% change any structure properties by altering this code, but the passable
% arguments are as follows:
%
%ARGUMENTS:
%   DeltaDop is the doping you would like to have on either side, in cm^-2
%      note that the resulting doping will be X2, as this doping conc. is
%      applied to both sides.
%      DeltaDop < 0 is n-type doping
%   BackGate is the desired potential bias (V) that the doped back gate
%       (right side of conduction band graph) will be biased by. The fermi 
%       energy in the back gate region will be manually offset by each value.
%   FrontGate is the desired potential bias (V) that the metal front gate
%       will be biased by (left side of the device, the actual metal gate is
%       not simulated, but the resulting shift in the conduction band at the
%       device's edge is). The value of the conduction band with respect to the
%       fermi energy will be manually offset by each value
%      
%
%
%RETURNS: Data
%   Data.Bound_CondBG(Data.Bound_CondFG): Records the inputed back and
%           front gate voltages respectively
%   Data.WellConc: Carrier sheet density inside the QW in cm^-2
%   Data.Subx.Occ: [i,j,k] matrix. If 0, the x band is NOT occupied under the
%           boundary conditions BackGate(i), FrontGate(j), and thermal
%           broadening energy offset (k). Thermal broadening offset is user
%           defined via the "Broadening" array. The case of zero thermal
%           broadening offset is performed automatically, and is saved in
%           Data.Subx.Occ(:, :, 1)
%   Data.Subx.Conc: [i,j] matrix. The carrier concentration occupying the
%           x subband at BackGate(i) and FrontGate(j)
%   Data.Vtop / Data.Vbot: [i,j] matrix. Calculated top and bottom gate
%           potentials for BackGate(i) and FrontGate(j)
%   Data.Etop / Data.Ebot: [i,j] matrix. Calculated E field next to the top
%           and bottom gate for BackGate(i) and FrontGate(j)
%   Data.Breakdown: [i,j] matrix. Records the maximum ratio of the E field at the
%           edges of the device to the breakdown E field of GaAs
%           (4E5V/cm). A value > 1 implies the E field in the device
%           exceeds the breakdown field.



% Convert DeltaDop (cm^-2 to sheet doping cm^-3)
DeltaDopF = (DeltaDopF/100)/(2E-10);  %/100 to convert to m, divide by 2 A 
                                    % because the doping layer is 2A wide

disp("Front side doping conc. is: "+ DeltaDopF + " cm^-3")

DeltaDopB = (DeltaDopB/100)/(2E-10);  %/100 to convert to m, divide by 2 A 
                                    % because the doping layer is 2A wide

disp("Back side doping conc. is: "+ DeltaDopB + " cm^-3")



%How large are the arrays passed throdugh for BackGate and FrontGate arguments?
l_BG = length(BackGate);
l_FG = length(FrontGate);


%Keep record of the Front gate and Back Gate conditions passed
Data.Bound_CondBG = BackGate;
Data.Bound_condFG = FrontGate;  

for i = 1:l_BG

    for j = 1:l_FG

        initaquila                             %initialize
        aquila_control.mode=1;                 %1D-Simulation
        aquila_control.fix_doping=0;           %use doping levels instead of treating the doping
                                               %as a space charge




        %%%%%Prof. Grayson Degeneracy Cooler Structure%%%%%%% 
        %%%%%Similar to device HS1 https://doi.org/10.1063/1.4945090%%%%

%{    
%%%%   Symmetric QW Saturation for various barrier widths    %%%%%
        
(((((0.326% Al)))))
        60nm Front Barrier:
        54nm Back Barrier:
            Front Side: -5E11 
            Back Side:  -3E11
            Well Conc:  -5.77E11
        80nm Front Barrier:
        73nm Back Barrier:
            Front Side: -5.0E11 
            Back Side:  -2.2E11
            Well Conc:  -4.46E11
        100nm Front Barrier:
        92nm Back Barrier:
            Front Side: -5E11 
            Back Side:  -3E11
            Well Conc:  -3.64E11
        120nm Front Barrier:
        110nm Back Barrier:
            Front Side: ?E11 
            Back Side:  ?E11
            Well Conc:  -3.07E11
        

(((((0.25% Al)))))
        60nm Front Barrier:
        54nm Back Barrier:
            Front Side: -5E11 
            Back Side:  -2.7E11
            Well Conc:  -4.64E11
        80nm Front Barrier:
        73nm Back Barrier:
            Front Side: -5.5E11 
            Back Side:  -2.5E11
            Well Conc:  -3.60E11
        100nm Front Barrier:
        92nm Back Barrier:
            Front Side: -4.5E11 
            Back Side:  -2.2E11
            Well Conc:  -2.93E11
        120nm Front Barrier:
        110nm Back Barrier:
            Front Side: -4.3E11 
            Back Side:  -1.9E11
            Well Conc:  -2.49E11
        140nm Front Barrier:
        130nm Back Barrier:
            Front Side: -4.1E11 
            Back Side:  -1.3E11
            Well Conc:  -2.14E11
        
        
%}

        FBarr = 800;  %Front Barrier Width (A)
        BBarr = 800; %Back Barrier Width (A)
        
        add_mbox(1000,20,0,0);                  %1000 A GaAs Cap (surface)
        add_mbox(1350,50,0.25,0);              %1350 A AlGaAs
        add_mbox(2,1,0.25,0);                  %2 A AlGaAs to increase grid resolution
        add_mbox(2,1,0.25,DeltaDopF);           %2 A delta-doped AlGaAs
        add_mbox(2,1,0.25,0);                  %2 A AlGaAs
        add_mbox(FBarr,20,0.25,0);               %800 A AlGaAs spacer
        
        QW = 750;

        add_mbox(QW,5,0,0);                    %650 A GaAs quantum well
        
        add_mbox(BBarr,20,0.25,0);               %800 A AlGaAs spacer
        add_mbox(2,1,0.25,0);                  %2 A AlGaAs
        add_mbox(2,1,0.25, DeltaDopB);          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.25,0);                  %2 A AlGaAs to increase grid resolution
        add_mbox(9700,100,0.25,0);             %9700 A AlGaAs
        add_mbox(300, 50, 0, 0);                %300A GaAs cap 
        add_mbox(800, 10, 0.0, -7.5E17);        %bulk doped back gate
        add_mbox(300, 50, 0, 0);                %300A GaAs cap 

        
        
        add_bias([500, 2100 + FBarr], 3);                    %Make fermi energy between top gate and QW pinned to mid gap
        add_bias([2500 + FBarr + QW, 14800], 3);                  %Make fermi energy between bottom gate and QW pinned to mid gap
        add_qbox([2250+FBarr, 2350+FBarr + QW],5,3, GE +XE + LE);      %set quantum box onto quantum well
        add_pbox([1800, 4000 + FBarr + BBarr],CB);                    %Graph charge density in well
        add_pbox([0 15700],CB);                      %Graph charge density throughout structure


        add_bias([0,500], FrontGate(j));                             %Set front gate potential
        add_bias([12200+QW+FBarr+BBarr, 16000], BackGate(i));       %Set bottom gate potential
        add_boundary(LEFT, POTENTIAL, 0);            %Set Potential at surface = 0 (Bias is what creates front gate voltage)
        add_boundary(RIGHT, FIELD, 0);               %Set E field at the bulk of the device = 0

%3.6294e+11
%{
        %%%%% ETH Bottom Gate Experimental Comparison Structure %%%%%%%%
        %%%%%Based on device HS1 https://doi.org/10.1063/1.4945090%%%%
        add_mbox(100,50,0.0,0);                %100 A GaAs (surface)
        add_mbox(2000,10,0.328,0);             %2000 A AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
        add_mbox(2,1,0.328,DeltaDop);          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs
        add_mbox(800,20,0.328,0);              %800 A AlGaAs spacer
        
        add_mbox(270, 5, 0, 0);                 %270 A GaAs quantum well

        add_mbox(800,20,0.328,0);              %800 A AlGaAs spacer
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs
        add_mbox(2,1,0.328, DeltaDop/2.6);          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
        add_mbox(13000,100,0.328,0);            %13000 A AlGaAs
        add_mbox(300, 50, 0, 0);                %300A GaAs cap (for ETH comparison)
        add_mbox(800, 10, 0.0, -7.5E17);        %bulk doped back gate
        add_mbox(300, 50, 0, 0);                %300A GaAs cap (for ETH comparison)
        

        add_bias([300, 2800], 4);               %Make fermi energy between top gate and QW pinned to mid gap
        add_bias([3200 , 16600], 4);            %Make fermi energy between QW and bottom gate pinned to mid gap
        add_bias([16600 40000], BackGate(i));       %Add desired voltage to bottom gate
        add_qbox([2900 3180],5,3,GE + XE + LE);          %set quantum box onto quantum well
        add_pbox([2900 3180],CB);              %Graph charge density in well
        add_pbox([0 20000],CB);                %Graph charge density throughout structure


        add_boundary(LEFT,POTENTIAL, FrontGate(j));  %Set top gate potential
        add_boundary(RIGHT, FIELD, 0);               %Set E field at the bulk of the device = 0
%}

%{
        %%%%%%% ETH TOP GATE COMPARISON STRUCTURE D170301B  %%%%%%%%%

        add_mbox(100,50,0.0,0);                %100 A GaAs (Surface)
        add_mbox(400,10,0.328,0);              %400 A AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
        add_mbox(2,1,0.328,DeltaDop);          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs
        add_mbox(1000,20,0.328,0);              %1000 A AlGaAs spacer
        
        add_mbox(270, 5, 0, 0);                %270 A GaAs quantum well

        add_mbox(1000,20,0.328,0);             %1000 A AlGaAs spacer
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs
        add_mbox(2,1,0.328,(DeltaDop*(18/110)));          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
        add_mbox(10000,100,0.328,0);           %10000 A AlGaAs
         
        
        add_qbox([1490 1790], 5, 3,GE + XE + LE);          %set quantum box onto quantum well
        add_pbox([1490 1790],CB);              %Graph charge density in well
        add_pbox([0 13000],CB);                %Graph charge density throughout structure


        add_boundary(LEFT,POTENTIAL, FrontGate(j));  %Set top gate potential
        add_boundary(RIGHT, FIELD, 0);               %Set E field at the bulk of the device = 0
%}

        




        disp("FrontGate: " + FrontGate(j) + "   BackGate: " + BackGate(i));
        runstructure;                          %%%%Perform Simulation%%%%%
        


        
       
        
        
        
        %Introduce subband occupation variables. Initialize as 0
        Data.Sub1.Occ(i,j) = 0;
        Data.Sub2.Occ(i,j) = 0;
        Data.Sub3.Occ(i,j) = 0;
        
        



        %OPTIONAL: Thermal Broadening: Calculate subband occupation if the 
        %fermi energy is offset via a specific amount. This predicts which 
        %subbands a carrier would occupy if its energy was shifted from the 
        %fermi energy via thermal broadening.

        %TO SKIP THERMAL BROADENING CALCULATIONS: Set "Broadening" = []; 

        Kb = 8.61733E-5;   %Boltzmann Constant in eV/K
        T = 4;  %The termperature used to calculate THERMAL BROADENING ONLY
                %to change the environmental temperature, see initaquila.
        Broadening = [3, -3];  %Array: For each element "a", the simulation will
                          %add a*Kb*T to the fermi level then calculate
                          %which subbands are occupied. 
                           


        %Calculate subband occupation
        %Start with no thermal broadening considerations (Broad = 0), 
        %then iterate through user defined Broadening
        Broad = [0, Broadening]; 
        for k = 1:length(Broad)

            if (aquila_subbands.ge(1).E(1) - ((Broad(k)*Kb*T) + aquila_control.Efermi)) < 0
                Data.Sub1.Occ(i,j,k) = 1; 
            else, Data.Sub1.Occ(i,j,k) = 0;
            end
            if (aquila_subbands.ge(1).E(2) - ((Broad(k)*Kb*T) + aquila_control.Efermi)) < 0
                Data.Sub2.Occ(i,j,k) = 1;
            else, Data.Sub2.Occ(i,j,k) = 0;
            end
            if (aquila_subbands.ge(1).E(3) - ((Broad(k)*Kb*T) + aquila_control.Efermi)) < 0
                Data.Sub3.Occ(i,j,k) = 1; 
            else, Data.Sub3.Occ(i,j,k) = 0;
            end
        end

        % How much charge is in the first and in the second subband
        chG=genqcharge(1,GE,1); %charge distribution for gamma electrons in the first subband
                               %of the first quantum box
        chX=genqcharge(1,XE,1); %charge distribution for X electrons in the first subband
                               %of the first quantum box
        chL=genqcharge(1,LE,1); %charge distribution for L electrons in the first subband
                               %of the first quantum box                       
        Data.Sub1.Conc(i,j) = (sum(chG.*aquila_subbands.structure(1).boxvol) + sum(chX.*aquila_subbands.structure(1).boxvol) + sum(chL.*aquila_subbands.structure(1).boxvol))*1e16; 
        %integrate to get the sheet density
        

        chG=genqcharge(1,GE,2); %charge distribution for gamma electrons in the second subband
                               %of the first quantum box
        chX=genqcharge(1,XE,2); %charge distribution for X electrons in the second subband
                               %of the first quantum box
        chL=genqcharge(1,LE,2); %charge distribution for L electrons in the second subband
                               %of the first quantum box                       
        Data.Sub2.Conc(i,j) = (sum(chG.*aquila_subbands.structure(1).boxvol) + sum(chX.*aquila_subbands.structure(1).boxvol) + sum(chL.*aquila_subbands.structure(1).boxvol))*1e16; 
        %integrate to get the sheet density


        chG=genqcharge(1,GE,3); %charge distribution for gamma electrons in the third subband
                               %of the first quantum box
        chX=genqcharge(1,XE,3); %charge distribution for X electrons in the third subband
                               %of the first quantum box
        chL=genqcharge(1,LE,3); %charge distribution for L electrons in the third subband
                               %of the first quantum box                       
        Data.Sub3.Conc(i,j) = (sum(chG.*aquila_subbands.structure(1).boxvol) + sum(chX.*aquila_subbands.structure(1).boxvol) + sum(chL.*aquila_subbands.structure(1).boxvol))*1e16; 
        %integrate to get the sheet density

        
        %Calculate total carrier density in QW
        Data.WellConc(i,j) = Data.Sub1.Conc(i,j) + Data.Sub2.Conc(i,j) + Data.Sub3.Conc(i,j);


        % What is the potential of the top and bottom gate?
        Data.Vbot(i,j) = (aquila_control.Efermi - aquila_material.ec(end) + phi(end));
        Data.Vtop(i,j) = (phi(1) - aquila_material.ec(1));


        % What is the effective E field at the edges of our material?
        % Note that breakdown E field in GaAs is 4E5 V/cm
        
        Data.Etop(i,j) = ((aquila_material.ec(2) - phi(2) - ...
            aquila_material.ec(1) + phi(1))/(aquila_structure.xpos(2) ...
            - aquila_structure.xpos(1)))*(-1.0e8);
                      % Mulitply by 1.0e8 to convert from V/A to V/cm
        Data.Ebot(i,j) = ((aquila_material.ec(end-100) - phi(end-100) - ...
            aquila_material.ec(end-101) + phi(end-101))/(aquila_structure.xpos(end - 100) ...
            - aquila_structure.xpos(end - 101)))*(-1.0e8);
                      % Mulitply by 1.0e8 to convert from V/A to V/cm
        


        %Record how close the E field at edges of the device get to the
        %breakdown E field for GaAs
        Data.Breakdown (i,j) = max([abs(Data.Etop(i,j)/4e5) , abs(Data.Ebot(i,j)/4e5)]);
        %Records the % of breakdown E field reached for run i,j

       

    end
end


%%%%%%Plot the Wavefunctions of the bottom 2 bands%%%%%
%We only have gamma electrons in one QBOX, this is, what the ge(1) stand for.
%We have to square the wavefunctions to get the probability distribution.
l=length(aquila_subbands.structure(1).xpos); %the size of one wavefunction
plot(aquila_subbands.structure(1).xpos,[aquila_subbands.ge(1).psi(1:l)' ...
      aquila_subbands.ge(1).psi(l+1:2*l)'].^2)








%%%%Plot subband energies, the fermi energy, and thermal broadening%%%%%%

%{
%Find the slice of xpos array which corresponds to our quantum well qbox
fi = find(aquila_structure.xpos == aquila_subbands.structure(1).xpos(1));
li = find(aquila_structure.xpos == aquila_subbands.structure(1).xpos(end));
%Calculate conduction band in our quantum well
cb= aquila_material.ec(fi:li)-phi(fi:li);



figure
plot(aquila_subbands.structure(1).xpos, cb); %Plot Conduction Band
hold on
w = length(aquila_subbands.structure(1).xpos);
ylim([aquila_subbands.ge.E(1)- 0.005, aquila_subbands.ge.E(3) + 0.001]); %Set y-axis limits to be more zoomed in on subbands
plot(aquila_subbands.structure(1).xpos, ones(w,1)*aquila_control.Efermi); %Graph fermi energy
%Graph Subband Energies
plot(aquila_subbands.structure(1).xpos, ones(w,1)*aquila_subbands.ge.E(1), 'Color', 'g', 'LineWidth', 1);
plot(aquila_subbands.structure(1).xpos, ones(w,1)*aquila_subbands.ge.E(2), 'Color', 'g', 'LineWidth', 1);
plot(aquila_subbands.structure(1).xpos, ones(w,1)*aquila_subbands.ge.E(3), 'Color', 'g', 'LineWidth', 1);
%Graph vertical black bar showing fermi energy +/- 3KbT
plot([aquila_subbands.structure(1).xpos(end-25),aquila_subbands.structure(1).xpos(end-25)], ...
    [aquila_control.Efermi - (3*Kb*T), aquila_control.Efermi + (3*Kb*T)], 'LineWidth', 2, 'color', 'black')

xlabel("Distance (A)");
ylabel("Energy (eV)");
title("T = " + T  + "K   V_t_o_p = " + FrontGate + "V   V_b_o_t = " + BackGate + "V");
hold off
%}
end