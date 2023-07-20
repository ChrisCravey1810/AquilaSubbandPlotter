function Data = calcbands(DeltaDop, Field, Pot)
% Calculate which subbands are occupied for varied of boundary conditions
% Quantum well structure is hard coded into this function, must manually
% change any structure properties by altering this code, but the passable
% arguments are as follows:
%
%ARGUMENTS:
%   DeltaDop is the doping you would like to have on either side, in cm^-2
%      note that the resulting doping will be X2, as this doping conc. is
%      applied to both sides.
%      DeltaDop < 0 is n-type doping
%   Field is the desired E field boundary condition at the bulk end of the
%      device (right side). This will eventually be converted into a 
%      bottom gate potential to simply graph the subband occupancy as a 
%      function of both gate voltages
%   Pot (TOP GATE VOLTAGE) is the desired Potential boundary condition at 
%      the surface of the device (left side)
%
%
%RETURNS: Data
%   Data.Bound_CondF(Data.Bound_CondP): The passed field or potential
%           boundary conditions repectively
%   Data.Well.Conc: Carrier sheet density inside the QW in cm^-2
%   Data.Subx.Occ: [i,j] matrix. If 0, the x band is NOT occupied under the
%           boundary conditions Field(i) and Pot(j)
%   Data.Subx.Conc: [i,j] matrix. The carrier concentration occupying the
%           x subband at Field(i) and Pot(j)
%   Data.Vtop / Data.Vbot: [i,j] matrix. Calculated top and bottom gate
%           potentials for Field(i) and Pot(j)
%   Data.Etop / Data.Ebot: [i,j] matrix. Calculated E field next to the top
%           and bottom gate for Field(i) and Pot(j)
%   Data.Breakdown: [i,j] matrix. Records how close the E field at the
%           edges of the device get to the breakdown E field of GaAs
%           (4e5V/cm)



% Convert DeltaDop (cm^-2 to sheet doping cm^-3)
DeltaDop = (DeltaDop/100)/(2E-10);  %/100 to convert to m, divide by 2 A 
                                    % because the doping layer is 2A wide
disp("Doping conc. on each side is: "+ DeltaDop + " cm^-3")



%Multiply desired potential by -1
%For some reason this program takes the opposite sign of potential passed
%Pot = Pot*-1;



%How large are the arrays passed through for Field and Potential arguments?
l_F = length(Field);
l_P = length(Pot);


%Keep record of the field and potential boundary conditions passed
Data.Bound_CondF = Field;
Data.Bound_condP = Pot;  

for i = 1:l_F

    for j = 1:l_P

        initaquila                             %initialize
        aquila_control.mode=1;                 %1D-Simulation
        aquila_control.fix_doping=0;           %use doping levels instead of treating the doping
                                               %as a space charge

        %%%Prof. Grayson Structure Degeneracy Cooler Structure%%% 

        
        add_mbox(1000,20,0,0);                 %1000 A GaAs Cap
        add_mbox(1350,50,0.328,0);             %1350 A AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
        add_mbox(2,1,0.328,DeltaDop);          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.328,0);                 %2 A AlGaAs
        add_mbox(800,20,0.328,0);              %800 A AlGaAs spacer
        
        add_mbox(650,5,0,0);                   %650 A GaAs quantum well

        add_mbox(800,20,0.328,0);               %800 A AlGaAs spacer
        add_mbox(2,1,0.328,0);                  %2 A AlGaAs
        add_mbox(2,1,0.328, DeltaDop);          %2 A delta-doped AlGaAs
        add_mbox(2,1,0.328,0);                  %2 A AlGaAs to increase grid resolution
        add_mbox(9700,100,0.328,0);             %9700 A AlGaAs
        add_mbox(300, 50, 0, 0);                %300A GaAs cap 
        add_mbox(800, 10, 0.0, -7.5E17);        %bulk doped back gate
        add_mbox(300, 50, 0, 0);                %300A GaAs cap 

        
        add_bias([500, 3000], 3);               %Make fermi energy between top gate and QW pinned to mid gap
        add_bias([4000, 14000], 3);             %Make fermi energy between bottom gate and QW pinned to mid gap
        add_bias([14000, 16000], Field(i));
        add_qbox([3100 3850],5,3, GE +XE + LE);          %set quantum box onto quantum well
        add_pbox([3100 3850],CB);              %Graph charge density in well
        add_pbox([0 15600],CB);                %Graph charge density throughout structure




        %%%%% ETH Bottom Gate Experimental Comparison Structure %%%%%%%%
%{        
        add_mbox(100,50,0.0,0);                %100 A GaAs
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
        add_mbox(9700,100,0.328,0);            %9700 A AlGaAs
        add_mbox(300, 50, 0, 0);                %300A GaAs cap (for ETH comparison)
        add_mbox(800, 10, 0.0, -7.5E17);        %bulk doped back gate
        add_mbox(300, 50, 0, 0);                %300A GaAs cap (for ETH comparison)
        

        add_bias([300, 2800], 4);               %Make fermi energy between top gate and QW pinned to mid gap
        add_bias([3200 , 13600], 4);            %Make fermi energy between QW and bottom gate pinned to mid gap
        add_bias([13600 40000], Field(i));       %Add desired voltage to bottom gate
        add_qbox([2900 3180],5,3,GE + XE + LE);          %set quantum box onto quantum well
        add_pbox([2900 3180],CB);              %Graph charge density in well
        add_pbox([0 20000],CB);                %Graph charge density throughout structure

%}


        %%%%%%% ETH TOP GATE COMPARISON STRUCTURE %%%%%%%%%
%{
        add_mbox(100,50,0.0,0);                %100 A GaAs
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
%}



        add_boundary(LEFT,POTENTIAL, Pot(j));  %Set top gate potential
        %add_boundary(RIGHT, POTENTIAL, 1)
        %add_boundary(RIGHT,FIELD, Field(i));   %Set E field at bottom gate (V/A)
        add_boundary(RIGHT, FIELD, 0);

        
        

        runstructure;                          %DO IT
        


        %Plot the Wavefunctions of the bottom 2 bands
        %We only have gamma electrons in one QBOX, this is, what the ge(1) stand for.
        %We have to square the wavefunctions to get the probability distribution.
        l=length(aquila_subbands.structure(1).xpos); %the size of one wavefunction
        plot(aquila_subbands.structure(1).xpos,[aquila_subbands.ge(1).psi(1:l)' ...
              aquila_subbands.ge(1).psi(l+1:2*l)'].^2)
       
        
        
        
        %Introduce subband variables. If 0, subband is unoccupied at 0T
        %If 1, subband is occupied at 0T
        
        
        Data.Sub1.Occ(i,j) = 0;
        Data.Sub2.Occ(i,j) = 0;
        Data.Sub3.Occ(i,j) = 0;
        
        
        %Calculate if subband energy is above fermi energy
        if (aquila_subbands.ge(1).E(1) - aquila_control.Efermi) < 0
            Data.Sub1.Occ(i,j) = 1; 
        end
        if (aquila_subbands.ge(1).E(2) - aquila_control.Efermi) < 0
            Data.Sub2.Occ(i,j) = 1;
        end
        if (aquila_subbands.ge(1).E(3) - aquila_control.Efermi) < 0
            Data.Sub3.Occ(i,j) = 1; 
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
        Data.Well.Conc(i,j) = Data.Sub1.Conc(i,j) + Data.Sub2.Conc(i,j) + Data.Sub3.Conc(i,j)


        % What is the potential of the top and bottom gate?
        %Data.Vbot(i,j) = (-0.71313 - aquila_material.ec(end) + phi(end));
        Data.Vbot(i,j) = (-0.7424 - aquila_material.ec(end) + phi(end));
        Data.Vtop(i,j) = (phi(1) - aquila_material.ec(1));


        % What is the effective E field at the edges of our material?
        % Note that breakdown E field in GaAs is 4E5 V/cm
        
        Data.Etop(i,j) = ((aquila_material.ec(2) - phi(2) - ...
            aquila_material.ec(1) + phi(1))/(aquila_structure.xpos(2) ...
            - aquila_structure.xpos(1)))*(-1.0e8);
                      % Mulitply by 1.0e8 to convert from V/A to V/cm
        Data.Ebot(i,j) = ((aquila_material.ec(end) - phi(end) - ...
            aquila_material.ec(end-1) + phi(end-1))/(aquila_structure.xpos(end) ...
            - aquila_structure.xpos(end - 1)))*(-1.0e8);
                      % Mulitply by 1.0e8 to convert from V/A to V/cm
        


        %Record how close the E field at edges of the device get to the
        %breakdown E field for GaAs
        Data.Breakdown (i,j) = max([abs(Data.Etop(i,j)/4e5) , abs(Data.Ebot(i,j)/4e5)]);
        %Records the % of breakdown E field reached for run i,j



        %{
        if abs(Data.Ebot(i,j)) > 4.0e5 || abs(Data.Etop(i,j)) > 4.0e5
            Data.Breakdown(i,j) = 1;
        else
            Data.Breakdown(i,j) = 0;
        end
        %}
        

    end
end
end