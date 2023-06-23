%Example for AQUILA 1.0
%1D-Simulation of a quantum wire substrate

%Copyright 1999 Martin Rother
%
%This file is part of AQUILA.
%
%AQUILA is free software; you can redistribute it and/or modify
%it under the terms of the BSD License as published by
%the Open Source Initiative according to the License Policy
%on MATLAB(R)CENTRAL.
%
%AQUILA is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%BSD License for more details.

initaquila                             %initialize
aquila_control.mode=1;                 %1D-Simulation
aquila_control.fix_doping=0;           %use doping levels instead of treating the doping
                                       %as a space charge



add_mbox(1020,20,0,0);                 %1020 A GaAs Cap
add_mbox(1350,50,0.328,0);             %1350 A AlGaAs
add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
add_mbox(2,1,0.328,-.8e19);           %2 A delta-doped AlGaAs
add_mbox(2,1,0.328,0);                 %2 A AlGaAs
add_mbox(800,20,0.328,0);              %800 A AlGaAs spacer
add_mbox(650,5,0,0);                   %650 A GaAs quantum well
add_mbox(800,20,0.328,0);              %800 A AlGaAs spacer
add_mbox(2,1,0.328,0);                 %2 A AlGaAs
add_mbox(2,1,0.328,-0.8e19);           %2 A delta-doped AlGaAs
add_mbox(2,1,0.328,0);                 %2 A AlGaAs to increase grid resolution
add_mbox(10000,100,0.328,0);           %10000 A AlGaAs

add_qbox([3175 3850],5,3,GE);          %set quantum box onto quantum well


add_boundary(LEFT,POTENTIAL, 0);        %free surface, no gate
%add_boundary(RIGHT, POTENTIAL, 1)
add_boundary(RIGHT,FIELD, 0);           %transition to bulk

add_pbox([3150 3850],CB);
add_pbox([0 14600],CB);

runstructure;                          %DO IT

%how do the wavefunctions look like
%We only have gamma electrons in one QBOX, this is, what the ge(1) stand for.
%We have to square the wavefunctions to get the probability distribution.
l=length(aquila_subbands.structure(1).xpos); %the size of one wavefunction
plot(aquila_subbands.structure(1).xpos,[aquila_subbands.ge(1).psi(1:l)' ...
      aquila_subbands.ge(1).psi(l+1:2*l)'].^2)


%what is the sheet density of the electrons in the channel
ch=sumcharge(phi);  %the total charge
ix=boxindex(aquila_structure.xpos,[3150 3850]); %the channel position
%sum the charge weighted by the area covered by each node and scale it
%from A^-2 to cm^-2. The output is negative because of the negative electron charge.
sum(ch(ix).*aquila_structure.boxvol(ix))*1e16  %the electron sheet density in the channel



%now compare the energy levels to the Fermi energy
%aquila_subbands.ge(1).E(1) - aquila_subbands.xe(1).E(1)
%aquila_subbands.le(1).E(1) - aquila_subbands.xe(1).E(1)


aquila_subbands.ge(1).E(1) - aquila_control.Efermi
aquila_subbands.ge(1).E(2) - aquila_control.Efermi
aquila_subbands.ge(1).E(3) - aquila_control.Efermi

%how much charge is in the first and in the second subband
ch=genqcharge(1,GE,1); %charge distribution for gamma electrons in the first subband
                       %of the first quantum box
sum(ch.*aquila_subbands.structure(1).boxvol)*1e16 %integrate to get the sheet density
ch=genqcharge(1,GE,2); %charge distribution for gamma electrons in the second subband
                       %of the first quantum box
sum(ch.*aquila_subbands.structure(1).boxvol)*1e16 %integrate to get the sheet density


