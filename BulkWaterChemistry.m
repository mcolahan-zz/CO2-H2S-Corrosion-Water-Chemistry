
    %WaterChemistry

function OutputArray = BulkWaterChemistry(Temperature, pCO2, pH2S, ConcArray)

    %Temperature (°C)
    %pCO2, pH2S are in Bar
    %ConcArray: 1 x nSpec row array with the species concencentrations (mol/L)
    %           which serves as the initial guess of species concentrations
    %           Mainly serves as the input of spectator ions (not those
    %           associated with CO2 and H2S. 
    %Species Numbers which coorespond to array indices
    nCO2 = 1;
    nH2CO3 = 2;
    nHCO3 = 3;
    nCO3 = 4;
    nH2S = 5;
    nHS = 6;
    nS = 7;
    nOH = 8;
    nH = 9;
    nFe = 10;
    nNa = 11;
    nCl = 12;
    
    nSpec = 12;

    % Charges of the species 1:nSpec
    ChargeArray = [0,0,-1,-2, 0, -1, -2, -1, 1, 2, 1, -1];

%{
    MW = zeros(1,nSpec); %g/mol
    MW(nCO2) = 44;
    MW(nH2CO3) = 62.03;
    MW(nHCO3) = 61.02;
    MW(nCO3) = 60.01;
    MW(nH2S) = 34.08;
    MW(nHS) = 33.07;
    MW(nS) = 32.06;
    MW(nOH) = 17.01;
    MW(nH) = 1.01;
    MW(nFe) = 55.845;
    MW(nNa) = 23;
    MW(nCl) = 35.45;
%}
    
    %Initialize:
    %Nodes which indicate which row of ConcArray stand for what location in
    %the point model
    nbulk = 1;
    
    
    xleft = 0;
    xright = 14;
    Ionic = 0;
    
    % Bisection Method to solve for pH
    while abs(xright - xleft) > 0.01
        
        ChargeBalance1 = EquilibriumCalc(xleft);
        root = (xleft + xright)/2;
        ChargeBalance2 = EquilibriumCalc(root);
        
        if (ChargeBalance1 * ChargeBalance2 < 0)
            xright = root;
        else
            xleft = root;
        end
        
    end
    
    OutputArray = ConcArray;


function EquilibriumCalc = EquilibriumCalc(pH)
    
    %Nested function within BulkWaterChemistry that calculates how far away
    %from neutral charge a solution is with a given pH. The Bisection loop
    %in BulkWaterChemistry inputs a guess pH, then this function and
    %subsequent functions calculate concentration values based on this pH
    %and then the sum of charges (conc * charge) is computed. The sum of
    %charges should converge to zero as a pH is narrowed in upon.
    
    ConcArray(nbulk, nH) = 10 ^ - pH;
    
    Ionic = IonicStrength();
    Concentrations;
    
    ChargeBalanceSum = 0;
    for i = 1:nSpec
        ChargeBalanceSum = ChargeBalanceSum + ConcArray(i) * ChargeArray(i);
    end
    
    EquilibriumCalc = ChargeBalanceSum;



function Ionic = IonicStrength()

    Ionic = 0;
    for j = 1:nSpec
        Ionic = Ionic + ConcArray(j) * (ChargeArray(j)^2);
    end
    Ionic = 0.5 * Ionic;
    
    
end


function Concentrations()

    EquilibriumConstants = ReactionConstants(Temperature,Ionic,pCO2);
    
    %ReactionConstants = [Kwa_H2O, kf_wa, kb_wa; K_CO2_sol, 0, 0; ...
    %  K_CO2_hy, kf_CO2_hy, kb_CO2_hy; K_CO2_ca, kf_ca, kb_ca; K_CO2_bi, kf_bi, ...
    %  kb_bi; K_H2S_sol, 0, 0; K_H2S_hs, kf_H2S, kb_H2S; K_H2S_hsa, kf_HS, kb_HS];
    
    %concentration of OH- in water, mol/L
    ConcArray(nbulk, nOH) = EquilibriumConstants(1,1) / ConcArray(nbulk, nH);
    
    %concentration of CO2 in water, mol/L
    ConcArray(nbulk, nCO2) = pCO2 * EquilibriumConstants(2,1);
    
    %concentration of H2CO3 in water,mol/L
    ConcArray(nbulk, nH2CO3) = EquilibriumConstants(3,1) * ConcArray(nbulk, nCO2);
    
    %concentration of HCO3- in water, mol/L
    ConcArray(nbulk, nHCO3) = EquilibriumConstants(4,1) * ConcArray(nbulk, nH2CO3) / ConcArray(nbulk, nH);
    
    %concentration of CO3= in water, mol/L
    ConcArray(nbulk, nCO3) = EquilibriumConstants(5,1) * ConcArray(nbulk, nHCO3) / ConcArray(nbulk, nH);
    
    % concentration of H2S in bulk solution, mol/L
    ConcArray(nbulk, nH2S) = pH2S * EquilibriumConstants(6,1);
    
    % concentratio of HS-, mol/L
    ConcArray(nbulk, nHS) = EquilibriumConstants(7,1) * ConcArray(nbulk, nH2S) / ConcArray(nbulk, nH);
    
    %concnetration of S=, kmol/m3
    ConcArray(nbulk, nS) = EquilibriumConstants(8,1) * ConcArray(nbulk, nHS) / ConcArray(nbulk, nH);
    

end


     
end
    
end

    

function ChemicalReactionConstants = ReactionConstants(t,Ionic,pCO2)
   
    %Equilibrium constant for H2O dissociation
    Kwa_H2O = 10 ^ (-(29.3868 - 0.0737549 * (t + 273.15) + 7.47881 * 10 ^ (-5)...
        * (t + 273.15) ^ 2));
    %Y.K. Kharaka, E.H. Perkins, W.D. Gunter, J.D. Debral, C.H.Bamford, ¡°Solmineq 88: A Computer Program for Geochemical Modeling of Water Rock Interactions¡± (Menlo Park, CA: Alberta Research Council, 1989).
    %forward reaction rate constant for H2O dissociation
    kb_wa = 7.85 * 10 ^ 10;
    %backward reaction constant for H2O dissociation
    kf_wa = Kwa_H2O * kb_wa;
    
    %Solubility constant for gaseous CO2 dissolution
    K_CO2_sol = 14.5 * 10 ^ (-(2.27 + 5.65 * 10 ^ (-3) * (t * 9 / 5 + 32) - 8.06 * 10 ^ (-6) ...
        * (t * 9 / 5 + 32) ^ 2 + 0.075 * Ionic)) / 1.00258;
    %J.E. Oddo, M.B. Tomson, ¡°Simplified Calculation of CaCO3 Saturation at High Temperatures and Pressures in Brine Solutions,¡±SPE of AIME (Richardson, TX: Society of Petroleum Engineers,1982), p. 1,583.
    %Temperature: 50-200 Celsius; Pressure: 1-1200 bar; Ionic strength: 0-1
    
    %Equilibrium constant for CO2 hydration: CO2+H2O==>H2CO3
    K_CO2_hy = 2.58 * 10 ^ (-3);
    %D.A. Palmer, R. van Eldik, Chem. Rev. 83 (1983): p. 651.
    %Temperature: 25 Celsius, K_CO2_hy = 2.31*10^(-3) at 300 Celsius, varies little with temperature
    %Forward reaction constant for CO2 hydration
    kf_CO2_hy = 10 ^ (329.85 - 110.541 * log10(t + 273.15) -...
        (17265.4 / (t + 273.15)));
    %Backward reaction constant for Co2 hydration
    kb_CO2_hy = kf_CO2_hy / K_CO2_hy;
    
    %Reaction constants for H2CO3 disociation: H2CO3==>HCO3-+H+
    K_CO2_ca = 387.6 * 10 ^ (-(6.41 - 1.594 * 10 ^ (-3) * (t * 9 / 5 + 32) + ...
        8.52 * 10 ^ (-6) * (t * 9 / 5 + 32) ^ 2 - 3.07 * 10 ^ (-5) * 14.5 - ...
        0.4772 * Ionic ^ 0.5 + 0.118 * Ionic));
    %J.E. Oddo, M.B. Tomson, ¡°Simplified Calculation of CaCO3 Saturation at High Temperatures and Pressures in Brine Solutions,¡±SPE of AIME (Richardson, TX: Society of Petroleum Engineers,1982), p. 1,583.
    %Temperature: 50-200 Celsius; Pressure: 1-1200 bar; Ionic strength: 0-1
    kf_ca = 10 ^ (5.71 + 0.0526 * t - 0.000294 * t ^ 2 + 0.000000791 * t ^ 3);
    kb_ca = kf_ca / K_CO2_ca;
    
    %Reaction constants for HCO3- disociation: HCO3- ==> CO32- + H+
    K_CO2_bi = 10 ^ (-(10.61 - 4.97 * 10 ^ (-3) * (t * 9 / 5 + 32) +...
        1.331 * 10 ^ (-5) * (t * 9 / 5 + 32) ^ 2 - 2.624 * 10 ^ (-5) * pCO2...
        * 14.5 - 1.166 * Ionic ^ 0.5 + 0.3466 * Ionic));
    %J.E. Oddo, M.B. Tomson, ¡°Simplified Calculation of CaCO3 Saturation at 
    %High Temperatures and Pressures in Brine Solutions,¡±SPE of AIME
    %(Richardson, TX: Society of Petroleum Engineers,1982), p. 1,583.
    %Temperature: 50-200 Celsius; Pressure: 1-1200 bar; Ionic strength: 0-1
    kf_bi = 1000000000;
    kb_bi = kf_bi / K_CO2_bi;
    
    %Gaseous H2S dissolution, H2S(g)<===> H2S(aq), use equations from Suleimenov (1994)
    K_H2S_sol = 10 ^ (-(634.27 + 0.2709 * (t + 273.15) - 0.00011132 * (t + 273.15)...
        ^ 2 - 16719 / (t + 273.15) - 261.9 * log(t + 273.15) / log(10)));

    %H2S first disociation, H2S<-->HS- + H+, use equations from Suleimenov (1994)
     K_H2S_hs = 10 ^ (782.43945 + 0.36126 * (t + 273.15) - 0.00016722 * (t + 273.15) ^ 2 ...
     - 20565.7315 / (t + 273.15) - 142.7417222 * log(t + 273.15));
     kf_H2S = 10000;
     kb_H2S = kf_H2S / K_H2S_hs;
     %calculated by the equation from Kharaka (1988)
    %K_H2S_hs = 10 ^ (-15.345 + 0.045676 * (t + 273.15) - 5.9666 * 10 ^ (-5) * (t + 273.15) ^ 2)
    %Y.K. Kharaka, E.H. Perkins, W.D. Gunter, J.D. Debral, C.H.Bamford, ¡°Solmineq 88: A Computer Program for Geochemical Modeling of Water Rock Interactions¡± (Menlo Park, CA: Alberta Research Council, 1989)
    %Temperature: 0-250 Celsius; Pressure: total pressure equal to the vapor pressure of water or 1 bar
    
    %HS- dissociation, HS-<---->S2-+H+
    K_H2S_hsa = 10 ^ (-23.93 + 0.030446 * (t + 273.15) - 2.4831 * 10 ^ (-5) ...
        * (t + 273.15) ^ 2);
    %Y.K. Kharaka, E.H. Perkins, W.D. Gunter, J.D. Debral, C.H.Bamford, ¡°Solmineq 88: A Computer Program for Geochemical Modeling of Water Rock Interactions¡± (Menlo Park, CA: Alberta Research Council, 1989)
    %Temperature: 0-250 Celsius; Pressure: total pressure equal to the vapor pressure of water or 1 bar
    kf_HS = 1;
    kb_HS = kf_HS / K_H2S_hsa;
    
    % Chemical reaction constants are in rows by reaction with the
    % equalibrium constant column 1, forward reaction rate in column 2 and
    % the backward reaction rate in 3rd column
    ChemicalReactionConstants = [Kwa_H2O, kf_wa, kb_wa; K_CO2_sol, 0, 0; K_CO2_hy, ...
        kf_CO2_hy, kb_CO2_hy; K_CO2_ca, kf_ca, kb_ca; K_CO2_bi, kf_bi, kb_bi; ...
        K_H2S_sol, 0, 0; K_H2S_hs, kf_H2S, kb_H2S; K_H2S_hsa, kf_HS, kb_HS];
end
