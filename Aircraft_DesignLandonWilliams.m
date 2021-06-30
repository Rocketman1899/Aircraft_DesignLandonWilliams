clc
clear all 
%% Specs

Pax = 225;              % Number of Passangers 
Seats_Abreast = 6;
W_Cargo =600;         % Weight of Cargo (lbs)
Range = 7400;
tofl = 10500/1000;       % Take off Field Length
AP = 140;               % Approach Speed (KTS)
Fuel_percentage = .350;  % Fuel left at landing
M_cruise = .85;          % Mach Number @ cruise
Alt =35000;             % Altitude at cruise
%Lambda = deg2rad(15);  % Sweep Angle
Ar = (6:.5:10);                 % Aspect Ratio
sigma = 0.953;
NE = 2;                 % Number of Engines
international = 0; 
SuperCrit = 0;
NAIS =1 

%% Chart Load
%  if SuperCrit == 0
%         deltaM_div = readtable('Fig.2_Conventional.xlsx');
%         deltaM_div = table2array(deltaM_div);
%     else
%         deltaM_div = readtable('Fig.2_SuperCritical.xlsx');
%         deltaM_div = table2array(deltaM_div);
%     end


    CL_maxLG = readtable('fig3_landing.csv');
    CL_maxLG = table2array(CL_maxLG); 

    CL_maxTo = readtable('fig3_takeoff.csv');
    CL_maxTo = table2array(CL_maxTo); 


    WeightRatio = readtable('fig.4.xlsx');
    WeightRatio = table2array(WeightRatio);


    if NE == 2
        TOFL = readtable('Fig.5_2Engines.xlsx');
        TOFL = table2array(TOFL);
    end


    if NE == 3
        TOFL = readtable('Fig.5_3Engines.csv');
        TOFL = table2array(TOFL);
    end



    if NE == 4
        TOFL = readtable('Fig.5_4Engines.csv');
        TOFL = table2array(TOFL);
    end

    AtmosphericVals = readtable('AtmosphericTable.xlsx');
    AtmosphericVals = table2array(AtmosphericVals);
    Speed_of_Sound  = AtmosphericVals(:,8); %Speed of Sound Values
    Altitude_Vals   =  AtmosphericVals(:,1); %Altitude Values
    Speed_of_Sound  = interp1(Altitude_Vals,Speed_of_Sound,30);   
    Kinematic_Visc = AtmosphericVals(30,10);

    % Cf Wing
    SkinFriction= readtable('Shevell_Cf.csv');
    SkinFriction= table2array(SkinFriction);

    %K_Fuselage ("Form Factor")
    FuselageFormFactor = readtable('Shevell_LDKf.csv');
    FuselageFormFactor = table2array(FuselageFormFactor);

    %Profile Drag 
    CD_pTakeoff = readtable('Fig6_Takeoff.xlsx');
    CD_pLanding = readtable('Fig6_Landing.xlsx');
    CD_pTakeoff = table2array(CD_pTakeoff);
    CD_pLanding = table2array(CD_pLanding);



format long

%% AR Loop
for ARcount = 1:length(Ar)
AR = Ar(ARcount);
lambda = deg2rad([10 15 20 25 30 35]);
for Lambdacount = 1:length(lambda)
%% Chart Loading 
Lambda = lambda(Lambdacount);
if Lambdacount == 1

    %CLmax Clean
    if rad2deg(Lambda) < 7.5 
    CLmaxClean = readtable('CLmaxClean_0.xlsx');
    CLmaxClean = table2array(CLmaxClean);
    end

    if rad2deg(Lambda) >= 7.5
    CLmaxClean = readtable('CLmaxClean_15.xlsx');
    CLmaxClean = table2array(CLmaxClean);
    end

    if rad2deg(Lambda) > 15
    CLmaxClean = readtable('CLmaxClean_15.xlsx');
    CLmaxClean = table2array(CLmaxClean);
    end
end  
%% THRUST REQUIRED LOOP 
TReq_JT9D_IC = 0;
WT_Adjustment = 0;
while TReq_JT9D_IC < 10000
%% Range Condition Loop
    count = 0;
    RangeCondition = 1;
    Fuel_Adjustment = 0;

    while RangeCondition == true 
%%  CL Convergence   
%%  CL Convergence   
        Cl_Condition = 1 ;
        cl = 0;
        clguess = .4500;    %Asumme a Value for cl 

         while  Cl_Condition == true
            
             % Use figure 2 to find deltaM_div
             
          
             if SuperCrit == 0
             DeltaM_div =-0.348*clguess + 0.191;
             else
             DeltaM_div =-0.0393 + 0.288*clguess + -0.393*clguess^2;
             end
            % Calculate M_div
             M_div = (M_cruise + 0.004) - DeltaM_div;
   
            % Use Figure 1 to calculate t/c
            
            %Figure 1 Lambda = 0
            if  Lambda == deg2rad(0)
                T_C = -0.635*M_div + 0.573;
            elseif Lambda == deg2rad(0) && SuperCrit == 1
                T_C = -38.239*M_div^3 + 89.121*M_div^2 - 69.895*M_div + 18.529;
            end
            
            % Figure 1 Lambda = 10
            if  Lambda == deg2rad(10)
                T_C =  -0.618*M_div + 0.565;
            elseif Lambda == deg2rad(10) && SuperCrit == 1
                T_C = -22.082*M_div^3 + 53.636*M_div^2 - 43.991*M_div + 12.253;
            end
            
            % Figure 1 Lambda = 15
            if  Lambda == deg2rad(15)
                T_C = -0.594*M_div + 0.552;
            elseif Lambda == deg2rad(15) && SuperCrit == 1
                T_C = -22.822*M_div^3 + 55.879*M_div^2 - 46.157*M_div + 12.938;
                
            end
            
            % Figure 1 Lambda = 20
            if Lambda == deg2rad(20)
                T_C =  -0.567*M_div + 0.538;
            elseif Lambda == deg2rad(20) && SuperCrit == 1
                T_C =  -14.76*M_div^3 + 37.802*M_div^2 - 32.747*M_div + 9.6624;
            end
            
            % Figure 1 Lambda = 25
            if Lambda == deg2rad(25)
                T_C =  -0.567*M_div + 0.538;
            elseif Lambda == deg2rad(25) && SuperCrit == 1
                T_C =  -19.337*M_div^3 + 49.806*M_div^2 - 43.268*M_div + 12.754;
            end
            
            % Figure 1 Lambda = 30
            if Lambda == deg2rad(30)
                T_C = -0.505*M_div + 0.506;
            elseif Lambda == deg2rad(30) && SuperCrit == 1
                T_C = -17.237*M_div^3 + 46.325*M_div^2 - 41.998*M_div + 12.913;
                
            end
            
            % Figure 1 Lambda = 35
            if Lambda == deg2rad(35)
                T_C = -0.469*M_div + 0.486;
            elseif Lambda == deg2rad(35) && SuperCrit == 1
                T_C = -40.421*M_div^3 + 108.03*M_div^2 - 96.966*M_div + 29.311;
            end
            
            % Figure 1 Lambda = 40
            if Lambda == deg2rad(40)
                T_C = -0.43*M_div + 0.466;
            elseif Lambda == deg2rad(40) && SuperCrit == 1
                T_C = -377.09*M_div^3 + 993.77*M_div^2 - 874.63*M_div + 257.2;
            end
            
            
            P = (cos(Lambda))^2*(T_C)^2*AR;
            % Figure 3 Cl_maxTo Cl_maxLG
            if P <.04 
              Cl_maxTo = 2.19 + 11*P + -23*P^2;
              Cl_maxLG = 1.18 + 12.9*P + -30.6*P^2;
            else
            Cl_maxTo = 0.558*log(P) + 3.440;     %Cl Max Take off
            Cl_maxLG =  0.5774*log(P) + 4.4069;
            end 
            
            %  Calculate Wing Loading
            WSLD = (AP/1.3)^2*((sigma*Cl_maxLG)/296);      
            % Calculate "All-out" Range From
            V_cruise = M_cruise*(576.4);
            RAO  = Range + 200 + .75*V_cruise;
            % Obtain an Initial Prediction for the fuel from figure 4
            W_RatioJT8D = 0.0234 + 1.03E-04*RAO + -5.48E-09*RAO^2;
            W_ratioJT9D = (W_RatioJT8D*(.61/.78)*1.05)+ Fuel_Adjustment; 
            % Wingloading Take off
            WL_To = WSLD/(1-((1-Fuel_percentage)*W_ratioJT9D)); 
            % Cruise Wingloading at Cruise
            WL_C = .965*WL_To;  
            cl =  WL_C /(1481*.2360*M_cruise^2);
            Cl_Condition = (clguess - cl) <= .001; 
            clguess = clguess+.0001;
         end 

%% Weight to Thrust | TOFL CALC
         W_T7Lo = interp1(TOFL(:,2),TOFL(:,1), tofl);
         WTTO = (W_T7Lo/WL_To)*(.953)*(Cl_maxTo);  %Weight to Thrust @ Take off
         VL_o = 1.2*(((296*WL_To) / (sigma*Cl_maxTo)))^.5;
         ML_o = VL_o/(661 * (.953^.5));
         ML_o7 = .7*ML_o;
         C = [0 0.15 .30 .45; 45500 39120 34829 31750];
         Val = interp1(C(1,:),C(2,:), ML_o7);
        
         if isnan(Val) == true
             Val = interp1(C(1,:),C(2,:), ML_o7,'linear','extrap');
         end
         WT =(Val/45550 * WTTO) + WT_Adjustment ; %Weight To Thrust Ratio 
%% WEIGHT CONVERGENCE 
        Wing_Weight =  (0.00945*(AR)^.8 * (1.35 + 1*Lambda)^.25*1*(3.75)^.5)/ ((T_C+ 0.03)^.4 * cos(Lambda)*(WL_To)^.695); %Change Top EqN
        kf = 11.5;
        l = 3.76*(Pax/Seats_Abreast) + 33.2;
        d = 1.75*Seats_Abreast +1.58*NAIS + 1;
      
        if international == 1
        l =l*1.10;
        d = d*1.10;
        end
        
        Fuselage_Weight = 0.6727*kf*l^0.6*d^0.72*3.75^0.3;
        LandingGear_Weight = 0.04;
        NacellePylon_Weight = 0.0555/WTTO;
        Tail_Weight = (0.17)*Wing_Weight;
        WingTail_Weight = Tail_Weight + Wing_Weight;
        Powerplant_Weight = 1/(3.58*WTTO);
        Fuel_Weight =1.0275*W_ratioJT9D;
        Payload_Weight = 215*Pax + W_Cargo;
        FixedEquipment_Weight = 132*Pax + 300*NE + 260*2 + 170*6; 
        Weight_Condition = 1;
        WTO_guess =200000;  
        iterationCount1 = 0;
        while Weight_Condition == true
            a = WingTail_Weight*WTO_guess^1.195;
            b = Fuselage_Weight*WTO_guess^.235;
            c = (LandingGear_Weight + NacellePylon_Weight + Powerplant_Weight + Fuel_Weight + .035 - 1)*WTO_guess;
            D_Stat = Payload_Weight + FixedEquipment_Weight;
            WTO = a+b+c+D_Stat;
            Weight_Condition = WTO >=0 ;
           
            if iterationCount1 > 20000000
               error('WEIGHT COULD NOT CONVERGE')      
            end
            
            WTO_guess = WTO_guess+1;
            iterationCount1 = iterationCount1 +1;
        end
        
        WTO = WTO_guess;
        S = WTO/ WL_To;
        b = sqrt(AR*S);
        MAC = S/b;
        T = WTO/WT;
        T_E = T/NE;
%% Drag Calculations 
        RN = (.5*Speed_of_Sound)/Kinematic_Visc;
        RN_Wing = RN*MAC;
        cf_Wing = interp1(SkinFriction(:,1),SkinFriction(:,2),RN_Wing);
        RN_Fuselage = RN*l;

        cf_Fuselage = interp1(SkinFriction(:,1),SkinFriction(:,2),RN_Fuselage);

        if isnan(cf_Fuselage) == true
            cf_Fuselage = interp1(SkinFriction(:,1),SkinFriction(:,2),RN_Fuselage, 'linear','extrap' );
        end

        S_WetWing = 2*(S - d*MAC)*1.02;
        Mo = 0.5;
        z = ((2-Mo^2)*cos(Lambda)) /sqrt(1 - Mo^2*(cos(Lambda))^2);
        K_Wing = 1 + z*(T_C) + 100*(T_C)^4;
        F_Wing = cf_Wing*S_WetWing*K_Wing;

        S_WetFuselage = 0.9*pi*d*l;
        K_Fuselage  = interp1(FuselageFormFactor(:,1),FuselageFormFactor(:,2), (l/d));
        if isnan(K_Fuselage) == true
            K_Fuselage= interp1(FuselageFormFactor(:,1),FuselageFormFactor(:,2), (l/d),'linear','extrap' );
        end

        F_Fuselage = K_Fuselage*cf_Fuselage*S_WetFuselage;

        F_Tail = 0.38*F_Wing;
        S_WetNacelles = 2.1*sqrt(T_E)*NE;
        F_Nacelles = 1.25*.0027*S_WetNacelles;
        F_Pylons = 0.2*F_Nacelles;
        F_Total = F_Wing+ F_Fuselage +F_Tail + F_Nacelles + F_Pylons;
        CD_0 = F_Total/S;
        Oswald_e = 1/(1.035+ 0.38*CD_0*pi*AR); 


%% Climb
        w_AvgClimb = ((1+0.965)/2)*(WTO);
        Climb_Velo = 1.3*(12.9/((F_Total*Oswald_e)^.25)) * sqrt(w_AvgClimb/(sigma*b));
        Tr_Climb = (sigma*F_Total*Climb_Velo^2)/296 + (94.1/(sigma*Oswald_e))*(w_AvgClimb/b)^2*(1/Climb_Velo^2);
        SFC25 = -1.54*M_cruise + 15.5;
        SFC15 = -10.5*M_cruise + 25.2;
        Ta_Climb =((T_E/45500)*(1000*(SFC25 + SFC15)));
        RateofClimb = (101*((NE*Ta_Climb)-Tr_Climb)/w_AvgClimb)*Climb_Velo;
        Time_Climb = Alt/RateofClimb;
        Range_Climb = Climb_Velo*(Time_Climb/60);
        Wf_climb = NE*Ta_Climb*.65*(Time_Climb/60);
%% Range
        W0 = WTO - Wf_climb;
        W1 = (1 - (W_ratioJT9D)) *WTO;
        CL_avg = ((W0 + W1)/(2*S))/ (1481*.2360*M_cruise^2);
        CD_i = CL_avg^2/(pi*AR*Oswald_e);
        CD = CD_0+CD_i+ .0010;
        LToDrag = CL_avg/CD;
        Tr = ((W0+W1)/2)/LToDrag;
        Tr_JT9D = (Tr*(45500/T_E))/ NE;
        SFC_35000 = -0.08*M_cruise^2 + 0.74*M_cruise;
        R_Cruise = (Climb_Velo/SFC_35000)*(LToDrag)*log((W0/W1));
        Range_Total = Range_Climb+R_Cruise;
        count = count + 1;
        RangeCondition = RAO-Range_Total>5;  %5/4
        Fuel_Adjustment = Fuel_Adjustment +.001; %5/5 Change .001
    end
%% Thrust on Top
        CL_IC  = (W0/S) / (1481*.2360*M_cruise^2);
        CDi_IC = (CL_IC)^2 / (pi * AR* Oswald_e) ;
        CD_IC  = CD_0 + CDi_IC + 0.0010;
        LD_IC  = CL_IC/CD_IC;
        TReq_IC = (W0/LD_IC)/NE;
        TReq_JT9D_IC = (Tr*(45500/T_E)); 

        if TReq_JT9D_IC > 10000
            fprintf('!!! NOT ENOUGH THRUST !!!\n')    
        end
        WT_Adjustment = WT_Adjustment + .01;
end

%% Climb Gradient

%Segment 1
CL_TO      = Cl_maxTo/(1.2)^2;
Seg1_CL_Ratio   = 1/(1.2^2);
DeltaCD_0  = interp1(CD_pTakeoff(:,1),CD_pTakeoff(:,2),Seg1_CL_Ratio);
Seg1_CD    = CD_0 +  DeltaCD_0 + .0145 +  (CL_TO/(pi * AR* Oswald_e));
Seg1_LDTO  = CL_TO/Seg1_CD;
Seg1_TReq  = WTO/Seg1_LDTO;
Seg1_TA    = (T_E/45500)*34500;
Grad_Seg1  = ((NE-1)*(Seg1_TA-Seg1_TReq )*100)/WTO;

    if NE == 2  && Grad_Seg1 < 0
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE == 3  && Grad_Seg1 < 0.3
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE >= 4 && Grad_Seg1 < 0.5
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end

% Segment 2 
Seg2_CD   = CD_0 +  DeltaCD_0 +  (CL_TO/(pi * AR* Oswald_e));
Seg2_LDTO = CL_TO/Seg2_CD;
Seg2_TReq  = WTO/Seg2_LDTO;
Grad_Seg2  = (((NE-1)*(Seg1_TA)-Seg2_TReq )*100)/WTO;

    if NE == 2  && Grad_Seg2 < 2.4
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE == 3  && Grad_Seg2 < 2.7
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE >= 4 && Grad_Seg2 < 3.0
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end
    
% Segment 3
Seg3_CL = interp1(CLmaxClean(:,1),CLmaxClean(:,2),T_C);
if isnan(Seg3_CL) == true
    Seg3_CL= interp1(CLmaxClean(:,1),CLmaxClean(:,2),T_C,'linear','extrap' );
end
       
Seg3_V  = 1.2*sqrt((296*WL_To)/.925*Seg3_CL );
Seg3_M  = Seg3_V/sqrt(1717*1.4*539.67);
Seg3_Cl = Seg3_CL/1.2^2;
Seg3_CD = CD_0 + Seg3_Cl^2/(pi*AR* Oswald_e);
Seg3_LD = Seg3_Cl/Seg3_CD;
Seg3_TReq =  WTO/Seg3_LD;
Seg3_TA    = (T_E/45500)*26500;
Grad_Seg3  = ((NE-1)*(Seg3_TA-Seg3_TReq )*100)/WTO;

    if NE == 2  && Grad_Seg3 < 1.2
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE == 3  && Grad_Seg3 < 1.5
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE >= 4 && Grad_Seg3 < 1.7
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end
 
%Approach
CL_ToAp = Cl_maxTo/(1.3)^2;
SegAp_CL_Ratio   = 1/(1.3^2);
SegAp_DeltaCD_0  = interp1(CD_pTakeoff(:,1),CD_pTakeoff(:,2),SegAp_CL_Ratio);
if  isnan(SegAp_DeltaCD_0) == true
    SegAp_DeltaCD_0= interp1(CD_pTakeoff(:,1),CD_pTakeoff(:,2),Seg3_CL_Ratio,'linear','extrap' );
end

SegAp_CD   = CD_0 +  DeltaCD_0 + (CL_ToAp/(pi * AR* Oswald_e));
SegAp_LD  = CL_ToAp/SegAp_CD;
SegAp_W    = WSLD*S;
SegAp_TReq = SegAp_W/SegAp_LD;
SegAp_V    = sqrt((296*WSLD)/(.953*CL_ToAp));
SegAp_M    = SegAp_V/659;
SegAp_TA   = (T_E/45500)*29500;
Grad_SegAp = ((((NE-1)*(SegAp_TA))-SegAp_TReq)*100)/SegAp_W;

    if NE == 2  && Grad_SegAp < 2.1
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE == 3  && Grad_SegAp < 2.4
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE >= 4 && Grad_SegAp < 2.7
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end





%Landing 
 SegLG_CL = Cl_maxLG/(1.3)^2;
 SegLG_CL_Ratio   = 1/(1.3^2);
 SegLG_DeltaCD_0  = interp1(CD_pTakeoff(:,1),CD_pTakeoff(:,2),SegLG_CL_Ratio);
if  isnan(SegLG_DeltaCD_0) == true
    SegLG_DeltaCD_0= interp1(CD_pTakeoff(:,1),CD_pTakeoff(:,2),Seg3_CL_Ratio,'linear','extrap' );
end

SegLG_CD   = CD_0 +  DeltaCD_0 + .0145 + (SegLG_CL/(pi * AR* Oswald_e));
SegLG_LD  = SegLG_CL/SegAp_CD;
SegLG_W    = WSLD*S;
SegLG_TReq = SegLG_W/SegLG_LD;
SegLG_V    = sqrt((296*WSLD)/(.953*SegLG_CL));
SegLG_M    = SegLG_V/659;
SegLG_TA   = (T_E/45500)*37200;
Grad_SegLG = ((((NE-1)*(SegLG_TA))-SegLG_TReq)*100)/SegLG_W;
    
    if NE == 2  && Grad_SegLG < 3.2
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE == 3  && Grad_SegLG < 3.2
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end


    if NE >= 4 && Grad_SegLG < 3.2
        fprintf('DOES NOT MEET MINIMUM REQUIREMENTS\n')
    end
    
%Direct
D_Statute = Range*1.15;
Ka = (7+0.015*D_Statute);
T_Cruise = ((D_Statute*1.02+20) - 1.15*(Range_Climb)) / V_cruise;
T_GroundManuever = D_Statute/(11.866+0.040669*D_Statute)/60;
V_block  = D_Statute/(T_GroundManuever + Time_Climb/60+ 0.10 + T_Cruise);
Time_block  = T_GroundManuever + (Time_Climb/60) + .10 + T_Cruise;
F_block  = Wf_climb + Tr_Climb*SFC_35000*(T_Cruise + .10);

%Flight Cost
Payload = ((205*Pax + 50*Pax) + (W_Cargo))/2000;
Dollar_BlockHr =  17.849*((V_cruise)*(WTO/10^5))^.3+40+.83;
CTM_Cr    =  Dollar_BlockHr / (V_block*Payload);
CTM_Fuel  =  (1.02*F_block*0.0438 + NE*2.15*Time_block*.135)/(D_Statute*Payload);   

% Insurance Cost 
Wa = WTO*(1-W_ratioJT9D)-(215*Pax+W_Cargo);
Ca = 2.4E6+ 87.5*Wa;
Ce = 590000 + 16*T_E;
Ct = Ca+ NE*Ce;
U = 630 + 4000/(1+1/(Time_block+.5));
CTM_Hull = (0.01*Ct)/(U*V_block*Payload);

%Direct Maintanence
KFH  = 4.9169*log10(Wa/1000)-6.425;
KFCA = .21256*log10(Wa/1000)^3.37375;
T_F  = Time_block - T_GroundManuever;
CTM_DirectEq = 8.60*(KFH*T_F+KFCA)/(V_block*Time_block*Payload);

%Airframe Material
CFHA = 1.5994*(Ca/10^6) + 3.4263;
CFCA = 1.9229*(Ca/10^6) + 2.2504;
CTM_Airframe = (CFHA*T_F+ CFHA) / (V_block*Time_block*Payload);

%Enigne Labor Cost
KFHE = NE*(T_E./10^3)./(0.82715*(T_E./10^3) + 13.639);
KFCE = .2*NE;
CTM_Engine = 8.60*(KFHE .*T_F+KFCE)/(V_block*Time_block*Payload);

%Cost of Materials
CFHE = (28.2353.*(Ce./10^6)-6.5175).*NE;
CFCE = (3.6698.*(Ce./10^6)+1.3685).*NE;
CTM_EnMaterial = (CFHE *T_F+CFCE)/(V_block*Time_block*Payload);

CTM_Total = 2*(CTM_DirectEq+CTM_Airframe+CTM_Engine+CTM_EnMaterial);
CTM_Depreciation = 1/(V_block*Payload)*(Ct+.06*(Ct-NE*Ce)+.3*NE*Ce)/(14*U);

Total_DOC = CTM_Cr+ CTM_Fuel+ CTM_Hull+ CTM_Total + CTM_Depreciation;
DOC_Payload = Total_DOC*Payload/Pax;


% if rad2deg(Lambda) == 10
%     WTO_saved10(ARcount) = WTO;
%     Total_DOC_saved10(ARcount) = Total_DOC;
% end
% 
% if rad2deg(Lambda) > 10 &&  rad2deg(Lambda)< 20
%     WTO_saved15(ARcount) = WTO;
%     Total_DOC_saved15(ARcount) = Total_DOC;
% end
% 
% if rad2deg(Lambda) == 20
%     WTO_saved20(ARcount) = WTO;
%     Total_DOC_saved20(ARcount) = Total_DOC;
% end
% 
% if rad2deg(Lambda) == 25
%     WTO_saved25(ARcount) = WTO;
%     Total_DOC_saved25(ARcount) = Total_DOC;
% end
% 
% if rad2deg(Lambda) > 25 &&  rad2deg(Lambda)< 35
%     WTO_saved30(ARcount) = WTO;
%     Total_DOC_saved30(ARcount) = Total_DOC;
% end
% 
% if rad2deg(Lambda) == 35
%     WTO_saved35(ARcount) = WTO;
%     Total_DOC_saved35(ARcount) = Total_DOC;
% end
% 
% 
% end
% %% Data Storage
% 
% end
% figure(1)
% plot(Ar,WTO_saved10,'og')
% hold on
% plot(Ar,WTO_saved20,'ob')
% plot(Ar,WTO_saved35,'oc')
% plot(Ar,WTO_saved15,'om')
% plot(Ar,WTO_saved30,'oy')
% plot(Ar,WTO_saved25,'ok')
% title('Take off Weight vs. Aspect Ratio')
% xlabel('Aspect Ratio')
% ylabel('Takeoff Weight')
% legend('25 \lambda', '30 \lambda', '10 \lambda', '20 \lambda', '15 \lambda', '35 \lambda')
% hold off
% 
% figure(2)
% plot(Ar,Total_DOC_saved10,'og')
% hold on
% plot(Ar,Total_DOC_saved20,'ob')
% plot(Ar,Total_DOC_saved35,'oc')
% plot(Ar,Total_DOC_saved15,'om')
% plot(Ar,Total_DOC_saved30,'oy')
% plot(Ar,Total_DOC_saved25,'ok')
% title('Direct Operating Cost vs. Aspect Ratio')
% xlabel('Aspect Ratio')
% ylabel('Direct Operating Cost')
% legend('25 \lambda', '30 \lambda', '10 \lambda', '20 \lambda', '30 \lambda', '35 \lambda')
% hold off
if rad2deg(Lambda) == (0)
WTO_saved0(ARcount) = WTO ;
Total_DOC_saved0(ARcount) = Total_DOC;
end

if rad2deg(Lambda) == (10)
WTO_saved10(ARcount) = WTO;
Total_DOC_saved10(ARcount) = Total_DOC;
end

if (Lambda) == deg2rad(15)
WTO_saved15(ARcount) = WTO;
Total_DOC_saved15(ARcount) = Total_DOC;
end

if rad2deg(Lambda) == (20)
WTO_saved20(ARcount) = WTO;
Total_DOC_saved20(ARcount) = Total_DOC;
end

if rad2deg(Lambda) == (25)
WTO_saved25(ARcount) = WTO;
Total_DOC_saved25(ARcount) = Total_DOC;
end

if (Lambda) == deg2rad(30)
WTO_saved30(ARcount) = WTO;
Total_DOC_saved30(ARcount) = Total_DOC;
end

if rad2deg(Lambda) == (35)
WTO_saved35(ARcount) = WTO;
Total_DOC_saved35(ARcount) = Total_DOC;
end

if rad2deg(Lambda) ==(40)
WTO_saved40(ARcount) = WTO;
Total_DOC_saved40(ARcount) = Total_DOC;
end
end
end
%% Data Storage

figure(1)

plot(Ar,WTO_saved10,'-og')
hold on 
plot(Ar,WTO_saved15,'-ob')
plot(Ar,WTO_saved20,'-oc')
plot(Ar,WTO_saved25,'-om')
plot(Ar,WTO_saved30,'-oy')
plot(Ar,WTO_saved35,'-ok')
% plot(Ar,WTO_saved40,'-+r')
title('Take off Weight vs. Aspect Ratio')
xlabel('Aspect Ratio')
ylabel('Takeoff Weight')
legend( '10 \lambda', '15 \lambda', '20 \lambda', '25 \lambda', '30 \lambda', '35 \lambda')   
hold off 

figure(2)
plot(Ar,Total_DOC_saved10,'-og')
hold on 
plot(Ar,Total_DOC_saved15,'-ob')
plot(Ar,Total_DOC_saved20,'-oc')
plot(Ar,Total_DOC_saved25,'-om')
plot(Ar,Total_DOC_saved30,'-oy')
plot(Ar,Total_DOC_saved35,'-ok')
% plot(Ar,Total_DOC_saved40,'-+r')
title('Direct Operating Cost vs. Aspect Ratio')
xlabel('Aspect Ratio')
ylabel('Direct Operating Cost')
legend( '10 \lambda', '15 \lambda', '20 \lambda', '25 \lambda', '30 \lambda', '35 \lambda)') 
hold off 
