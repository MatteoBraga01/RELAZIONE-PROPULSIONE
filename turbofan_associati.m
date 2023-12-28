%% TURBOFAN A FLUSSI ASSOCIATI
% fuzione che determina parametri di interesse, in funzione del BPR e della velocità di volo,
% per poi essere plottati in un ulteriore script
%ATTENZIONE: analisi in regime non supersonico

% INPUT

    % M=mach di volo
    % T_a=la temperatura ambientale  nelle condizioni di volo[K]
    % p_a=la pressione ambientale  nelle condizioni di volo[Pa]
    % ma= portata d'aria totalen processata [kg/s]
    
    % FAN:
    % BPR: bypass rateo
    
    % COMPRESSORE
    % B_c: rapporto di compressione del compressore
    
    % CAMERA DI COMBUSTIONE:
    % T_max: temperatura massima camera di combustione [K]
    
    % POST BRUCIATORE
    % T_max_PB: temperatura massima del postbruciatore [K]

    % RENDIMENTI:
    % n_presa:         rendimento pneumatico presa
    % n_fan:           rendimento adiabatico fan
    % n_m_fan:         rendimento meccanico fan
    % n_compressore:   rendimento adiabatico compressore
    % n_m_compressore: rendimento meccanico compressore
    % n_combustione:   rendimento pneumatico camera di combustione
    % n_turbina:       rendimento adiabatico turbina
    % n_m_turbina:     rendimento meccanico turbina (uguale per entrambe)
    % n_ugello         rendimento adiabatico ugello (uguale per entrambi)
    % n_combustione_PB:   rendimento pneumatico post bruciatore
%OUTPUT:

    %f: rapporto combustibile aria totale
    %f_1: rapporto combustibile aria non afterburner
    %f_2: rapporto combustibile aria afterburner
    %T: trust [N]
    %TSFC: trust specific fuel consumption [kg/N*s]
    %n_p: rendimento propulsivo
    %n_th: rendimento termodinamico
    %n_o: rendimento globale
    %B_f: rapporto di compressione del fan

%indicazioni: 
% pedice _01 fa riferimento alle condizioni totali
% pedice _1  fa riferimento alle condizioni dinamiche

function [f,f_1,f_2, T, TSFC, n_p, n_th, n_o,B_f] = turbofan_associati(M, T_a, p_a, ma, BPR, B_c, T_max, ...
                                   n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                   n_combustione,n_turbina,n_m_turbina, n_ugello, ...
                                   T_max_PB, n_combustione_PB)

%definizioni costanti:
R=287;   %costante gas aria  [J/kg*K] 
cp=1004; %calore specifico aria (a press con)[J/kg*K] 
gamma=1.4; %per gas biatomico []

%STEP 0: calcolo del rapporto di compressione in funzione del BPR


%STEP1: calcolo velocità di volo
v0=M*sqrt(gamma*R*T_a);  %[m/s] 

%STEP 2: 1-2: PRESA DINAMICA
    %calcolo temperature e pressioni totali in ingresso
    T_01=T_a*(1+(gamma-1)/2*M^2); 
    p_01=p_a*(1+(gamma-1)/2*M^2)^(gamma/(gamma-1));

    %calcolo temperature e pressioni totali in uscita
    T_02=T_01;
    p_02=n_presa*p_01;   

    B_f=0;
    e=inf;
    
%ciclo iterativo che ha come fine quello di ricavare il rapporto di compressione del fan
%dato il valore del BPR
while e>0.5
    
    B_f= B_f+0.000001;
    
    %STEP 3: 2-21/13: FAN    
    
        %calcolo temperatura in uscita:
        T_021_ideale=T_02*B_f^((gamma-1)/gamma);   %temperatura ideale
        T_021=T_02+(T_021_ideale-T_02)/n_fan;        %temperatura reale 
        T_013=T_021;
    
        %calcolo pressione in uscita :
        p_021=B_f*p_02;   
        p_013=p_021;
    
    %STEP 4: 21-3: COMPRESSORE  
        %calcolo temperatura in uscita:
        T_03_ideale=T_021*B_c^((gamma-1)/gamma);   %temperatura ideale
        T_03=T_021+(T_03_ideale-T_021)/n_compressore;        %temperatura reale 
        
        %calcolo pressione in uscita :
        p_03=B_c*p_021;   
        
    %STEP 5: 3-4: CAMERA DI COMBUSTIONE
        Hi=43000000;         %potere calorifico inferiore del combustibile JP5 [J/kg]
        T_04=T_max;          %temperatura in uscita dal combustore fissata
        
        %calcolo rapporto combustibile/ flusso prinicpale
        f_1=(cp*(T_04-T_03)) / (Hi-cp*T_04);  %rapporto tra mf e ma;
        %calcolo pressione in uscita :
        p_04=n_combustione*p_03;
        
    %STEP 5: 4-41: HIGH PRESSURE TURBINE
        %calcolo temperatura in uscita
        T_041=T_04-(T_03-T_021)/((1+f_1)*n_m_compressore*n_m_turbina); 
        % calcolo pressione in uscita
        T_041_ideale=T_04-(T_04-T_041)/n_turbina;  
        p_041=p_04*(T_041_ideale/T_04)^(gamma/(gamma-1));  
    
    %STEP 6: 41-5: LOW PRESSURE TURBINE
        %calcolo temperatura in uscita
        T_05=T_041-(T_021-T_02)*(1+BPR)/((1+f_1)*n_m_fan*n_m_turbina);  
        % calcolo pressione in uscita
        T_05_ideale=T_041-(T_041-T_05)/n_turbina;  
        p_05=p_041*(T_05_ideale/T_041)^(gamma/(gamma-1)); 
    
        e=abs(p_013-p_05);
        
end



%STEP 7: 5-6: MIXER
    %calcolo temperatura statica in uscita
    T_06= ((1+f_1)*T_05+BPR*T_013)/(1+f_1+BPR);

    p_06=p_05;

%STEP 8: 6-7 valutare presenza POSTBRUCIATORE
    if nargin>16
        Hi=43000000;         %potere calorifico inferiore del combustibile JP5 [J/kg]
        T_07=T_max_PB;          %temperatura in uscita dal combustore fissata
        
        %calcolo rapporto combustibile postbruciatore/ flusso 
        f_2=(1+BPR+f_1)*cp*(T_07-T_06)/(Hi-cp*T_07);  %rapporto tra mf e ma;
        %calcolo pressione in uscita :
        p_07=n_combustione_PB*p_06;
        %rapporto comb flusso totale
        f=f_1+f_2;
    else
        T_07=T_06;
        p_07=p_06;
        f=f_1;
        f_2=0;
    end

%STEP 9: 7-9: UGELLO 
    %calcolo temperatura statica in uscita
    T_9_ideale=T_07*(p_a/p_07)^((gamma-1)/gamma) ;  
    T_9=T_07-n_ugello*(T_07-T_9_ideale);

    %velocità efflusso ugello:
    ve=sqrt(2*cp*n_ugello*(T_07-T_9_ideale));  

%SPINTA E SPINTA SPECIFICA 
T=ma*((1+f+BPR)*ve-(1+BPR)*v0)  %[N] 

%TSFC
TSFC=f*ma/T %[kg/N*s] 

%RENDIMENTI
    %rendimento propulsivo
    n_p=(2*v0*((1+f+BPR)*ve-(1+BPR)*v0)) / ((1+f+BPR)*ve^2-(1+BPR)*v0^2); 
    %rendimento termodinamico
    n_th=((1+f+BPR)*ve^2-(1+BPR)*v0^2) / (2*f*Hi); 
    %rendimento globale
    n_o=n_th*n_p;   

end


