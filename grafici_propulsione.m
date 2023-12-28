%% GRAFICI PER RELAZIONE PROPULSIONE
%% rendimento propulsivo vs portata d'aria

%fissati la spinta e la velocità di volo V0 si pone attenzione
%sull'andamento del rendimento propulsivo vs portata d'aria

clc; clear

T=500*10^3; %N
V0=700*10^3/3600; %m/s
m_a=linspace(1,10000,10000);%portata d'aria in ingresso in kg/s
n_P=zeros(length(m_a),1);
for i=1:length(m_a)
    n_P(i)=1/(1+(T/(2*m_a(i)*V0)));
end

plot(m_a, n_P, LineWidth=3)
grid on;
xlabel('Portata aria [kg/s]', 'FontSize',14);
ylabel('Rendimento Propulsivo Unitario [/]', FontSize=14)

%% PLOT GRAFICI FLUSSI SEPARATI
clc; clear;
M=[0.1,0.2,0.4,0.6,0.8,1];
%studiati a temperatura e bpressione a 10 km di quota
T_a= 223; 
p_a= 0.26*10^5;
ma=60; %valore approssimativo, dipende da diversi altri fattori
B_f=1.5;
B_c=20;
BPR= [0.1:0.1:10];
T_max=1400;

%rendimenti, presi valori verosimili
n_presa=0.99;
n_fan=0.95; 
n_m_fan=0.99; 
n_compressore=0.9;
n_m_compressore=0.99;
n_combustione=0.95;
n_turbina= 0.94;
n_m_turbina=0.99;
n_ugello=0.99;

% inizializzazione valori
TSFC_mat=[];
n_p_mat=[];
n_th_mat=[];
n_o_mat=[];

for i=1:length(M)

    TSFC_vect=[];
    n_p_vect=[];
    n_th_vect=[];
    n_o_vect=[];

    for j=1:length(BPR)
        [f, T, TSFC, n_p, n_th, n_o] = turbofan_separati(M(i), T_a, p_a, ma, B_f, BPR(j), B_c, T_max, ...
                                       n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                       n_combustione,n_turbina,n_m_turbina, n_ugello);
        TSFC_vect=[TSFC_vect;TSFC];
        n_p_vect=[n_p_vect;n_p];
        n_th_vect=[n_th_vect;n_th];
        n_o_vect=[n_o_vect;n_o];
    end
    TSFC_mat=[TSFC_mat, TSFC_vect];
    n_p_mat=[n_p_mat, n_p_vect];
    n_th_mat=[n_th_mat, n_th_vect];
    n_o_mat=[n_o_mat, n_o_vect];
    
end



%% PLOT
map= colormap(parula(length(M)+1));

% TSFC
figure(1)
for i=1:length(M)
    plot(BPR(:),TSFC_mat(:,i),LineWidth=2, Color=map(i,:))
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("TSFC[kg/N*s]",FontSize=15)
legend("M=0.1","M=0.2","M=0.4","M=0.6","M=0.8","M=1",Orientation="horizontal")

%%
%RENDIMENTI
tiledlayout(1,3)

nexttile()
for i=1:length(M)
    plot(BPR(:),n_th_mat(:,i),LineWidth=2, Color=map(i,:))
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("\eta_t_h",FontSize=15)
legend("M=0.1","M=0.2","M=0.4","M=0.6","M=0.8","M=1",Orientation="vertical")

nexttile()
for i=1:length(M)
    plot(BPR(:),n_p_mat(:,i),LineWidth=2, Color=map(i,:))
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("\eta_p",FontSize=15)


nexttile()
for i=1:length(M)
    plot(BPR(:),n_o_mat(:,i),LineWidth=2, Color=map(i,:))
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("\eta_o",FontSize=15)

%% PLOT GRAFICI FLUSSI ASSOCIATI (usare grafici sopra)
clc; clear;
M=[0.1,0.2,0.4,0.6,0.8,0.9];
%studiati a temperatura e bpressione a 10 km di quota
T_a= 223; 
p_a= 0.26*10^5;
ma=60; %valore approssimativo, dipende da diversi altri fattori
B_c=20;
BPR= [0.1:0.1:10];
T_max=1400;

%rendimenti, presi valori verosimili
n_presa=0.99;
n_fan=0.95; 
n_m_fan=0.99; 
n_compressore=0.9;
n_m_compressore=0.99;
n_combustione=0.95;
n_turbina= 0.94;
n_m_turbina=0.99;
n_ugello=0.99;

% inizializzazione valori
TSFC_mat=[];
n_p_mat=[];
n_th_mat=[];
n_o_mat=[];
B_f_mat=[];

for i=1:length(M)

    TSFC_vect=[];
    n_p_vect=[];
    n_th_vect=[];
    n_o_vect=[];
    B_f_vect=[];

    for j=1:length(BPR)
        [f, f_1, f_2, T, TSFC, n_p, n_th, n_o, B_f] = turbofan_associati(M(i), T_a, p_a, ma, BPR(j), B_c, T_max, ...
                                   n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                   n_combustione,n_turbina,n_m_turbina, n_ugello);
        TSFC_vect=[TSFC_vect;TSFC];
        n_p_vect=[n_p_vect;n_p];
        n_th_vect=[n_th_vect;n_th];
        n_o_vect=[n_o_vect;n_o];
        B_f_vect=[B_f_vect,B_f];
    end
    TSFC_mat=[TSFC_mat, TSFC_vect];
    n_p_mat=[n_p_mat, n_p_vect];
    n_th_mat=[n_th_mat, n_th_vect];
    n_o_mat=[n_o_mat, n_o_vect];
    B_f_mat=[B_f_mat, B_f_vect];

end
%% rapporto BPR/ B_f (a mach 0.4)
clc; clear;
M=0.4;
%studiati a temperatura e bpressione a 10 km di quota
T_a= 223; 
p_a= 0.26*10^5;
ma=60; %valore approssimativo, dipende da diversi altri fattori
B_c=20;
BPR= [0.2:0.1:15];
T_max=1450;

%rendimenti, presi valori verosimili
n_presa=0.90;
n_fan=0.88; 
n_m_fan=0.99; 
n_compressore=0.9;
n_m_compressore=0.99;
n_combustione=0.94;
n_turbina= 0.92;
n_m_turbina=0.99;
n_ugello=0.98;

% inizializzazione valori
B_f_vect=[];

for j=1:length(BPR)
    [f, f_1, f_2, T, TSFC, n_p, n_th, n_o, B_f] = turbofan_associati(M, T_a, p_a, ma, BPR(j), B_c, T_max, ...
                               n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                               n_combustione,n_turbina,n_m_turbina, n_ugello);
    B_f_vect=[B_f_vect,B_f];
end
 %%
figure()

plot(BPR(:),B_f_vect(:),LineWidth=2)
hold on;

grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("B_f",FontSize=15)
legend("M=0.4");
%% confronto turbofan associati, separati a M=0.5
clc; clear;

M=0.5;
%studiati a temperatura e bpressione a 10 km di quota
T_a= 223; 
p_a= 0.26*10^5;
ma=60; %valore approssimativo, dipende da diversi altri fattori
B_c=20;
BPR= [0.1:0.1:10];
T_max=1400;

%rendimenti, presi valori verosimili
n_presa=0.99;
n_fan=0.95; 
n_m_fan=0.99; 
n_compressore=0.9;
n_m_compressore=0.99;
n_combustione=0.95;
n_turbina= 0.94;
n_m_turbina=0.99;
n_ugello=0.99;

TSFC_vect_AS=[];
n_p_vect_AS=[];
n_th_vect_AS=[];
T_AS=[];

TSFC_vect_SEP=[];
n_p_vect_SEP=[];
n_th_vect_SEP=[];
T_SEP=[];

for j=1:length(BPR)
        [f, f_1, f_2, T, TSFC, n_p, n_th, n_o, B_f] = turbofan_associati(M, T_a, p_a, ma, BPR(j), B_c, T_max, ...
                                   n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                   n_combustione,n_turbina,n_m_turbina, n_ugello);
        TSFC_vect_AS=[TSFC_vect_AS;TSFC];
        n_p_vect_AS=[n_p_vect_AS;n_p];
        n_th_vect_AS=[n_th_vect_AS;n_th];
        T_AS=[T_AS;T];
        
        [f, T, TSFC, n_p, n_th, n_o] = turbofan_separati(M, T_a, p_a, ma, B_f, BPR(j), B_c, T_max, ...
                                       n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                       n_combustione,n_turbina,n_m_turbina, n_ugello);
        TSFC_vect_SEP=[TSFC_vect_SEP;TSFC];
        n_p_vect_SEP=[n_p_vect_SEP;n_p];
        n_th_vect_SEP=[n_th_vect_SEP;n_th];
        T_SEP=[T_SEP;T];
end

%%
% TSFC
tiledlayout(2,2)

nexttile(2)
for i=1:length(M)
    plot(BPR(:),TSFC_vect_AS(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),TSFC_vect_SEP(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("TSFC[kg/N*s]",FontSize=15)
legend("associati","separati",Orientation="vertical")

% n_th
nexttile(4)
for i=1:length(M)
    plot(BPR(:),n_th_vect_AS(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),n_th_vect_SEP(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("\eta_t_h",FontSize=15)
legend("associati","separati",Orientation="vertical")

% n_P
nexttile(3)
for i=1:length(M)
    plot(BPR(:),n_p_vect_AS(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),n_p_vect_SEP(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("\eta_p",FontSize=15)
legend("associati","separati",Orientation="vertical")

%T
nexttile(1)
for i=1:length(M)
    plot(BPR(:),T_AS(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),T_SEP(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("T [N]",FontSize=15)
legend("associati","separati",Orientation="vertical")

%% confronto con e senza post bruciatore
clc; clear;

M=0.8;
%studiati a temperatura e bpressione a 10 km di quota
T_a= 223; 
p_a= 0.26*10^5;
ma=60; %valore approssimativo, dipende da diversi altri fattori
B_c=20;
BPR= [0.1:0.1:2];
T_max=1400;
T_max_PB=1600;

%rendimenti, presi valori verosimili
n_presa=0.99;
n_fan=0.95; 
n_m_fan=0.99; 
n_compressore=0.9;
n_m_compressore=0.99;
n_combustione=0.95;
n_turbina= 0.94;
n_m_turbina=0.99;
n_ugello=0.99;
n_PB=0.98;

TSFC_vect_AS=[];
T_AS=[];
f_AS=[];

f_1_v=[];
f_2_v=[];

for j=1:length(BPR)
        [f, f_1, f_2, T, TSFC, n_p, n_th, n_o, B_f] = turbofan_associati(M, T_a, p_a, ma, BPR(j), B_c, T_max, ...
                                   n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                   n_combustione,n_turbina,n_m_turbina, n_ugello);
        TSFC_vect_AS=[TSFC_vect_AS;TSFC];
        T_AS=[T_AS;T];
        f_AS=[f_AS,f];
end 

TSFC_vect_PB=[];
T_PB=[];
f_PB=[];

for j=1:length(BPR)
        [f, f_1, f_2, T, TSFC, n_p, n_th, n_o] = turbofan_associati(M, T_a, p_a, ma, BPR(j), B_c, T_max, ...
                                       n_presa,n_fan, n_m_fan, n_compressore,n_m_compressore, ...
                                       n_combustione,n_turbina,n_m_turbina, n_ugello,T_max_PB,n_PB);
        TSFC_vect_PB=[TSFC_vect_PB;TSFC];
        T_PB=[T_PB;T];
        f_PB=[f_PB,f];
        f_1_v=[f_1_v;f_1];
        f_2_v=[f_2_v;f_2];
end

%%
% TSFC
figure(1)
for i=1:length(M)
    plot(BPR(:),TSFC_vect_AS(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),TSFC_vect_PB(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("TSFC[kg/Ns]",FontSize=15)
legend("post bruciatore OFF","post bruciatore ON",Orientation="vertical")

% T
figure(2)
for i=1:length(M)
    plot(BPR(:),T_AS(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),T_PB(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("Trust[N]",FontSize=15)
legend("post bruciatore OFF","post bruciatore ON",Orientation="vertical")

% f
figure(3)
for i=1:length(M)
    plot(BPR(:),f_1_v(:,i),LineWidth=2)
    hold on;
    plot(BPR(:),f_2_v(:,i),LineWidth=2)
    hold on;
end
grid on,grid minor
xlabel("BPR","FontSize",15);
ylabel("rapporto aria combustibile",FontSize=15)
legend("f_1","f_2",Orientation="vertical")


