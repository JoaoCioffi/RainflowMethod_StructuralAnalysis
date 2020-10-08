clear all
close all
clc
%%
%
Data = load('Pontos Rainflow.txt'); %Cycle Data 
Size_Data = size(Data);
%
%--------------------------------------------------------------------------
%
%
Horas_Voo = 6;
%                             // Plot //  
time_pass = (Horas_Voo/Size_Data(1,1)); %passo de tempo considerando um ciclo GAG de ~6h
t = (0 + time_pass):time_pass:Horas_Voo; %tempo do ciclo GAG 
s = Data; 
%
figure(1)
hold on; grid on; grid minor;
plot(t,s,'LineWidth',0.5,'Color','b')
xlabel('Time [h]','FontSize',12,'FontWeight','bold','Color','r');
ylabel('\sigma [ksi]','FontSize',12,'FontWeight','bold','Color','r');
title('_._:_:_. Ground-Air-Ground Fatigue Spectrum _._:_:_.'); 
%
%--------------------------------------------------------------------------
%
%Picos:
Peaks = findpeaks(Data); %Peaks values in Data array
%
%Vales:
vls = islocalmin(Data); %logical value for local min (valley)
find_valley = find(vls==1); %find logical occurence for local min (Valley positions)
Valleys = Data(find_valley,1); %Valley values in Data array
%
Peaks_Valleys = [Peaks Valleys]; %Matriz Picos e Vales
%
%Tabulação das Posições: Picos x Vales : 
fprintf('\n')
fprintf('\t -------------------------------------------------------------')
fprintf('\n')
fprintf('\t\t\t .::. Tabulação Considerando Picos e Vales .::.')
fprintf('\n')
fprintf('\n')
Tabulation_Peaks_Valleys = array2table(Peaks_Valleys,...
    'VariableNames',{'Picos','Vales'});
disp(Tabulation_Peaks_Valleys);
%
%% RAINFLOW - IMPLEMENTAÇÃO 
%Função Rainflow: Baseado no modelo ASTM E1049-85(2017), porém considera o rainflow "completo"
[c,hist,edges,rmm,idx] = rainflow(Data); %Compute cycle counts for Data. Display the matrix of cycle counts.
c(:,1)=c(:,1).*2; %Aproximação para o método de ocorrências, i.e, considera o rainflow "simplificado"
%
%Tabulação Rainflow : 
fprintf('\n')
fprintf('\t -------------------------------------------------------------')
fprintf('\n')
fprintf('\t\t\t .::. Tabulação do Método Rainflow .::.')
fprintf('\n')
fprintf('\n')
fprintf('\n')
Rainflow_Sample = array2table(c,'VariableNames',{'Cicle_Count','Range','Mean','Start','End'});
Size_Rainflow_Sample = size(Rainflow_Sample);
disp(Rainflow_Sample);
%
figure(2)
histogram('BinEdges',edges','BinCounts',sum(hist,2))
xlabel('Stress Range')  
ylabel('Cycle Counts')
%
%% CÁLCULO DE VIDA ÚTIL:
%
%--------------------------------------------------------------------------
%
%                   // Outputs Método Rainflow //
ll = 1:1:Size_Rainflow_Sample(1,1);
%
Cycle_Count_Table = Rainflow_Sample(ll,1); %Forma Tabular
Cycle_Count_Array = Cycle_Count_Table{:,:}; %Converte a Forma Tabular em Array
%
Range_Table = Rainflow_Sample(ll,2); %Forma Tabular
Range_Array = Range_Table{:,:}; %Converte a Forma Tabular em Array
%
Mean_Table = Rainflow_Sample(ll,3); %Forma Tabular
Mean_Array = Mean_Table{:,:}; %Converte a Forma Tabular em Array
%
Start_Table = Rainflow_Sample(ll,4); %Forma Tabular
Start_Array = Start_Table{:,:}; %Converte a Forma Tabular em Array
%
End_Table = Rainflow_Sample(ll,5); %Forma Tabular
End_Array = End_Table{:,:}; %Converte a Forma Tabular em Array
%
%--------------------------------------------------------------------------
%
%                            // Tensões //
Sigma_Max_Array = (Mean_Array + Range_Array); %Forma Array
Sigma_Max_Table = array2table(Sigma_Max_Array); %Converte a Forma Array em Tabular
%
Sigma_Min_Array = (Mean_Array - Range_Array); %Forma Array
Sigma_Min_Table = array2table(Sigma_Min_Array); %Converte a Forma Array em Tabular
%
Sigma_Med_Array = 0.5*(Sigma_Max_Array + Sigma_Min_Array); %Forma Array
Sigma_Med_Table = array2table(Sigma_Med_Array); %Converte a Forma Array em Tabular
%
%
R_Array = Sigma_Min_Array./Sigma_Max_Array; %Forma Array
R_Table = array2table(R_Array); %Converte a Forma Array em Tabular
%
%--------------------------------------------------------------------------
%
%                              // SWT // 
Sigma_eq_SWT_Array = Sigma_Max_Array.*sqrt((1-R_Array)./2); %Forma Array
Size_Sigma_eq_SWT_Array = size(Sigma_eq_SWT_Array);
Sigma_eq_SWT_Table = array2table(Sigma_eq_SWT_Array); %Converte a Forma Array em Tabular
%
%--------------------------------------------------------------------------
%
%             // Palmgren Miner: Número de Ciclos p/ falhar // 
%
N = Cycle_Count_Array;
Array_Pgm = [Sigma_Min_Array Sigma_Max_Array R_Array Sigma_eq_SWT_Array N]; %Array dados Palmgren Miner
fprintf('\n')
fprintf('\t -------------------------------------------------------------')
fprintf('\n')
fprintf('\t\t\t .::. Tabulação Palmgren Miner .::.')
fprintf('\n')
fprintf('\n')
fprintf('\n')
Tabulation_Pgm_Painel_Inteiro = array2table(Array_Pgm,'VariableNames',{'S_min','S_max','R','S_eq_SWT','N'});
disp(Tabulation_Pgm_Painel_Inteiro);
fprintf('\n')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
%              log(N)=11.10-3.97log(S-15.8): Al 2024-T3
%
Sigma_u = 72; %ultimate tensile strength [ksi] for Al 2024-T3
FS = 1.5; %Fator de Segurança
Limit_Clipping = Sigma_u*FS;
Limit_Truncation = 15.8;
%
i=1:1:Size_Sigma_eq_SWT_Array(1,1);
%
%Clipping: Desprezar tensões maiores que Sigma_u*FS
%Truncamento: Desprezar tensões menores de 15.8 ksi (valor retirado da equação da curva S-N do material)
%
j = find(Sigma_eq_SWT_Array(i,1) > Limit_Truncation & Sigma_eq_SWT_Array(i,1) < Limit_Clipping); %posições de Sigma_eq_SWT_Array que atendem às restrições de clipagem e truncamento
%
Sigma_eq_SWT_Array_Remastered = Sigma_eq_SWT_Array(j,Size_Sigma_eq_SWT_Array(1,2));
N_Remastered = N(j,1);
%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
%Painel Inteiriço:
Nf_Painel_Inteiro = 10 .^ (11.1 - (3.97.*log10(Sigma_eq_SWT_Array_Remastered - 15.8)) ); %Curva S Vs. log(N): Al 2024-T3
N_Nf_Painel_Inteiro = N_Remastered./Nf_Painel_Inteiro; %Razão N/Nf Palmgren Miner
%
D_Painel_Inteiro = sum(N_Nf_Painel_Inteiro);
Number_of_Cycles_to_Fail_Painel_Inteiro = real(1/D_Painel_Inteiro); %Número de Ciclos GAG (ou Voos) até a falha
Hours_to_Fail_Painel_Inteiro = real(Number_of_Cycles_to_Fail_Painel_Inteiro * Horas_Voo); %Quantidade de horas GAG (ou horas de voo) até a falha
%
fprintf('\t -------------------------------------------------------------')
fprintf('\n')
fprintf('\t\t\t .::. Vida Útil dos Painéis .::.')
fprintf('\n')
fprintf('\t\t\t\t\t\t    .')
fprintf('\n')
fprintf('\t\t\t\t\t\t    .')
fprintf('\n')
fprintf('\t\t\t\t\t\t    .')
fprintf('\n')
fprintf('\t //////////')
fprintf('\n')
fprintf('>> Painel Inteiro:')
fprintf('\n')
fprintf('\n')
disp(Number_of_Cycles_to_Fail_Painel_Inteiro);
fprintf('(Ciclos GAG | Vôos)')
fprintf('\n')
fprintf('\n')
fprintf('. . . . . . . . . . . . .')
fprintf('\n')
fprintf('\t\t ou')
fprintf('\n')
fprintf('. . . . . . . . . . . . .')
fprintf('\n')
fprintf('\n')
disp(Hours_to_Fail_Painel_Inteiro);
fprintf('(Horas)')
fprintf('\n')
%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
%Painel Alívio Circular:
Kt_Circular = 3; %Fator de Concentração de Tensão (Kt ~3)
Nf_Painel_Circular = 10 .^ (11.1 - (3.97.*log10((Kt_Circular.*Sigma_eq_SWT_Array_Remastered) - 15.8)) ); %Curva S Vs. log(N): Al 2024-T3
N_Nf_Painel_Circular = N_Remastered./Nf_Painel_Circular; %Razão N/Nf Palmgren Miner
%
D_Painel_Circular = sum(N_Nf_Painel_Circular);
Number_of_Cycles_to_Fail_Painel_Circular = real(1/D_Painel_Circular); %Número de Ciclos GAG (ou Voos) até a falha
Hours_to_Fail_Painel_Circular = real(Number_of_Cycles_to_Fail_Painel_Circular * Horas_Voo); %Quantidade de horas GAG (ou horas de voo) até a falha
%
fprintf('\n')
fprintf('\t //////////')
fprintf('\n')
fprintf('>> Painel Alívio Circular:')
fprintf('\n')
fprintf('\n')
disp(Number_of_Cycles_to_Fail_Painel_Circular);
fprintf('(Ciclos GAG | Vôos)')
fprintf('\n')
fprintf('\n')
fprintf('. . . . . . . . . . . . .')
fprintf('\n')
fprintf('\t\t ou')
fprintf('\n')
fprintf('. . . . . . . . . . . . .')
fprintf('\n')
fprintf('\n')
disp(Hours_to_Fail_Painel_Circular);
fprintf('(Horas)')
fprintf('\n')
%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
%Painel Alívio Quadrado:
Kt_Quadrado = 7; %Fator de Concentração de Tensão (Kt 6~8)
Nf_Painel_Quadrado = 10 .^ (11.1 - (3.97.*log10((Kt_Quadrado.*Sigma_eq_SWT_Array_Remastered) - 15.8)) ); %Curva S Vs. log(N): Al 2024-T3
N_Nf_Painel_Quadrado = N_Remastered./Nf_Painel_Quadrado; %Razão N/Nf Palmgren Miner
%
D_Painel_Quadrado = sum(N_Nf_Painel_Quadrado);
Number_of_Cycles_to_Fail_Painel_Quadrado = real(1/D_Painel_Quadrado); %Número de Ciclos GAG (ou Voos) até a falha
Hours_to_Fail_Painel_Quadrado = real(Number_of_Cycles_to_Fail_Painel_Quadrado * Horas_Voo); %Quantidade de horas GAG (ou horas de voo) até a falha
%
fprintf('\n')
fprintf('\t //////////')
fprintf('\n')
fprintf('>> Painel Alívio Quadrado:')
fprintf('\n')
fprintf('\n')
disp(Number_of_Cycles_to_Fail_Painel_Quadrado);
fprintf('(Ciclos GAG | Vôos)')
fprintf('\n')
fprintf('\n')
fprintf('. . . . . . . . . . . . .')
fprintf('\n')
fprintf('\t\t ou')
fprintf('\n')
fprintf('. . . . . . . . . . . . .')
fprintf('\n')
fprintf('\n')
disp(Hours_to_Fail_Painel_Quadrado);
fprintf('(Horas)')
fprintf('\n')
%
% userview = memory;
% [userview,systemview] = memory
