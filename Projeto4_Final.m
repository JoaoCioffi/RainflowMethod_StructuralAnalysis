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
%%
ksi_to_MPa = 6.8947; %Sigma values' conversion factor [ksi -> MPa]
%
%Passando de ocorrências para excedências:
occurrences = [N Sigma_eq_SWT_Array]; %Tabela de ocorrências (N x Sigma_eq) => output do método Rainflow
sum_N_Excedencias = occurrences; 
%
 for cont = 1:length(occurrences)
     termo = length(occurrences) - cont+ 1;
     sum_N_Excedencias(termo,1) = sum( occurrences(termo:length(occurrences)));
 end
%
exceedances = sum_N_Excedencias;
Size_exceedances = size(exceedances);
N_exceedances_ascending_order = sort(exceedances(:,1)); %Sort N in ascending order
[val,pos]=intersect(exceedances,N_exceedances_ascending_order); %gives common val and its position between 'exceedances' and 'exceedances_ascending_order' arrays
Sigma_eq_Values = exceedances(pos,2); %gives 'Sigma_eq' values position in 'exceedances' array
exceedances_remastered = [N_exceedances_ascending_order Sigma_eq_Values]; %rewrite 'exceedances' array in ascending order of Sigma_eq Vs. Relative N from the original
%
Delta_N = real(N_exceedances_ascending_order);
Delta_N_remastered = sort(Delta_N,'descend');
Delta_S_ksi = real(Sigma_eq_Values);
Delta_S_ksi_remastered = sort(Delta_S_ksi);
Delta_S_MPa = ksi_to_MPa.*Delta_S_ksi;
exceedances_remastered_real = [Delta_N_remastered Delta_S_ksi_remastered Delta_S_MPa]; %Tabela de excedências considerando apenas valores reais de Sigma_eq e de N, fazendo (/\N=N) & (/\S=Sigma_eq)
%
%% **** OBS: Exportar para Excell [(Delta_N_remastered) (Delta_S_ksi_remastered)] **** %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\n')
fprintf('\t -------------------------------------------------------------')
fprintf('\n')
fprintf('\t\t\t .::. Tabulação Excedências .::.')
fprintf('\n')
fprintf('\n')
fprintf('\n')
Exceedances_Tab = array2table(exceedances_remastered,'VariableNames',{'N','S_eq_SWT'});
Exceedances_Tab_real = array2table(exceedances_remastered_real,'VariableNames',{'Delta_N','Delta_S_ksi','Delta_S_MPa'});
disp(Exceedances_Tab_real);
fprintf('\n')
%