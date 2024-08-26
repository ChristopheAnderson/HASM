clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WEAK
Weak1 = xlsread('../Simulation1/Weak_Sph_negative_data_64_.xlsx');
Weak2 = xlsread('../Simulation1/Weak_Sph_symetric_data_64_.xlsx');
Weak3 = xlsread('../Simulation1/Weak_Sph_positive_data_64_.xlsx');

%% MODERATE
Moderate1 = xlsread('../Simulation1/moderate_Sph_negative_data_64_.xlsx');
Moderate2 = xlsread('../Simulation1/moderate_Sph_symetric_data_64_.xlsx');
Moderate3 = xlsread('../Simulation1/moderate_Sph_positive_data_64_.xlsx');

%% STRONG
Strong1 = xlsread('../Simulation1/strong_Sph_negative_data_64_.xlsx');
Strong2 = xlsread('../Simulation1/strong_Sph_symetric_data_64_.xlsx');
Strong3 = xlsread('../Simulation1/strong_Sph_positive_data_64_.xlsx');


% All data
ds_list = {Weak1,Weak2,Moderate1,Moderate2,Moderate3,Strong1,Strong2,Strong3};
clear('Weak1','Weak2','Moderate1','Moderate2','Moderate3','Strong1','Strong2','Strong3');
ds_names = {'WN_64_.xlsx','WS_64_.xlsx','WP_64_.xlsx','MN_64_.xlsx','MS_64_.xlsx','MP_64_.xlsx','SN_64_.xlsx','SS_64_.xlsx','SP_64_.xlsx'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0;
%para = xlsread('30_parameters.xlsx');


for q = 1:1    %18
    
    raw_data = ds_list{q};
    n=size(raw_data,1);
    % Specify the file name and sheet name
    filename = ds_names{q};
    sheetname = 'Sheet1';
    output = zeros(1, 2);
    outputs = cell(1000, 1); % initialize an empty list object
    
    
    bar=waitbar(1,'Program is runing , please wait');
    min=3;max=10;
 
    %for j = 3:10    %1002
    for j = min:max    %1002
        niveau=(j-min+1)/(max-min);
        waitbar( niveau , bar , strcat(num2str(niveau*100),'/100 Program is runing , please wait'));
        
        HD = raw_data(:,[1 2 j]);
        train = HD(1:ceil(n*0.7),1:3); % Train matrix
        test  = HD(ceil(n*0.7)+1:n,1:3); % Test matrix
        ch = train(:,1:2);  % matrix of coordinates for the locations where the attribute are known
        zh = train(:,3);    % Attribute
        zht = test(:,3);    % Test attribute
        ck = test(:,1:2);   % matrix of coordinates for the locations where the attribute are not known
        %range=para(:,3);
        dmax=12.72792;
        range=1.59099;
        covmodel='sphericalC';nhmax=length(ch);nsmax=0;order=NaN;cs=ck;covparam=[1 range];

        % Extract the column with the measurements
        measurements = raw_data(:,j);

        % Calculate the sample mean and standard deviation
        means = mean(measurements);
        standard_deviation = std(measurements);

        % Calculate the number of observations and the confidence level
        n = length(measurements);
        confidence_level = 0.95;

        % Calculate the margin of error
        t_value = tinv(confidence_level + (1 - confidence_level)/2, n-1);
        margin_of_error = t_value * standard_deviation / sqrt(n);

        % Calculate the lower and upper bounds of the interval data
        lower_bound = means - margin_of_error;
        upper_bound = means + margin_of_error;

        % Store the interval data in a matrix
        interval_data = [lower_bound upper_bound];

        a= repmat(lower_bound, length(ck), 1);
        b= repmat(upper_bound, length(ck), 1);

        %clear('n','para','range','confidence_level','HD','interval_data','mean','raw_data','standard_deviation','t_value','test','train','margin_of_error','lower_bound','upper_bound','measurements');
        %clc;

        zkBME=BMEintervalMode(ck,ch,cs,zh,a,b,covmodel,covparam,nhmax,nsmax,dmax,order);
        clear('a','b','ch','ck','cs','means','nhmax','nsmax','range','confidence_level','HD','interval_data','mean','standard_deviation','t_value','test','train','margin_of_error','lower_bound','upper_bound','measurements');
    %clc;
        BME_MAE = mean(abs(zkBME-zht));
        BME_RMSE = sqrt(mean(power((zkBME-zht),2)));
        output = [BME_MAE BME_RMSE];
        outputs{j} = output;
        
        pause(0.02);
        
    end;
    
    close(bar);
    
    % Convert the cell array to a single matrix
    result = cell2mat(outputs);

    % Define the new column names
    new_col_names = {'Biais', 'RMSE'};

    % Convert the matrix to a table
    A_table = array2table(result);

    % Rename the columns of the table
    A_table.Properties.VariableNames = new_col_names;

    % Write the matrix to an Excel file
    
    repertoire ='../Output_BME/BME_';
    nom_fichier=strcat(repertoire,filename);
    writetable(A_table, nom_fichier);
    
end;

