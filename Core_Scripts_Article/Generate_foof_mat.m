% Calculate sprectral power density and create mat for fooof toolbox

AllSubjects = {'01HC'; '02ZB'; '03JF'; '04BS'; '05MB'; '06NJ'; '07LB'; '08CB'; '09FN'; '10OL'; '11BD'; '12TG'; '14CD'...
    ;'15CV' ; '16LK'; '17ML'; '18LW'; '19HJ'; '20JJ'; '21IH'; '22HK';'24MK'; '25DS'; '26DB'; '27JS'; '28MD'; '29AN'...
    ; '30LM';'32DP'; '33NM'; '34AP'; '35NM'; '36DF'; '37SA'; '39NL'; '40EE'; '41RL'; '42AP'...
    ;'43AC'; '44LM'; '45AK'; '46SR'; '47OK'; '48KD'; '49LL'; '50CV'; '51EG'; '52NF'; '53EI'; ...
    '55BB' ; '56AB'; '57LL'; '58BP'; '59SC'; '60HD'; '61JF'; '62IB'; '63EG'; '64LC';...
    '67VF'; '68CY'; '69AM'; '70MA'; '71KN'; '72IG'; '73GF'; '74LW';'75CC'; '76LC';...
    '77JD'; '78BG'; '79CR'; '80MB';'81LB'; '82RD'; '83MB'; '84RD'; '85AA'; ...
    '87MD';'88SN'; '89AA'; '90ED'; '91MC'; '92GG'; '93RS'; '94CC'; '95FD'; '96CH';...
    '97AJ'; '98JB'; '99VB'; '100CP'; '101MS'; '102OB'; '103JZ'; '104MR'; '105OG';...
    '107AS'; '108JC'; '109JF'; '110SD'; '111MB'; '112GC'};

pathEEG = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\EEG\';
pathT = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\T\';
SavePath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\foof_aperiodic\data\';
PSDPath ='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\foof_aperiodic\power_spectrum\';

for i_sub=1:length(AllSubjects)
    
    sujet=AllSubjects{i_sub,1};
    
    
    FileName=strcat(pathT, filesep,'T_',sujet, '.mat');
    load(FileName);
    
    % Get start and end of each break
    StartNap = find(T.Start(:,1)==1);
    EndNap = find(T.End(:,1)==1);
    
    FileName=strcat(pathEEG,sujet,'_EEG.edf'); 
    [header,signalHeader,signalCell] = blockEdfLoad(FileName);
    sampRate = signalHeader.samples_in_record;
    EEG = pop_biosig(FileName,  'blockrange',[StartNap EndNap] );
    EEG = eeg_checkset(EEG);
    
    % For now, get O1 channel
    for i =1:length(EEG.chanlocs)
        
        if strcmp(EEG.chanlocs(i).labels, 'O1-A2')
            a = i;
        end
    end
    
    data = EEG.data(a,:)
    
    save([SavePath, 'data_ch_O1', sujet, '.mat'], 'data')
    
    % Calculate a power spectrum with Welch's method
    [psd, freqs] = pwelch(data, 500, [], [], sampRate);
    
    % Save the power spectra out to mat files
    save([PSDPath, 'PSD_', sujet, '.mat'], 'freqs', 'psd');
end