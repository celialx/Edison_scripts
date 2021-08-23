function trl = LSedison_trialfun_arousals(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr   = ft_read_header(cfg.dataset);

start_Arousal=cfg.start_Arousal;
% ArousalFlag=T.Stade(1:end-1)~=0 & T.Stade(2:end)==0 & ~isnan(T.Stade(1:end-1));
% start_Arousal=T.TimeID(ArousalFlag)*hdr.Fs;
% Beg=T.TimeID(find(~isnan(T.Stade)))*hdr.Fs; Beg=Beg(1);
% diffArousals=diff([Beg ; start_Arousal])/hdr.Fs;
% start_Arousal(diffArousals<60)=[];

trl = [];
for i=1:length(start_Arousal)
    
        % add this to the trl definition
        begsample     = start_Arousal(i) - cfg.trialdef.prestim*hdr.Fs;
        endsample     = start_Arousal(i) + cfg.trialdef.poststim*hdr.Fs - 1;
        offset        = -cfg.trialdef.prestim*hdr.Fs;
        trl(end+1, :) = [round([begsample endsample offset])];
end
