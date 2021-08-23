function trl = LSedison_trialfun_arousals_perm(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr   = ft_read_header(cfg.dataset);

start_drop=(cfg.AllDrop-1)*hdr.Fs;

trl = [];
for i=1:length(start_drop)
    
        % add this to the trl definition
        begsample     = start_drop(i) - cfg.trialdef.prestim*hdr.Fs;
        endsample     = start_drop(i) + cfg.trialdef.poststim*hdr.Fs - 1;
        offset        = -cfg.trialdef.prestim*hdr.Fs;
        trl(end+1, :) = [round([begsample endsample offset])];
end
