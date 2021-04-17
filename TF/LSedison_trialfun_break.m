function trl = LSedison_trialfun_break(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr   = ft_read_header(cfg.dataset);

T=cfg.T;
Beg_Task=(T.Epoch(find(T.Start(:,1)))-1)*30*hdr.Fs+1;
End_Task=T.Epoch(find(T.End(:,1)))*30*hdr.Fs;

trl = [];

% add this to the trl definition
begsample     = Beg_Task;
endsample     = End_Task;
offset        = 0;
trl = [round([begsample endsample offset])];

