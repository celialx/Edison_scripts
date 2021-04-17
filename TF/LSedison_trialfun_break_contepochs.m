function trl = LSedison_trialfun_break_contepochs(cfg)

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

lengthepochs=cfg.lengthepochs;

hdr   = ft_read_header(cfg.dataset);

T=cfg.T;
Beg_Task=(T.Epoch(find(T.Start(:,1)))-1)*30*hdr.Fs+1;
End_Task=T.Epoch(find(T.End(:,1)))*30*hdr.Fs;


% add this to the trl definition
begsample     = Beg_Task;
endsample     = End_Task;
offset        = 0;
NumEpochs=floor((endsample-begsample)/(lengthepochs*hdr.Fs));
trl = [round([begsample begsample+lengthepochs*hdr.Fs offset])];
for k=2:(NumEpochs)
    trl = [trl ; [round([begsample+(k-1)*lengthepochs*hdr.Fs begsample+(k+1-1)*lengthepochs*hdr.Fs offset])]];
end
