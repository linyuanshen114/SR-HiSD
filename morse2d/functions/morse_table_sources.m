function [original_file,npss_file] = morse_table_sources(ctx,allow_latest)
%MORSE_TABLE_SOURCES Find a pair from one run, with matching experiment settings.
if ~isempty(ctx.options.InputDir)
    roots={char(ctx.options.InputDir)};
elseif allow_latest
    roots={ctx.run_root};
    listing=dir(fullfile(ctx.task_dir,'results','*')); listing=listing([listing.isdir]);
    [~,order]=sort([listing.datenum],'descend');
    for k=order
        roots{end+1}=fullfile(listing(k).folder,listing(k).name); %#ok<AGROW>
    end
else
    roots={ctx.run_root};
end
for k=1:numel(roots)
    original_file=fullfile(roots{k},'original_methods','result.mat');
    npss_file=fullfile(roots{k},'npss','result.mat');
    if ~isfile(original_file) || ~isfile(npss_file), continue; end
    A=load(original_file,'run_info'); B=load(npss_file,'run_info');
    a=A.run_info; b=B.run_info;
    if strcmp(a.mode,ctx.mode) && strcmp(b.mode,ctx.mode) && ...
            isequal(a.quick,ctx.quick) && isequal(b.quick,ctx.quick) && ...
            isequal(a.Ns,b.Ns) && a.trials_per_N==b.trials_per_N
        return;
    end
end
error('Review:MorseTableInputs', ...
    'No matching original_methods + npss pair in one run. Run both with the same RunId or supply InputDir.');
end
