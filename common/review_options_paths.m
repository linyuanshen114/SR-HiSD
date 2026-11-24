function options = review_options_paths(options,base_dir)
%REVIEW_OPTIONS_PATHS Resolve paths before run() changes the current folder.
for name={'InputDir','OutputDir'}
    key=name{1};
    if ~isfield(options,key) || isempty(options.(key)); continue; end
    file=java.io.File(char(options.(key)));
    if ~file.isAbsolute(); file=java.io.File(base_dir,char(options.(key))); end
    options.(key)=char(file.getCanonicalPath());
end
end
