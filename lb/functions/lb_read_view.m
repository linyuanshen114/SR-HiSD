function values=lb_read_view(data_dir,view_path,names)
%LB_READ_VIEW Read selected variables directly from the immutable HDF5 store.
% No compatibility files, source staging, catalogs or cache directories.
persistent loaded_file descriptor
manifest=fullfile(data_dir,'compact_views.json');
if isempty(loaded_file) || ~strcmp(loaded_file,manifest)
    descriptor=jsondecode(fileread(manifest)); loaded_file=manifest;
end
views=descriptor.views; if iscell(views), views=[views{:}]; end
index=find(strcmp({views.path},view_path));
assert(isscalar(index),'LBReview:View','Unknown or duplicate data view: %s',view_path);
variables=views(index).variables;
if nargin<3 || isempty(names), names=fieldnames(variables); end
if ischar(names), names={names}; end
values=struct();
for j=1:numel(names)
    assert(isfield(variables,names{j}),'LBReview:Variable','Missing %s in %s.',names{j},view_path);
    values.(names{j})=decode(variables.(names{j}),fullfile(data_dir,'arrays.h5'));
end
end

function value=decode(d,store)
switch d.kind
    case 'array'
        raw=double(h5read(store,['/objects/' d.sha256])); raw=raw(:);
        assert(strcmp(lb_numeric_sha(raw),d.sha256),'LBReview:DataHash','Numeric object hash mismatch.');
        value=reshape(cast(raw,d.class),shape(d.shape));
    case 'numeric'; value=reshape(cast(d.data(:),d.class),shape(d.shape));
    case 'char'; value=char(d.data);
    case 'struct'
        value=struct(); names=fieldnames(d.fields);
        for j=1:numel(names), value.(names{j})=decode(d.fields.(names{j}),store); end
    case 'empty_struct'; value=struct([]);
    case 'cell'
        items=d.data; if ~iscell(items), items=num2cell(items); end
        value=cell(shape(d.shape));
        for j=1:numel(items), value{j}=decode(items{j},store); end
    otherwise; error('LBReview:Descriptor','Unsupported data kind: %s',d.kind);
end
end

function dims=shape(value)
dims=double(value(:).'); if isempty(dims), dims=[1 1]; elseif numel(dims)==1, dims=[dims 1]; end
end
