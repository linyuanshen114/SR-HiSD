function [V_lbsd, v_shared, alignment] = lb_align_initial_frame(V_original, v_sr, innp)
%BCC_ALIGN_INITIAL_DIRECTIONS Freeze SR v1 as the LBSD fourth initial column.
% Only the original first three near-zero columns are reorthogonalized. They
% are not replaced by exact translation derivatives. This is initialization,
% not a projection used during any later original-space LBSD iteration.

if ~isnumeric(V_original) || ~isreal(V_original) || ...
        ndims(V_original) ~= 2 || size(V_original,2) < 4 || ...
        ~isnumeric(v_sr) || ~isreal(v_sr) || ~isvector(v_sr) || ...
        size(V_original,1) ~= numel(v_sr) || ~isa(innp,'function_handle')
    error('bcc:InvalidInitialDirections', ...
        'Expected a real original-space frame with at least four columns and one SR vector.');
end
v_shared = v_sr(:);
if any(~isfinite(v_shared)) || any(any(~isfinite(V_original(:,1:4))))
    error('bcc:InvalidInitialDirections', 'Initial directions must be finite.');
end
% Match the existing benchmark normalize_rms convention, but reject zero
% directions rather than silently dividing them by the normalization floor.
v_shared = v_shared - mean(v_shared);
vnorm = direction_norm(v_shared, innp);
if vnorm <= 1e-14
    error('bcc:DegenerateSharedDirection', 'The zero-mean SR direction is degenerate.');
end
v_shared = v_shared / max(vnorm,1e-14);

old_v4 = V_original(:,4) - mean(V_original(:,4));
old_norm = direction_norm(old_v4, innp);
if old_norm <= 1e-14
    error('bcc:InvalidInitialDirections', 'The original fourth direction is degenerate.');
end
old_v4 = old_v4 / old_norm;

Q = zeros(numel(v_shared),3,'like',v_shared);
rank_tol = 1e-10;
for j = 1:3
    q = V_original(:,j) - mean(V_original(:,j));
    qnorm = direction_norm(q, innp);
    if qnorm <= 1e-14
        error('bcc:InitialFrameRankDeficient', 'Original near-zero column %d is degenerate.',j);
    end
    q = q / qnorm;
    % Keep the shared vector fixed. Project the other columns against it,
    % never the reverse, including on the second reorthogonalization pass.
    for pass = 1:2
        q = q - v_shared * innp(v_shared,q);
        q = q - Q(:,1:j-1) * innp(Q(:,1:j-1),q);
    end
    qnorm = direction_norm(q, innp);
    if qnorm <= rank_tol
        error('bcc:InitialFrameRankDeficient', ...
            'Near-zero columns lose rank after removing the shared SR direction (column %d).',j);
    end
    Q(:,j) = q / qnorm;
end
V_lbsd = [Q,v_shared];
entry = solver_entry_frame(V_lbsd,innp);
gram = innp(V_lbsd,V_lbsd) - eye(4);
entry_difference = direction_norm(entry(:,4)-v_shared,innp);
entry_overlap = real(innp(entry(:,4),v_shared));
if max(abs(gram(:))) > 1e-12 || entry_difference > 1e-12 || ...
        abs(entry_overlap-1) > 1e-12
    error('bcc:InitialAlignmentFailure', ...
        'The shared fourth direction does not survive original LBSD entry normalization.');
end

alignment = struct();
alignment.protocol = 'same_srkick_alignedv_blocktol1e8';
alignment.old_v4_abs_overlap = abs(real(innp(old_v4,v_shared)));
alignment.gram_error = max(abs(gram(:)));
alignment.entry_v4_rms_difference = entry_difference;
alignment.entry_v4_abs_overlap = abs(entry_overlap);
alignment.entry_v4_signed_overlap = entry_overlap;
alignment.first_three_source = ...
    'original_near_zero_columns_1_to_3_reorthogonalized_against_shared_sr_v1';
alignment.shared_direction_mean_error = abs(mean(v_shared));
alignment.near_zero_frame_mean_error = max(abs(mean(Q,1)));
alignment.rank_tolerance = rank_tol;
alignment.entry_alignment_tolerance = 1e-12;
alignment.scope = 'initial_directions_only_not_trajectory_equivalence';
end

function value = direction_norm(v,innp)
square = innp(v,v);
if ~isscalar(square) || ~isreal(square) || ~isfinite(square) || square < 0
    error('bcc:InvalidInitialInnerProduct', 'The supplied inner product is not positive and finite.');
end
value = sqrt(square);
end

function V = solver_entry_frame(V,innp)
% Deliberately match lb_lbsd entry GS, including its norm floor.
for j = 1:size(V,2)
    V(:,j) = V(:,j) - V(:,1:j-1) * innp(V(:,1:j-1),V(:,j));
    V(:,j) = V(:,j) / max(direction_norm(V(:,j),innp),1e-14);
end
end
