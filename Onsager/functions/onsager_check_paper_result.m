function paper = onsager_check_paper_result(route)
%ONSAGER_CHECK_PAPER_RESULT Check Figure 2's path and printed energy values.
% These constants are article quantities, not saved numerical reference fields.
paper_indices=[5 4 2 0];
paper_energies=[28.648747 25.115151 23.384271 4.608069];
printed_energy_tolerance=0.5e-6+1e-12;
assert(route.completed && isequal([route.nodes.index],paper_indices), ...
    'Onsager:PaperPath','The completed path must have indices 5 -> 4 -> 2 -> 0.');
energy_errors=abs([route.nodes.energy]-paper_energies);
assert(all(energy_errors<=printed_energy_tolerance), ...
    'Onsager:PaperEnergy','Completed state energies disagree with the article at six-decimal precision.');
paper=struct('passed',true,'indices',paper_indices,'printed_energies',paper_energies, ...
    'energy_errors',energy_errors,'energy_tolerance',printed_energy_tolerance);
end
