"""Original approved Table S2 renderer; pure functions only."""
def tex_escape(value: str) -> str:
    return value.replace('\\', '\\textbackslash{}').replace('&', '\\&').replace('%', '\\%').replace('_', '\\_').replace('#', '\\#')

def signed(value: str) -> str:
    number = float(value)
    return f'{number:+.3f}'

def classification(row: dict[str, str]) -> str:
    if row['shift_class'] == 'L20_enriched_supported_candidate':
        return 'supported L20 increase'
    if row['stable_contact_core'] == 'True':
        return 'stable contact core'
    if row['shift_class'].startswith('directional_but_uncertain'):
        return 'directional, unresolved'
    return 'no resolved shift'

def build_table(rows: list[dict[str, str]]) -> str:
    lines = ['\\begin{landscape}', '\\scriptsize', '\\setlength{\\tabcolsep}{3pt}', '\\begin{longtable}{@{}lllcS[table-format=1.3]S[table-format=1.3]S[table-format=+1.3,retain-explicit-plus=true]cS[table-format=+1.3,retain-explicit-plus=true]l@{}}', '\\caption{Trajectory-level L20-minus-L10 residue-contact estimates for supported shifts and stable contact cores in the exploratory candidate set.}', '\\label{tab:residue_contact_shifts}\\\\', '\\toprule', 'Enzyme & Residue & Contact class & $n_{10}/n_{20}$ & \\multicolumn{1}{c}{Mean L10} & \\multicolumn{1}{c}{Mean L20} & \\multicolumn{1}{c}{$\\Delta$} & 95\\% interval & \\multicolumn{1}{c}{Direction-adjusted $\\Delta$} & Classification \\\\', '\\midrule', '\\endfirsthead', '\\multicolumn{10}{l}{\\tablename\\ \\thetable\\ continued}\\\\', '\\toprule', 'Enzyme & Residue & Contact class & $n_{10}/n_{20}$ & \\multicolumn{1}{c}{Mean L10} & \\multicolumn{1}{c}{Mean L20} & \\multicolumn{1}{c}{$\\Delta$} & 95\\% interval & \\multicolumn{1}{c}{Direction-adjusted $\\Delta$} & Classification \\\\', '\\midrule', '\\endhead', '\\midrule', '\\multicolumn{10}{r}{Continued on next page}\\\\', '\\endfoot', '\\bottomrule', '\\endlastfoot']
    protein_names = {'CUT1': 'TfCut1'}
    for row in rows:
        protein = protein_names.get(row['protein_short'], row['protein_short'])
        interval = f'[{signed(row['delta_ci95_low'])}, {signed(row['delta_ci95_high'])}]'
        values = [protein, row['residue_label'], row['landmark_class'], f'{row['n_l10']}/{row['n_l20']}', f'{float(row['mean_contact_l10']):.3f}', f'{float(row['mean_contact_l20']):.3f}', signed(row['mean_delta_l20_minus_l10']), interval, signed(row['direction_adjusted_delta']), classification(row)]
        lines.append(' & '.join((tex_escape(value) for value in values)) + ' \\\\')
    lines.extend(['\\end{longtable}', '\\footnotesize', '$n_{10}$ and $n_{20}$ denote retained-trajectory counts for L10 and L20, respectively. $\\Delta$ is the L20-minus-L10 difference in mean residue contact fraction. The 29 candidates shown comprise eight with supported L20 increases and 21 classified as stable contact cores. Each included comparison has at least three retained trajectories at each PET length. Complete results for all 95 candidates are provided in Supplementary Data~S2. Values are based on retained trajectories over 20--100 ns using the 0.45 nm protein-heavy-atom/PET-heavy-atom contact definition. Candidate selection and classification criteria are given in the Supporting Methods subsection \\textit{Exploratory Residue-Contact Shift Analysis}. The 95\\% intervals are percentile bootstrap intervals based on 20,000 whole-trajectory resamples performed separately within the L10 and L20 groups. Direction-adjusted differences give equal weight to presentation directions represented at both lengths.', '\\end{landscape}', ''])
    return '\n'.join(lines)
