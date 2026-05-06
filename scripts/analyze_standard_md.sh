#!/usr/bin/env bash
set -eo pipefail

if [[ -f "/usr/local/gromacs/bin/GMXRC" ]]; then
    # shellcheck disable=SC1091
    source /usr/local/gromacs/bin/GMXRC
fi
set -u

# One-click standard MD analysis for protein + two ligands/DES components.
# Per replica it performs PBC centering, protein-fit trajectory generation,
# protein/ligand RMSD, RMSF, radius of gyration, protein-ligand H-bond analysis,
# frame extraction, and PNG plots.

: "${GMX_BIN:=gmx}"
: "${BASE_DIR:=md_run}"
: "${REPS:=auto}"
: "${SYSTEM_GROUP:=System}"
: "${CENTER_GROUP:=Protein}"
: "${FIT_GROUP:=Backbone}"
: "${PROTEIN_GROUP:=Protein}"
: "${LIG1_GROUP:=}"
: "${LIG2_GROUP:=}"
: "${OUT_BASE:=standard_md_analysis}"
: "${PDB_DT_PS:=1000}"
: "${RUN_HBOND:=yes}"
: "${HBOND_CMD:=auto}"
: "${PBC_MODE:=mol}"

usage() {
    cat << EOF
Usage:
  $0 [rep1 rep2 rep3]

Examples:
  $0
  $0 rep1 rep2 rep3
  LIG1_GROUP=CA1 LIG2_GROUP=LYS1 $0 rep1
  RUN_HBOND=no $0
  PBC_MODE=cluster $0

Outputs:
  ${OUT_BASE}/repX/xtc/md_center.xtc
  ${OUT_BASE}/repX/xtc/md_fit.xtc
  ${OUT_BASE}/repX/pdb/protein_ligands_dt${PDB_DT_PS}ps.pdb
  ${OUT_BASE}/repX/xvg/*.xvg
  ${OUT_BASE}/repX/png/*.png
  ${OUT_BASE}/repX/summary.tsv
  ${OUT_BASE}/combined_summary.tsv

Notes:
  - PBC_MODE=mol is conservative and keeps molecules whole after centering.
    For oligomers that appear to split across periodic boundaries, rerun with
    PBC_MODE=cluster and inspect the extracted PDB/movie before interpreting
    large protein RMSD/Rg values as real dissociation.
  - Protein RMSD is calculated after Backbone least-squares fitting.
  - Ligand pose RMSD is calculated on the protein-fitted trajectory with -fit none,
    so it reflects ligand drift relative to the protein-fitted frame.
  - Ligand self-fit RMSD is also reported for ligand internal conformational change.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

for cmd in "${GMX_BIN}" python3 awk sed grep; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] Missing command: ${cmd}" >&2
        exit 1
    }
done

if [[ $# -gt 0 ]]; then
    REPS="$*"
elif [[ "${REPS}" == "auto" ]]; then
    reps=()
    for d in "${BASE_DIR}"/rep*; do
        [[ -d "${d}" && -f "${d}/md.tpr" && -f "${d}/md.xtc" ]] || continue
        reps+=("$(basename "${d}")")
    done
    if [[ ${#reps[@]} -eq 0 ]]; then
        echo "[ERROR] No replica directories with md.tpr/md.xtc found under ${BASE_DIR}." >&2
        exit 1
    fi
    REPS="${reps[*]}"
fi

mkdir -p "${OUT_BASE}"
COMBINED="${OUT_BASE}/combined_summary.tsv"
printf "rep\tmetric\tgroup\tn_points\tx_start\tx_end\tmean\tsd\tmin\tmax\tlast\n" > "${COMBINED}"

infer_ligands_from_topology() {
    local top="${BASE_DIR}/topol.top"
    [[ -f "${top}" ]] || return 1
    mapfile -t top_ligs < <(awk '
        BEGIN { in_mol = 0 }
        /^\[ molecules \]/ { in_mol = 1; next }
        /^\[/ { if ($0 !~ /^\[ molecules \]/) in_mol = 0 }
        in_mol && $0 !~ /^;/ && NF >= 2 { print $1 }
    ' "${top}" | awk '
        $1 !~ /^(Protein|SOL|WAT|HOH|NA|CL|NA\+|CL-|K|K\+|MG|MG2\+|CA2\+|ZN|ZN2\+|Water_and_ions|System)$/ &&
        $1 !~ /^Protein_chain_/ &&
        $1 !~ /^topol_protein_/ { print }
    ')
    if [[ ${#top_ligs[@]} -ge 2 ]]; then
        LIG1_GROUP="${LIG1_GROUP:-${top_ligs[0]}}"
        LIG2_GROUP="${LIG2_GROUP:-${top_ligs[1]}}"
        echo "[INFO] Inferred ligand groups from topology: ${LIG1_GROUP}, ${LIG2_GROUP}"
        return 0
    fi
    return 1
}

infer_ligands_from_gro() {
    local gro="$1"
    [[ -f "${gro}" ]] || return 1
    mapfile -t gro_ligs < <(awk '
        NR > 2 && NF != 3 {
            resn = substr($0,6,5); gsub(/ /, "", resn)
            if (resn != "" &&
                resn !~ /^(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|HIE|HID|HIP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|ASH|GLH|LYN|CYM|CYX|ACE|NME)$/ &&
                resn !~ /^(SOL|WAT|HOH|NA|CL|NA\+|CL-|K|K\+|MG|MG2\+|CA2\+|ZN|ZN2\+)$/) seen[resn] = 1
        }
        END { for (r in seen) print r }
    ' "${gro}" | sort)
    if [[ ${#gro_ligs[@]} -ge 2 ]]; then
        LIG1_GROUP="${LIG1_GROUP:-${gro_ligs[0]}}"
        LIG2_GROUP="${LIG2_GROUP:-${gro_ligs[1]}}"
        echo "[INFO] Inferred ligand groups from ${gro}: ${LIG1_GROUP}, ${LIG2_GROUP}"
        return 0
    fi
    return 1
}

make_analysis_index() {
    local gro="$1"
    local out_idx="$2"
    local complex_name="Protein_${LIG1_GROUP}_${LIG2_GROUP}"
    printf "q\n" | "${GMX_BIN}" make_ndx -f "${gro}" -o "${out_idx}" >/dev/null 2>&1 || true

    awk -v lig1="${LIG1_GROUP}" -v lig2="${LIG2_GROUP}" -v group_name="${complex_name}" '
        function flush_group(name, arr, n,    i) {
            print "[ " name " ]"
            for (i = 1; i <= n; i++) {
                printf "%8d", arr[i]
                if (i % 15 == 0 || i == n) printf "\n"
            }
            print ""
        }
        BEGIN {
            split("ALA ARG ASN ASP CYS GLN GLU GLY HIS HIE HID HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASH GLH LYN CYM CYX ACE NME", aa, " ")
            for (i in aa) protein_res[aa[i]] = 1
        }
        NR <= 2 { next }
        NF == 3 { next }
        {
            resn = substr($0, 6, 5)
            atomn = substr($0, 16, 5) + 0
            gsub(/ /, "", resn)
            if (protein_res[resn] || resn == lig1 || resn == lig2) {
                atoms[++n] = atomn
            }
        }
        END {
            if (n > 0) flush_group(group_name, atoms, n)
        }
    ' "${gro}" >> "${out_idx}"
}

run_trjconv_center() {
    local tpr="$1" xtc="$2" out="$3" idx="$4"
    case "${PBC_MODE}" in
        mol)
            printf "%s\n%s\n" "${CENTER_GROUP}" "${SYSTEM_GROUP}" | \
                "${GMX_BIN}" trjconv -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -center -pbc mol -ur compact >/dev/null
            ;;
        nojump)
            printf "%s\n%s\n" "${CENTER_GROUP}" "${SYSTEM_GROUP}" | \
                "${GMX_BIN}" trjconv -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -center -pbc nojump -ur compact >/dev/null
            ;;
        cluster)
            printf "%s\n%s\n%s\n" "${CENTER_GROUP}" "${CENTER_GROUP}" "${SYSTEM_GROUP}" | \
                "${GMX_BIN}" trjconv -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -pbc cluster -center -ur compact >/dev/null
            ;;
        *)
            echo "[ERROR] Unsupported PBC_MODE=${PBC_MODE}; use mol, nojump, or cluster." >&2
            exit 1
            ;;
    esac
}

run_trjconv_fit() {
    local tpr="$1" xtc="$2" out="$3" idx="$4"
    printf "%s\n%s\n" "${FIT_GROUP}" "${SYSTEM_GROUP}" | \
        "${GMX_BIN}" trjconv -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -fit rot+trans >/dev/null
}

extract_pdb() {
    local tpr="$1" xtc="$2" out="$3" idx="$4"
    local complex_name="Protein_${LIG1_GROUP}_${LIG2_GROUP}"
    printf "%s\n" "${complex_name}" | \
        "${GMX_BIN}" trjconv -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -dt "${PDB_DT_PS}" >/dev/null || \
    printf "%s\n" "${SYSTEM_GROUP}" | \
        "${GMX_BIN}" trjconv -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -dt "${PDB_DT_PS}" >/dev/null
}

run_rmsd() {
    local tpr="$1" xtc="$2" idx="$3" fit_group="$4" rms_group="$5" out="$6" mode="$7"
    if [[ "${mode}" == "none" ]]; then
        if printf "%s\n%s\n" "${rms_group}" "${rms_group}" | \
            "${GMX_BIN}" rms -s "${tpr}" -f "${xtc}" -n "${idx}" -tu ns -fit none -o "${out}" >/dev/null 2>"${out%.xvg}.log"; then
            return 0
        fi
        echo "[WARN] ${GMX_BIN} rms -fit none failed for ${rms_group}; falling back to ligand self-fit RMSD. See ${out%.xvg}.log" >&2
    fi
    printf "%s\n%s\n" "${fit_group}" "${rms_group}" | \
        "${GMX_BIN}" rms -s "${tpr}" -f "${xtc}" -n "${idx}" -tu ns -o "${out}" >/dev/null
}

run_rmsf() {
    local tpr="$1" xtc="$2" idx="$3" group="$4" out="$5" avg_pdb="$6" bfac_pdb="$7" residue_flag="$8"
    if [[ "${residue_flag}" == "res" ]]; then
        printf "%s\n" "${group}" | \
            "${GMX_BIN}" rmsf -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -res -ox "${avg_pdb}" -oq "${bfac_pdb}" >/dev/null
    else
        printf "%s\n" "${group}" | \
            "${GMX_BIN}" rmsf -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" -ox "${avg_pdb}" -oq "${bfac_pdb}" >/dev/null
    fi
}

run_gyrate() {
    local tpr="$1" xtc="$2" idx="$3" group="$4" out="$5"
    printf "%s\n" "${group}" | "${GMX_BIN}" gyrate -s "${tpr}" -f "${xtc}" -n "${idx}" -o "${out}" >/dev/null
}

run_hbond_pair() {
    local tpr="$1" xtc="$2" idx="$3" lig="$4" out="$5" log="$6"
    [[ "${RUN_HBOND}" == "yes" ]] || return 0
    local cmd="${HBOND_CMD}"
    if [[ "${cmd}" == "auto" ]]; then
        if "${GMX_BIN}" hbond-legacy -h >/dev/null 2>&1; then
            cmd="hbond-legacy"
        else
            cmd="hbond"
        fi
    fi
    if [[ "${cmd}" == "hbond-legacy" ]]; then
        printf "%s\n%s\n" "${PROTEIN_GROUP}" "${lig}" | \
            "${GMX_BIN}" hbond-legacy -s "${tpr}" -f "${xtc}" -n "${idx}" -num "${out}" >"${log}" 2>&1 || \
            echo "[WARN] hbond-legacy failed for ${lig}; see ${log}" >&2
    else
        "${GMX_BIN}" hbond -s "${tpr}" -f "${xtc}" -n "${idx}" \
            -r "group \"${PROTEIN_GROUP}\"" -t "group \"${lig}\"" -num "${out}" >"${log}" 2>&1 || \
            echo "[WARN] hbond failed for ${lig}; see ${log}" >&2
    fi
}

plot_and_summarize_rep() {
    local rep="$1" xvg_dir="$2" png_dir="$3" out_root="$4"
    python3 - "${rep}" "${xvg_dir}" "${png_dir}" "${out_root}" "${COMBINED}" "${LIG1_GROUP}" "${LIG2_GROUP}" << 'PYPLOT'
import math, os, re, statistics as st, sys
rep, xvg_dir, png_dir, out_root, combined, lig1, lig2 = sys.argv[1:]
try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

def read_xvg(path):
    t, y = [], []
    if not os.path.exists(path):
        return t, y
    for line in open(path, encoding='utf-8', errors='replace'):
        s=line.strip()
        if not s or s[0] in '#@':
            continue
        p=re.split(r'\s+', s)
        if len(p) >= 2:
            t.append(float(p[0]))
            y.append(float(p[1]))
    return t, y

def as_ns(x):
    if x and max(x) > 1000:
        return [v / 1000.0 for v in x]
    return x

def split_repeated_index(x, y):
    segments = []
    sx, sy = [], []
    last = None
    for xi, yi in zip(x, y):
        if last is not None and xi <= last and sx:
            segments.append((sx, sy))
            sx, sy = [], []
        sx.append(xi)
        sy.append(yi)
        last = xi
    if sx:
        segments.append((sx, sy))
    return segments

def fmt(x):
    return 'NA' if isinstance(x, float) and math.isnan(x) else f'{x:.6g}'

def stats(metric, group, x, y):
    if not y:
        return None
    return [rep, metric, group, str(len(y)), fmt(x[0] if x else float('nan')), fmt(x[-1] if x else float('nan')), fmt(st.mean(y)), fmt(st.pstdev(y) if len(y)>1 else 0.0), fmt(min(y)), fmt(max(y)), fmt(y[-1])]

def save_plot(t, y, title, ylabel, out, xlabel='Time (ns)'):
    if plt is None or not y:
        return
    if xlabel == 'Time (ns)':
        t = as_ns(t)
    plt.figure(figsize=(8.2,4.8), dpi=150)
    plt.plot(t, y, lw=1.1)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out)
    plt.close()

def save_rmsf_plot(x, y, group, title, out, csv_out):
    if not y:
        return
    segments = split_repeated_index(x, y)
    with open(csv_out, 'w', encoding='utf-8') as f:
        f.write('rep,group,segment,index,rmsf_nm\n')
        for i, (sx, sy) in enumerate(segments, 1):
            for xi, yi in zip(sx, sy):
                f.write(f'{rep},{group},{i},{fmt(xi)},{fmt(yi)}\n')
    if plt is None:
        return
    plt.figure(figsize=(8.2,4.8), dpi=150)
    colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#56B4E9', '#E69F00']
    for i, (sx, sy) in enumerate(segments, 1):
        label = f'chain/segment {i}' if len(segments) > 1 else None
        plt.plot(sx, sy, lw=1.1, marker='o', markersize=2.4, color=colors[(i - 1) % len(colors)], label=label)
    plt.title(title)
    plt.xlabel('Residue/atom index')
    plt.ylabel('RMSF (nm)')
    plt.grid(alpha=0.25)
    if len(segments) > 1:
        plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(out)
    plt.close()

def diagnostic_rows(rows):
    by_key = {}
    for r in rows:
        by_key[(r[1], r[2])] = r
    out = []
    for lig in (lig1, lig2):
        pose = by_key.get(('rmsd_ligand_pose', lig))
        self_fit = by_key.get(('rmsd_ligand_self', lig))
        hbond = by_key.get(('hbond', lig))
        if pose:
            mean_pose = float(pose[6]); max_pose = float(pose[9]); last_pose = float(pose[10])
            if mean_pose >= 1.0 or max_pose >= 2.0 or last_pose >= 1.0:
                out.append([rep, lig, 'ligand_pose_rmsd_high', f'mean={mean_pose:.3g}; max={max_pose:.3g}; last={last_pose:.3g}', 'possible ligand drift/detachment relative to protein'])
        if pose and self_fit:
            mean_pose = float(pose[6]); mean_self = float(self_fit[6])
            if mean_pose >= 1.0 and mean_self <= 0.3:
                out.append([rep, lig, 'rigid_body_ligand_drift', f'pose_mean={mean_pose:.3g}; self_mean={mean_self:.3g}', 'ligand conformation is stable but its position relative to protein changes'])
        if hbond:
            last_hbond = float(hbond[10]); mean_hbond = float(hbond[6])
            if mean_hbond <= 0.5 or last_hbond == 0:
                out.append([rep, lig, 'hbond_loss_or_low_occupancy', f'mean={mean_hbond:.3g}; last={last_hbond:.3g}', 'check minimum distance/contact occupancy before declaring stable binding'])
    protein = by_key.get(('rmsd_protein', 'Protein'))
    rg = by_key.get(('rg_protein', 'Protein'))
    if protein and (float(protein[6]) >= 0.5 or float(protein[9]) >= 2.0):
        out.append([rep, 'Protein', 'protein_rmsd_high', f'mean={float(protein[6]):.3g}; max={float(protein[9]):.3g}; last={float(protein[10]):.3g}', 'check PBC treatment and oligomer stability'])
    if rg and (float(rg[8]) > 0 and float(rg[9]) / float(rg[8]) >= 1.5):
        out.append([rep, 'Protein', 'rg_large_excursion', f'min={float(rg[8]):.3g}; max={float(rg[9]):.3g}; last={float(rg[10]):.3g}', 'possible oligomer separation or PBC artifact'])
    return out

series = [
    ('rmsd_protein', 'Protein', 'rmsd_protein.xvg', 'RMSD (nm)', 'Protein backbone RMSD'),
    ('rg_protein', 'Protein', 'rg_protein.xvg', 'Rg (nm)', 'Protein radius of gyration'),
    ('rmsd_ligand_pose', lig1, f'rmsd_{lig1}_pose.xvg', 'RMSD (nm)', f'{lig1} pose RMSD after protein fit'),
    ('rmsd_ligand_pose', lig2, f'rmsd_{lig2}_pose.xvg', 'RMSD (nm)', f'{lig2} pose RMSD after protein fit'),
    ('rmsd_ligand_self', lig1, f'rmsd_{lig1}_self.xvg', 'RMSD (nm)', f'{lig1} self-fit RMSD'),
    ('rmsd_ligand_self', lig2, f'rmsd_{lig2}_self.xvg', 'RMSD (nm)', f'{lig2} self-fit RMSD'),
    ('hbond', lig1, f'hbond_protein_{lig1}.xvg', 'H-bond count', f'Protein-{lig1} hydrogen bonds'),
    ('hbond', lig2, f'hbond_protein_{lig2}.xvg', 'H-bond count', f'Protein-{lig2} hydrogen bonds'),
]
rows=[]
for metric, group, fname, ylabel, title in series:
    t, y = read_xvg(os.path.join(xvg_dir, fname))
    t_ns = as_ns(t)
    row = stats(metric, group, t_ns, y)
    if row:
        rows.append(row)
    save_plot(t, y, f'{rep} {title}', ylabel, os.path.join(png_dir, fname.replace('.xvg','.png')))

for group, fname in [('Protein','rmsf_protein_residue.xvg'), (lig1, f'rmsf_{lig1}_atom.xvg'), (lig2, f'rmsf_{lig2}_atom.xvg')]:
    x, y = read_xvg(os.path.join(xvg_dir, fname))
    row = stats('rmsf', group, x, y)
    if row:
        rows.append(row)
    save_rmsf_plot(x, y, group, f'{rep} {group} RMSF', os.path.join(png_dir, fname.replace('.xvg','.png')), os.path.join(out_root, fname.replace('.xvg','.csv')))

header='rep\tmetric\tgroup\tn_points\tx_start\tx_end\tmean\tsd\tmin\tmax\tlast\n'
summary=os.path.join(out_root, 'summary.tsv')
with open(summary, 'w', encoding='utf-8') as f:
    f.write(header)
    for r in rows:
        f.write('\t'.join(r)+'\n')
with open(combined, 'a', encoding='utf-8') as f:
    for r in rows:
        f.write('\t'.join(r)+'\n')
diag = diagnostic_rows(rows)
diag_path = os.path.join(out_root, 'diagnostics.tsv')
with open(diag_path, 'w', encoding='utf-8') as f:
    f.write('rep\tgroup\tflag\tvalues\tinterpretation\n')
    for r in diag:
        f.write('\t'.join(r)+'\n')
print(f'[OK] {summary}')
PYPLOT
}

analyze_rep() {
    local rep="$1"
    local rep_dir="${BASE_DIR}/${rep}"
    local tpr="${rep_dir}/md.tpr"
    local xtc="${rep_dir}/md.xtc"
    local gro="${rep_dir}/npt.gro"
    [[ -f "${tpr}" && -f "${xtc}" ]] || {
        echo "[WARN] Skipping ${rep}: missing ${tpr} or ${xtc}." >&2
        return 0
    }
    [[ -f "${gro}" ]] || gro="${tpr}"

    if [[ -z "${LIG1_GROUP}" || -z "${LIG2_GROUP}" ]]; then
        infer_ligands_from_topology || infer_ligands_from_gro "${gro}" || true
    fi
    if [[ -z "${LIG1_GROUP}" || -z "${LIG2_GROUP}" ]]; then
        echo "[ERROR] Could not infer ligand names. Set LIG1_GROUP and LIG2_GROUP." >&2
        exit 1
    fi

    local out_root="${OUT_BASE}/${rep}"
    local xvg_dir="${out_root}/xvg"
    local png_dir="${out_root}/png"
    local xtc_dir="${out_root}/xtc"
    local pdb_dir="${out_root}/pdb"
    mkdir -p "${xvg_dir}" "${png_dir}" "${xtc_dir}" "${pdb_dir}"

    local idx="${out_root}/analysis.ndx"
    make_analysis_index "${gro}" "${idx}"

    echo "[INFO] ${rep}: PBC centering and protein-fit trajectory"
    run_trjconv_center "${tpr}" "${xtc}" "${xtc_dir}/md_center.xtc" "${idx}"
    run_trjconv_fit "${tpr}" "${xtc_dir}/md_center.xtc" "${xtc_dir}/md_fit.xtc" "${idx}"
    extract_pdb "${tpr}" "${xtc_dir}/md_fit.xtc" "${pdb_dir}/protein_ligands_dt${PDB_DT_PS}ps.pdb" "${idx}"

    echo "[INFO] ${rep}: RMSD/RMSF/Rg/H-bond analysis for ${LIG1_GROUP}/${LIG2_GROUP}"
    run_rmsd "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${FIT_GROUP}" "${FIT_GROUP}" "${xvg_dir}/rmsd_protein.xvg" "normal"
    run_rmsd "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG1_GROUP}" "${LIG1_GROUP}" "${xvg_dir}/rmsd_${LIG1_GROUP}_pose.xvg" "none"
    run_rmsd "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG2_GROUP}" "${LIG2_GROUP}" "${xvg_dir}/rmsd_${LIG2_GROUP}_pose.xvg" "none"
    run_rmsd "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG1_GROUP}" "${LIG1_GROUP}" "${xvg_dir}/rmsd_${LIG1_GROUP}_self.xvg" "normal"
    run_rmsd "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG2_GROUP}" "${LIG2_GROUP}" "${xvg_dir}/rmsd_${LIG2_GROUP}_self.xvg" "normal"

    run_rmsf "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${PROTEIN_GROUP}" "${xvg_dir}/rmsf_protein_residue.xvg" "${pdb_dir}/protein_avg.pdb" "${pdb_dir}/protein_bfac.pdb" "res"
    run_rmsf "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG1_GROUP}" "${xvg_dir}/rmsf_${LIG1_GROUP}_atom.xvg" "${pdb_dir}/${LIG1_GROUP}_avg.pdb" "${pdb_dir}/${LIG1_GROUP}_bfac.pdb" "atom"
    run_rmsf "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG2_GROUP}" "${xvg_dir}/rmsf_${LIG2_GROUP}_atom.xvg" "${pdb_dir}/${LIG2_GROUP}_avg.pdb" "${pdb_dir}/${LIG2_GROUP}_bfac.pdb" "atom"

    run_gyrate "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${PROTEIN_GROUP}" "${xvg_dir}/rg_protein.xvg"
    run_hbond_pair "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG1_GROUP}" "${xvg_dir}/hbond_protein_${LIG1_GROUP}.xvg" "${out_root}/hbond_${LIG1_GROUP}.log"
    run_hbond_pair "${tpr}" "${xtc_dir}/md_fit.xtc" "${idx}" "${LIG2_GROUP}" "${xvg_dir}/hbond_protein_${LIG2_GROUP}.xvg" "${out_root}/hbond_${LIG2_GROUP}.log"

    plot_and_summarize_rep "${rep}" "${xvg_dir}" "${png_dir}" "${out_root}"
}

for rep in ${REPS}; do
    analyze_rep "${rep}"
done

echo "Done. Combined summary: ${COMBINED}"
