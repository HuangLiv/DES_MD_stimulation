#!/usr/bin/env bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "/usr/local/gromacs/bin/GMXRC" ]]; then
    # shellcheck disable=SC1091
    source /usr/local/gromacs/bin/GMXRC
fi
set -u

# Dynamic weak-binding analysis for protein + two DES/ligand components.
# Outputs support weak, reversible, hydrogen-bond-level interaction analysis.

: "${GMX_BIN:=gmx}"
: "${BASE_DIR:=md_run}"
: "${REPS:=auto}"
: "${SYSTEM_GROUP:=Protein}"
: "${LIG1_GROUP:=}"
: "${LIG2_GROUP:=}"
: "${RUN_INTERACTION_ENERGY:=no}"
: "${SKIP_HBOND:=yes}"
: "${OUT_BASE:=dynamic_binding_analysis}"
: "${MDP_FILE:=${SCRIPT_DIR}/mdp/md.mdp}"

usage() {
    cat << EOF
Usage:
  $0 [rep1 rep2 rep3]

Examples:
  $0
  $0 rep1 rep2 rep3
  RUN_INTERACTION_ENERGY=yes $0 rep1
  LIG1_GROUP=CA1 LIG2_GROUP=LYS1 $0 rep1

Outputs:
  ${OUT_BASE}/repX/xvg/       raw XVG curves
  ${OUT_BASE}/repX/png/       plots
  ${OUT_BASE}/repX/summary.tsv
  ${OUT_BASE}/combined_summary.tsv

Notes:
  - Interaction energy is disabled by default because it requires clean rerun
    temperature-coupling and energy-group/index definitions.
  - If enabled, interaction energy is Coul-SR + LJ-SR from mdrun -rerun with energygrps.
    It is a short-range nonbonded interaction-energy proxy, not absolute ΔG.
  - For surface weak binding, interpret minimum distance and contact occupancy
    as the primary dynamic metrics; COM distance alone is often misleading.
    Hydrogen bonds are handled by the standard analysis script.
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
printf "rep\tligand\tn_frames\tt_start_ns\tt_end_ns\tcom_mean_nm\tcom_sd_nm\tcom_last_nm\tmindist_mean_nm\tmindist_sd_nm\tmindist_last_nm\tcontact_lt_0.35\tcontact_lt_0.45\tcontact_lt_0.60\thbond_mean\thbond_occupancy\tinteraction_mean_kjmol\tinteraction_sd_kjmol\tinteraction_last_kjmol\n" > "${COMBINED}"

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

index_has_group() {
    local idx_out="$1"
    local group="$2"
    echo "${idx_out}" | grep -Eq "[[:space:]][0-9]+[[:space:]]+${group}[[:space:]]*:"
}

make_analysis_index() {
    local tpr="$1"
    local out_root="$2"
    local idx_out
    idx_out=$(printf "q\n" | "${GMX_BIN}" make_ndx -f "${tpr}" -o "${out_root}/analysis_auto.ndx" 2>/dev/null || true)
    if [[ -z "${LIG1_GROUP}" || -z "${LIG2_GROUP}" ]]; then
        infer_ligands_from_topology || infer_ligands_from_gro "${REP_DIR}/npt.gro" || true
    fi
    if [[ -z "${LIG1_GROUP}" || -z "${LIG2_GROUP}" ]]; then
        echo "[ERROR] Could not infer two ligand groups. Set LIG1_GROUP and LIG2_GROUP." >&2
        return 1
    fi
    if ! index_has_group "${idx_out}" "${SYSTEM_GROUP}"; then
        echo "[WARN] ${SYSTEM_GROUP} was not listed by make_ndx output; selections may still resolve it from the TPR." >&2
    fi
    if ! index_has_group "${idx_out}" "${LIG1_GROUP}" || ! index_has_group "${idx_out}" "${LIG2_GROUP}"; then
        echo "[WARN] Ligand group(s) not listed by make_ndx; energy rerun may be skipped." >&2
    fi
}

write_energy_mdp() {
    local src="${MDP_FILE}"
    local out="$1"
    [[ -f "${src}" ]] || {
        echo "[WARN] MDP_FILE not found: ${src}; skipping interaction-energy MDP generation." >&2
        return 1
    }
    local lig="$2"
    awk -v sysgrp="${SYSTEM_GROUP}" -v liggrp="${lig}" '
        BEGIN { saw_energygrps = 0; saw_nsteps = 0; saw_continuation = 0; saw_genvel = 0 }
        /^[[:space:]]*nsteps[[:space:]]*=/ { print "nsteps                  = 0"; saw_nsteps = 1; next }
        /^[[:space:]]*continuation[[:space:]]*=/ { print "continuation            = yes"; saw_continuation = 1; next }
        /^[[:space:]]*gen[-_]vel[[:space:]]*=/ { print "gen_vel                 = no"; saw_genvel = 1; next }
        /^[[:space:]]*energygrps[[:space:]]*=/ { print "energygrps              = " sysgrp " " liggrp; saw_energygrps = 1; next }
        { print }
        END {
            if (!saw_nsteps) print "nsteps                  = 0"
            if (!saw_continuation) print "continuation            = yes"
            if (!saw_genvel) print "gen_vel                 = no"
            if (!saw_energygrps) print "energygrps              = " sysgrp " " liggrp
        }
    ' "${src}" > "${out}"
}

extract_interaction_energy() {
    local rep="$1" lig="$2" out_root="$3" xvg_dir="$4"
    local xtc="${BASE_DIR}/${rep}/md.xtc"
    local gro="${BASE_DIR}/${rep}/npt.gro"
    local cpt="${BASE_DIR}/${rep}/npt.cpt"
    local mdp="${out_root}/energy_${lig}.mdp"
    local rerun="${out_root}/energy_${lig}"
    local energy_xvg="${xvg_dir}/protein_${lig}_interaction_energy_raw.xvg"
    local sum_xvg="${xvg_dir}/protein_${lig}_interaction_energy_sum.xvg"

    [[ "${RUN_INTERACTION_ENERGY}" == "yes" ]] || return 0
    [[ -f "${BASE_DIR}/topol.top" && -f "${gro}" ]] || {
        echo "[WARN] Missing topol.top or ${gro}; skipping interaction-energy rerun for ${rep}/${lig}." >&2
        return 0
    }

    if ! write_energy_mdp "${mdp}" "${lig}"; then
        echo "[WARN] Could not prepare interaction-energy MDP for ${rep}/${lig}; set MDP_FILE or RUN_INTERACTION_ENERGY=no." >&2
        return 0
    fi
    grompp_cmd=("${GMX_BIN}" grompp -f "${mdp}" -c "${gro}" -p "${BASE_DIR}/topol.top" -n "${out_root}/analysis_auto.ndx" -o "${rerun}.tpr" -maxwarn 2)
    [[ -f "${cpt}" ]] && grompp_cmd+=( -t "${cpt}" )
    if ! "${grompp_cmd[@]}" >/dev/null 2>"${out_root}/energy_${lig}_grompp.log"; then
        echo "[WARN] grompp failed for interaction energy ${rep}/${lig}; see ${out_root}/energy_${lig}_grompp.log" >&2
        return 0
    fi
    if ! "${GMX_BIN}" mdrun -s "${rerun}.tpr" -rerun "${xtc}" -deffnm "${rerun}" -ntmpi 1 -ntomp 1 -nb cpu -pme cpu -bonded cpu -update cpu >/dev/null 2>"${out_root}/energy_${lig}_rerun.log"; then
        echo "[WARN] mdrun -rerun failed for interaction energy ${rep}/${lig}; see ${out_root}/energy_${lig}_rerun.log" >&2
        return 0
    fi

    if printf "Coul-SR:${SYSTEM_GROUP}-${lig}\nLJ-SR:${SYSTEM_GROUP}-${lig}\n0\n" | "${GMX_BIN}" energy -f "${rerun}.edr" -o "${energy_xvg}" >/dev/null 2>"${out_root}/energy_${lig}_extract.log"; then
        :
    elif printf "Coul-SR:${lig}-${SYSTEM_GROUP}\nLJ-SR:${lig}-${SYSTEM_GROUP}\n0\n" | "${GMX_BIN}" energy -f "${rerun}.edr" -o "${energy_xvg}" >/dev/null 2>>"${out_root}/energy_${lig}_extract.log"; then
        :
    else
        echo "[WARN] Could not extract Coul-SR/LJ-SR for ${rep}/${lig}; see ${out_root}/energy_${lig}_extract.log" >&2
        return 0
    fi

    python3 - "${energy_xvg}" "${sum_xvg}" << 'PYENERGY'
import re, sys
src, out = sys.argv[1:]
rows=[]
for line in open(src, encoding='utf-8', errors='replace'):
    s=line.strip()
    if not s or s[0] in '#@':
        continue
    p=re.split(r'\s+', s)
    if len(p) >= 3:
        rows.append((float(p[0]), sum(float(x) for x in p[1:])))
with open(out, 'w', encoding='utf-8') as f:
    f.write('@    title "Protein-ligand short-range interaction energy proxy"\n')
    f.write('@    xaxis  label "Time (ps)"\n')
    f.write('@    yaxis  label "Coul-SR + LJ-SR (kJ/mol)"\n')
    f.write('@TYPE xy\n')
    for t, y in rows:
        f.write(f'{t:12.3f} {y:12.5f}\n')
PYENERGY
}

analyze_rep() {
    local rep="$1"
    REP_DIR="${BASE_DIR}/${rep}"
    local tpr="${REP_DIR}/md.tpr"
    local xtc="${REP_DIR}/md.xtc"
    [[ -f "${tpr}" && -f "${xtc}" ]] || {
        echo "[WARN] Skipping ${rep}: missing md.tpr or md.xtc." >&2
        return 0
    }

    local out_root="${OUT_BASE}/${rep}"
    local xvg_dir="${out_root}/xvg"
    local png_dir="${out_root}/png"
    mkdir -p "${xvg_dir}" "${png_dir}"

    make_analysis_index "${tpr}" "${out_root}"
    echo "[INFO] ${rep}: analyzing ${LIG1_GROUP}/${LIG2_GROUP}"

    for lig in "${LIG1_GROUP}" "${LIG2_GROUP}"; do
        "${GMX_BIN}" distance -s "${tpr}" -f "${xtc}" -select "com of group \"${SYSTEM_GROUP}\" plus com of group \"${lig}\"" -oall "${xvg_dir}/protein_${lig}_comdist.xvg" >/dev/null
        "${GMX_BIN}" pairdist -s "${tpr}" -f "${xtc}" -ref "group \"${SYSTEM_GROUP}\"" -sel "group \"${lig}\"" -type min -o "${xvg_dir}/protein_${lig}_mindist.xvg" >/dev/null

        if [[ "${SKIP_HBOND}" != "yes" ]]; then
            if ! "${GMX_BIN}" hbond -s "${tpr}" -f "${xtc}" -r "group \"${SYSTEM_GROUP}\"" -t "group \"${lig}\"" -num "${xvg_dir}/protein_${lig}_hbond.xvg" >/dev/null 2>"${out_root}/hbond_${lig}.log"; then
                echo "[WARN] gmx hbond selection mode failed for ${rep}/${lig}; see ${out_root}/hbond_${lig}.log. Set SKIP_HBOND=yes to silence." >&2
            fi
        fi

        extract_interaction_energy "${rep}" "${lig}" "${out_root}" "${xvg_dir}"
    done

    python3 - "${rep}" "${xvg_dir}" "${png_dir}" "${out_root}" "${COMBINED}" "${LIG1_GROUP}" "${LIG2_GROUP}" << 'PYSUMMARY'
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
            t.append(float(p[0])/1000.0)
            y.append(float(p[1]))
    return t, y

def mean(xs): return st.mean(xs) if xs else float('nan')
def sd(xs): return st.pstdev(xs) if len(xs) > 1 else 0.0 if xs else float('nan')
def last(xs): return xs[-1] if xs else float('nan')
def fmt(x): return 'NA' if isinstance(x, float) and math.isnan(x) else f'{x:.5g}'

def plot_xy(t, y, title, ylabel, out):
    if plt is None or not t or not y:
        return
    plt.figure(figsize=(8.6, 4.8), dpi=140)
    plt.plot(t, y, lw=1.1)
    plt.title(title)
    plt.xlabel('Time (ns)')
    plt.ylabel(ylabel)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out)
    plt.close()

summary_path=os.path.join(out_root, 'summary.tsv')
with open(summary_path, 'w', encoding='utf-8') as sf, open(combined, 'a', encoding='utf-8') as cf:
    header='rep\tligand\tn_frames\tt_start_ns\tt_end_ns\tcom_mean_nm\tcom_sd_nm\tcom_last_nm\tmindist_mean_nm\tmindist_sd_nm\tmindist_last_nm\tcontact_lt_0.35\tcontact_lt_0.45\tcontact_lt_0.60\thbond_mean\thbond_occupancy\tinteraction_mean_kjmol\tinteraction_sd_kjmol\tinteraction_last_kjmol\n'
    sf.write(header)
    for lig in [lig1, lig2]:
        tc, com = read_xvg(os.path.join(xvg_dir, f'protein_{lig}_comdist.xvg'))
        tm, md = read_xvg(os.path.join(xvg_dir, f'protein_{lig}_mindist.xvg'))
        th, hb = read_xvg(os.path.join(xvg_dir, f'protein_{lig}_hbond.xvg'))
        te, en = read_xvg(os.path.join(xvg_dir, f'protein_{lig}_interaction_energy_sum.xvg'))
        n=len(md) or len(com)
        t0=(tm or tc or [float('nan')])[0]
        t1=(tm or tc or [float('nan')])[-1]
        occ35=sum(v < 0.35 for v in md)/len(md) if md else float('nan')
        occ45=sum(v < 0.45 for v in md)/len(md) if md else float('nan')
        occ60=sum(v < 0.60 for v in md)/len(md) if md else float('nan')
        hb_occ=sum(v > 0 for v in hb)/len(hb) if hb else float('nan')
        row=[rep, lig, str(n), fmt(t0), fmt(t1), fmt(mean(com)), fmt(sd(com)), fmt(last(com)), fmt(mean(md)), fmt(sd(md)), fmt(last(md)), fmt(occ35), fmt(occ45), fmt(occ60), fmt(mean(hb)), fmt(hb_occ), fmt(mean(en)), fmt(sd(en)), fmt(last(en))]
        line='\t'.join(row)+'\n'
        sf.write(line); cf.write(line)
        plot_xy(tc, com, f'{rep} Protein-{lig} COM distance', 'Distance (nm)', os.path.join(png_dir, f'protein_{lig}_comdist.png'))
        plot_xy(tm, md, f'{rep} Protein-{lig} minimum distance', 'Distance (nm)', os.path.join(png_dir, f'protein_{lig}_mindist.png'))
        plot_xy(th, hb, f'{rep} Protein-{lig} hydrogen bonds', 'H-bond count', os.path.join(png_dir, f'protein_{lig}_hbond.png'))
        plot_xy(te, en, f'{rep} Protein-{lig} interaction energy proxy', 'Coul-SR + LJ-SR (kJ/mol)', os.path.join(png_dir, f'protein_{lig}_interaction_energy.png'))
print(f'[OK] {summary_path}')
PYSUMMARY
}

for rep in ${REPS}; do
    analyze_rep "${rep}"
done

echo "Done. Combined summary: ${COMBINED}"
