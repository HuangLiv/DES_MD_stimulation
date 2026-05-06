#!/usr/bin/env bash
#
# Protein + two-ligand GROMACS MD workflow
# ========================================
#
# This script builds and simulates a protein complex containing two ligands.
# It expects the MDP files to be stored in an mdp/ folder located in the same
# directory as this script:
#
#   run_md_prot_lig2_acpype.sh
#   mdp/ions.mdp
#   mdp/em.mdp
#   mdp/em1.mdp
#   mdp/nvt.mdp
#   mdp/npt.mdp
#   mdp/md.mdp
#
# Basic usage:
#
#   bash run_md_prot_lig2_acpype.sh
#   bash run_md_prot_lig2_acpype.sh 7f01.pdb ca.sdf lys.sdf
#
# Positional arguments:
#
#   1. protein PDB file
#   2. ligand 1 structure file, such as SDF or MOL2
#   3. ligand 2 structure file, such as SDF or MOL2
#
# If no positional arguments are provided, the default inputs are:
#
#   PROTEIN_PDB=7f01.pdb
#   LIG1_SDF=ca.sdf
#   LIG2_SDF=lys.sdf
#
# Ligand naming:
#
#   Ligand residue/molecule names are inferred automatically from file names.
#   For example:
#
#     ca.sdf  -> CA1
#     lys.sdf -> LYS1
#
#   The script rewrites the ACPYPE-generated ligand .gro/.itp files so that
#   the two ligands do not both remain named MOL. You can override the inferred
#   names explicitly:
#
#     LIG1_NAME=CIT LIG2_NAME=LYG bash run_md_prot_lig2_acpype.sh 7f01.pdb ca.sdf lys.sdf
#
# Common environment overrides:
#
#   LIG1_CHARGE=-3      ligand 1 net charge for ACPYPE/Antechamber
#   LIG2_CHARGE=1       ligand 2 net charge for ACPYPE/Antechamber
#   TEMP_K=310          simulation temperature
#   SALT_CONC=0.15      NaCl concentration in mol/L
#   N_REPLICAS=3        number of independent replicas
#   USE_GPU=yes         use GPU acceleration for equilibration/production
#   PRODUCTION_MODE=protein_hold
#                        main manuscript mode: restrained protein trimer, free ligands
#   GMX_BIN=gmx         set to gmx_mpi if using MPI GROMACS
#
# Restraint policy:
#
#   This system is a small 3 x 27-residue trimer with DES components bound on a
#   solvent-exposed surface rather than in a deep binding pocket. The default
#   production mode therefore keeps the protein trimer restrained and leaves the
#   ligands free, so the analysis focuses on dynamic surface H-bond/contact
#   exchange instead of artificial strong-pocket binding.
#
#     PRODUCTION_MODE=protein_hold             protein restrained, ligands free
#     PRODUCTION_MODE=protein_hold_ligand_weak protein restrained, ligands weakly restrained
#     PRODUCTION_MODE=free                     all components free; use as detachment control
#     PRODUCTION_MODE=pose_refine              protein + ligand restrained; pose refinement only
#     PRODUCTION_MODE=custom                   honor manual restraint variables
#
#   For manuscript-quality claims, compare protein_hold with free and/or
#   protein_hold_ligand_weak, and report contact occupancy/H-bond persistence
#   rather than claiming tight binding.
#
#     LIGAND_POSRES_IN_EQ=yes
#     PRODUCTION_MODE=protein_hold
#     LIGAND_POSRES_EQ_FC=1000
#     LIGAND_POSRES_MD_FC=25
#
# Important notes:
#
#   - The ligand input files must already have protonation states consistent
#     with the intended pH. ACPYPE -n sets the net charge, but it does not fix
#     an incorrectly protonated ligand structure.
#   - The ligand coordinates are assumed to already be in the correct binding
#     pose frame relative to the protein.
#   - Existing prep/, ligands/, and md_run/ outputs may be overwritten or mixed
#     with new outputs. Run in a clean working directory when starting a new
#     system.
#
set -eo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDP_DIR="${SCRIPT_DIR}/mdp"
if [[ -f "/usr/local/gromacs/bin/GMXRC" ]]; then
    # shellcheck disable=SC1091
    source /usr/local/gromacs/bin/GMXRC
fi
set -u
# ============ User-configurable parameters ============
: "${PROTEIN_PDB:=7f01.pdb}"
: "${LIG1_SDF:=ca.mol2}"
: "${LIG2_SDF:=lys.mol2}"
: "${LIG1_NAME:=}"
: "${LIG2_NAME:=}"

# Ligand net charges for ACPYPE/Antechamber at physiological pH.
# CA is modeled as citrate-like CA3-, and LYS as zwitterionic/protonated lysine-like LYS+.
# The SDF files must already contain the matching protonation states and hydrogens.
: "${LIG1_CHARGE:=-3}"
: "${LIG2_CHARGE:=1}"

: "${FF:=amber99sb-ildn}"
: "${WATER:=tip3p}"
# For elongated proteins, triclinic boxes with larger clearance reduce
# self-image artifacts around surface-bound ligands.
: "${BOX_TYPE:=triclinic}"
: "${BOX_DIST:=1.8}"
: "${SALT_CONC:=0.15}"
: "${TEMP_K:=310}"

# Replicate + reproducibility controls
: "${N_REPLICAS:=3}"
: "${BASE_SEED:=20260421}"
# deterministic: all replicas use the same seed
# ensemble: each replica uses BASE_SEED + replica_index
: "${REPLICA_SEED_MODE:=ensemble}"

# GPU acceleration controls
: "${USE_GPU:=yes}"
: "${NTMPI:=1}"
: "${NTOMP:=8}"

# Production-mode controls for weak surface-contact simulations.
: "${PROTEIN_POSRES_IN_MD:=no}"
: "${PRODUCTION_MODE:=protein_hold}"
: "${LIGAND_POSRES_IN_EQ:=yes}"
: "${LIGAND_POSRES_IN_MD:=no}"
: "${LIGAND_POSRES_EQ_FC:=1000}"
: "${LIGAND_POSRES_MD_FC:=25}"

# Energy-minimization controls. The production workflow uses steepest descent
# first, then regenerates em1.tpr from em.gro and runs conjugate gradients.
: "${EM_CG_EMTOL:=100}"
: "${EM_CG_EMSTEP:=0.0005}"
: "${EM_CG_NSTEPS:=50000}"
: "${EM_CG_FALLBACK_STEEP:=yes}"
: "${EM_FINAL_FALLBACK_STEEP:=yes}"
: "${EM_RESCUE_STEEP_EMTOL:=300}"
: "${EM_RESCUE_STEEP_EMSTEP:=0.001}"
: "${EM_RESCUE_STEEP_NSTEPS:=20000}"

# Use gmx_mpi if needed
: "${GMX_BIN:=gmx}"

case "${PRODUCTION_MODE}" in
    protein_hold)
        PROTEIN_POSRES_IN_MD=yes
        LIGAND_POSRES_IN_MD=no
        PRODUCTION_MODE_DESCRIPTION="protein trimer restrained; ligands free"
        ;;
    protein_hold_ligand_weak|stability)
        PROTEIN_POSRES_IN_MD=yes
        LIGAND_POSRES_IN_MD=yes
        PRODUCTION_MODE_DESCRIPTION="protein trimer restrained; ligands weakly restrained"
        ;;
    free)
        PROTEIN_POSRES_IN_MD=no
        LIGAND_POSRES_IN_MD=no
        PRODUCTION_MODE_DESCRIPTION="protein and ligands unrestrained; detachment/control mode"
        ;;
    pose_refine)
        PROTEIN_POSRES_IN_MD=yes
        LIGAND_POSRES_IN_MD=yes
        PRODUCTION_MODE_DESCRIPTION="protein and ligands restrained; restrained pose-refinement mode"
        ;;
    custom)
        PRODUCTION_MODE_DESCRIPTION="manual restraint variables honored"
        ;;
    *)
        echo "[ERROR] Unsupported PRODUCTION_MODE=${PRODUCTION_MODE}; use protein_hold, protein_hold_ligand_weak, free, pose_refine, or custom." >&2
        exit 1
        ;;
esac

require_yes_no() {
    local name="$1" value="$2"
    case "${value}" in
        yes|no) ;;
        *)
            echo "[ERROR] ${name} must be yes or no, got: ${value}" >&2
            exit 1
            ;;
    esac
}

require_yes_no PROTEIN_POSRES_IN_MD "${PROTEIN_POSRES_IN_MD}"
require_yes_no LIGAND_POSRES_IN_EQ "${LIGAND_POSRES_IN_EQ}"
require_yes_no LIGAND_POSRES_IN_MD "${LIGAND_POSRES_IN_MD}"
require_yes_no EM_CG_FALLBACK_STEEP "${EM_CG_FALLBACK_STEEP}"
require_yes_no EM_FINAL_FALLBACK_STEEP "${EM_FINAL_FALLBACK_STEEP}"
# =====================================================

usage() {
    cat << EOF
Usage:
  $0 [protein.pdb] [ligand1.sdf/mol2] [ligand2.sdf/mol2]

Environment overrides:
  LIG1_NAME, LIG2_NAME      Optional ligand residue/molecule names.
  LIG1_CHARGE, LIG2_CHARGE  Ligand net charges for ACPYPE.
  TEMP_K, SALT_CONC         Simulation temperature and salt concentration.
  PRODUCTION_MODE           protein_hold, protein_hold_ligand_weak, free, pose_refine, or custom.
  LIGAND_POSRES_EQ_FC       Ligand restraint force constant during equilibration.
  LIGAND_POSRES_MD_FC       Ligand restraint force constant during production.

Examples:
  $0
  PRODUCTION_MODE=protein_hold $0 7f01.pdb ca.mol2 lys.mol2
  PRODUCTION_MODE=free $0 7f01.pdb ca.mol2 lys.mol2
  PRODUCTION_MODE=protein_hold_ligand_weak $0 7f01.pdb ca.mol2 lys.mol2
  LIG1_NAME=CIT LIG2_NAME=LYG $0 7f01.pdb ca_pH74.mol2 lys_pH74.mol2
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -gt 0 ]]; then PROTEIN_PDB="$1"; fi
if [[ $# -gt 1 ]]; then LIG1_SDF="$2"; fi
if [[ $# -gt 2 ]]; then LIG2_SDF="$3"; fi
if [[ $# -gt 3 ]]; then
    echo "[ERROR] Too many positional arguments." >&2
    usage >&2
    exit 1
fi

make_ligand_name() {
    local input="$1"
    local fallback="$2"
    local stem
    stem="$(basename "${input}")"
    stem="${stem%.*}"
    stem="$(printf "%s" "${stem}" | tr '[:lower:]' '[:upper:]' | sed 's/[^A-Z0-9]//g')"
    if [[ -z "${stem}" ]]; then
        stem="${fallback}"
    fi
    printf "%.5s" "${stem}"
}

is_reserved_ligand_name() {
    case "$1" in
        ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|HIE|HID|HIP|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|ASH|GLH|LYN|CYM|CYX|ACE|NME|\
        SOL|WAT|HOH|NA|CL|K|CA|MG|ZN|FE|CU|MN|CO|NI|CD|HG|Protein|System|Water_and_ions)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

make_safe_ligand_name() {
    local name="$1"
    local fallback="$2"
    local safe
    safe="$(make_ligand_name "${name}" "${fallback}")"
    if is_reserved_ligand_name "${safe}"; then
        safe="$(printf "%.4s1" "${safe}")"
        echo "[WARN] Ligand name ${name} conflicts with a force-field/residue/ion name; using ${safe} instead." >&2
    fi
    printf "%s" "${safe}"
}

LIG1_NAME="${LIG1_NAME:-$(make_ligand_name "${LIG1_SDF}" "LIG1")}"
LIG2_NAME="${LIG2_NAME:-$(make_ligand_name "${LIG2_SDF}" "LIG2")}"
LIG1_NAME="$(make_safe_ligand_name "${LIG1_NAME}" "LIG1")"
LIG2_NAME="$(make_safe_ligand_name "${LIG2_NAME}" "LIG2")"
if [[ "${LIG1_NAME}" == "${LIG2_NAME}" ]]; then
    echo "[ERROR] Ligand names must be distinct. Got ${LIG1_NAME} for both ligands." >&2
    echo "        Set LIG1_NAME and LIG2_NAME explicitly." >&2
    exit 1
fi
TC_GROUP_NAME="Protein_${LIG1_NAME}_${LIG2_NAME}"

for cmd in "${GMX_BIN}" acpype obabel awk sed grep; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] Missing required command: ${cmd}" >&2
        exit 1
    }
done

run_mdrun() {
    local deffnm="$1"
    local stage="md"
    if [[ $# -ge 2 ]]; then
        stage="$2"
        shift 2 || true
    else
        shift || true
    fi

    if [[ "${USE_GPU}" == "yes" && "${stage}" != "em" ]]; then
        # GROMACS does not allow -reprod with GPU kernels.
        "${GMX_BIN}" mdrun \
            -deffnm "${deffnm}" \
            -ntmpi "${NTMPI}" \
            -ntomp "${NTOMP}" \
            -pin on \
            -nb gpu \
            -pme gpu \
            -bonded gpu \
            -update gpu \
            -v \
            "$@"
    else
        # EM/EM1: use conservative CPU-only settings to avoid rare CG segfaults.
        "${GMX_BIN}" mdrun \
            -deffnm "${deffnm}" \
            -ntmpi 1 \
            -ntomp 1 \
            -pin off \
            -nb cpu \
            -pme cpu \
            -bonded cpu \
            -update cpu \
            -v \
            "$@"
    fi
}

run_acpype() {
    # acpype/underlying scripts may source GMXRC and expect unset vars to be tolerated
    set +u
    acpype "$@"
    set -u
}

make_nvt_mdp_for_seed() {
    local seed="$1"
    local out_mdp="$2"
    local extra_define="${3:-}"

    awk -v seed="${seed}" -v extra_define="${extra_define}" -v temp_k="${TEMP_K}" -v tc_group="${TC_GROUP_NAME}" '
BEGIN {
    saw_define = 0
    saw_seed = 0
}
/^/ {
    raw = $0
    code = raw
    comment = ""
    cpos = index(raw, ";")
    if (cpos > 0) {
        code = substr(raw, 1, cpos - 1)
        comment = substr(raw, cpos)
    }
    sub(/[[:space:]]+$/, "", code)
}
/^[[:space:]]*define[[:space:]]*=/ {
    if (extra_define != "") {
        if (index(code, extra_define) == 0) {
            code = code " " extra_define
        }
    }
    if (comment != "") {
        print code " " comment
    } else {
        print code
    }
    saw_define = 1
    next
}
/^[[:space:]]*gen[-_]seed[[:space:]]*=/ {
    if (saw_seed == 0) {
        print "gen-seed                = " seed
        saw_seed = 1
    }
    next
}
/^[[:space:]]*ref[-_]t[[:space:]]*=/ {
    print "ref_t                   = " temp_k "   " temp_k
    next
}
/^[[:space:]]*tc[-_]grps[[:space:]]*=/ {
    print "tc-grps                 = " tc_group " Water_and_ions"
    next
}
/^[[:space:]]*gen[-_]temp[[:space:]]*=/ {
    print "gen_temp                = " temp_k
    next
}
{ print }
END {
    if (!saw_define && extra_define != "") {
        print "define                  = " extra_define
    }
    if (!saw_seed) {
        print "gen-seed                = " seed
    }
}
' "${MDP_DIR}/nvt.mdp" > "${out_mdp}"
}

make_mdp_with_extra_define() {
    local src_mdp="$1"
    local out_mdp="$2"
    local extra_define="${3:-}"

    awk -v extra_define="${extra_define}" -v temp_k="${TEMP_K}" -v tc_group="${TC_GROUP_NAME}" '
BEGIN {
    saw_define = 0
}
/^/ {
    raw = $0
    code = raw
    comment = ""
    cpos = index(raw, ";")
    if (cpos > 0) {
        code = substr(raw, 1, cpos - 1)
        comment = substr(raw, cpos)
    }
    sub(/[[:space:]]+$/, "", code)
}
/^[[:space:]]*define[[:space:]]*=/ {
    if (extra_define != "") {
        if (index(code, extra_define) == 0) {
            code = code " " extra_define
        }
    }
    if (comment != "") {
        print code " " comment
    } else {
        print code
    }
    saw_define = 1
    next
}
/^[[:space:]]*ref[-_]t[[:space:]]*=/ {
    print "ref_t                   = " temp_k "   " temp_k
    next
}
/^[[:space:]]*tc[-_]grps[[:space:]]*=/ {
    print "tc-grps                 = " tc_group " Water_and_ions"
    next
}
{ print }
END {
    if (!saw_define && extra_define != "") {
        print "define                  = " extra_define
    }
}
' "${src_mdp}" > "${out_mdp}"
}

make_em_override_mdp() {
    local src_mdp="$1"
    local out_mdp="$2"
    local integrator="$3"
    local emtol="$4"
    local emstep="$5"
    local nsteps="$6"

    awk -v integrator="${integrator}" -v emtol="${emtol}" -v emstep="${emstep}" -v nsteps="${nsteps}" '
BEGIN {
    saw_integrator = 0
    saw_emtol = 0
    saw_emstep = 0
    saw_nsteps = 0
}
/^[[:space:]]*integrator[[:space:]]*=/ { print "integrator              = " integrator; saw_integrator = 1; next }
/^[[:space:]]*emtol[[:space:]]*=/ { print "emtol                   = " emtol; saw_emtol = 1; next }
/^[[:space:]]*emstep[[:space:]]*=/ { print "emstep                  = " emstep; saw_emstep = 1; next }
/^[[:space:]]*nsteps[[:space:]]*=/ { print "nsteps                  = " nsteps; saw_nsteps = 1; next }
{ print }
END {
    if (!saw_integrator) print "integrator              = " integrator
    if (!saw_emtol) print "emtol                   = " emtol
    if (!saw_emstep) print "emstep                  = " emstep
    if (!saw_nsteps) print "nsteps                  = " nsteps
}
' "${src_mdp}" > "${out_mdp}"
}

run_em_steep_then_cg() {
    echo "[7/9] Energy minimization: steepest descent"
    "${GMX_BIN}" grompp \
        -f "${MDP_DIR}/em.mdp" \
        -c md_run/complex_solv_ions.gro \
        -p md_run/topol.top \
        -o md_run/em.tpr
    run_mdrun md_run/em em

    echo "[7.1/9] Secondary minimization: CG from first-step steep result (md_run/em.gro)"
    make_em_override_mdp "${MDP_DIR}/em1.mdp" "md_run/em1_cg.mdp" "cg" "${EM_CG_EMTOL}" "${EM_CG_EMSTEP}" "${EM_CG_NSTEPS}"
    "${GMX_BIN}" grompp \
        -f md_run/em1_cg.mdp \
        -c md_run/em.gro \
        -p md_run/topol.top \
        -o md_run/em1.tpr

    if run_mdrun md_run/em1 em; then
        return 0
    fi

    if [[ "${EM_CG_FALLBACK_STEEP}" != "yes" ]]; then
        echo "[ERROR] CG minimization failed and EM_CG_FALLBACK_STEEP=no." >&2
        return 1
    fi

    echo "[WARN] CG minimization failed, likely due to remaining bad contacts or water SETTLE failure." >&2
    echo "[WARN] Running rescue steepest-descent minimization, then retrying CG." >&2
    make_em_override_mdp "${MDP_DIR}/em1.mdp" "md_run/em1_rescue_steep.mdp" "steep" "${EM_RESCUE_STEEP_EMTOL}" "${EM_RESCUE_STEEP_EMSTEP}" "${EM_RESCUE_STEEP_NSTEPS}"
    "${GMX_BIN}" grompp \
        -f md_run/em1_rescue_steep.mdp \
        -c md_run/em.gro \
        -p md_run/topol.top \
        -o md_run/em1_rescue_steep.tpr
    run_mdrun md_run/em1_rescue_steep em

    "${GMX_BIN}" grompp \
        -f md_run/em1_cg.mdp \
        -c md_run/em1_rescue_steep.gro \
        -p md_run/topol.top \
        -o md_run/em1.tpr
    if run_mdrun md_run/em1 em; then
        return 0
    fi

    if [[ "${EM_FINAL_FALLBACK_STEEP}" != "yes" ]]; then
        echo "[ERROR] CG minimization failed after rescue steep and EM_FINAL_FALLBACK_STEEP=no." >&2
        return 1
    fi

    echo "[WARN] CG still failed after rescue steep." >&2
    echo "[WARN] Falling back to final steepest-descent minimization to EM_CG_EMTOL=${EM_CG_EMTOL}." >&2
    make_em_override_mdp "${MDP_DIR}/em1.mdp" "md_run/em1_final_steep.mdp" "steep" "${EM_CG_EMTOL}" "${EM_RESCUE_STEEP_EMSTEP}" "${EM_CG_NSTEPS}"
    "${GMX_BIN}" grompp \
        -f md_run/em1_final_steep.mdp \
        -c md_run/em1_rescue_steep.gro \
        -p md_run/topol.top \
        -o md_run/em1.tpr
    run_mdrun md_run/em1 em
}

for f in "${PROTEIN_PDB}" "${LIG1_SDF}" "${LIG2_SDF}"; do
    [[ -f "${f}" ]] || {
        echo "[ERROR] Required file not found: ${f}" >&2
        exit 1
    }
done

for f in ions.mdp em.mdp em1.mdp nvt.mdp npt.mdp md.mdp; do
    [[ -f "${MDP_DIR}/${f}" ]] || {
        echo "[ERROR] Required MDP file not found: ${MDP_DIR}/${f}" >&2
        exit 1
    }
done

cat << EOF
[INFO] Physiological-condition setup:
       Temperature: ${TEMP_K} K
       NaCl concentration: ${SALT_CONC} M
       Ligand 1: ${LIG1_SDF} -> ${LIG1_NAME}, net charge ${LIG1_CHARGE}
       Ligand 2: ${LIG2_SDF} -> ${LIG2_NAME}, net charge ${LIG2_CHARGE}
       Temperature-coupling group: ${TC_GROUP_NAME}
       Production mode: ${PRODUCTION_MODE} (${PRODUCTION_MODE_DESCRIPTION})
       Production restraints: protein=${PROTEIN_POSRES_IN_MD}, ligands=${LIGAND_POSRES_IN_MD}
       Ligand restraint force constants: EQ=${LIGAND_POSRES_EQ_FC}, MD=${LIGAND_POSRES_MD_FC} kJ mol^-1 nm^-2
[WARN] Ensure ${LIG1_SDF} and ${LIG2_SDF} use protonation states consistent with pH ~7.4.
       ACPYPE -n sets the net charge, but it does not by itself fix an incorrectly protonated SDF.
EOF

case "${PRODUCTION_MODE}" in
    protein_hold)
        cat << EOF
[INFO] PRODUCTION_MODE=protein_hold restrains the protein trimer during production
       and leaves ligands free. Interpret results as dynamic surface-contact
       simulations on a maintained trimer model, suitable for weak DES-protein
       H-bond/contact occupancy analysis.
EOF
        ;;
    protein_hold_ligand_weak|stability)
        cat << EOF
[INFO] PRODUCTION_MODE=${PRODUCTION_MODE} restrains the protein trimer and applies
       weak ligand restraints during production. Use this as a sensitivity/control
       run, not as sole evidence for spontaneous binding.
EOF
        ;;
    free)
        cat << EOF
[WARN] PRODUCTION_MODE=free leaves the protein and ligands unrestrained during production.
       For this small trimer and weak surface-binding system, ligand detachment and
       oligomer separation can be real outcomes. Use it as an exploratory/detachment
       control rather than the only production condition.
EOF
        ;;
    pose_refine)
        cat << EOF
[WARN] PRODUCTION_MODE=pose_refine restrains both protein and ligands during production.
       This is useful for restrained pose refinement, but it should not be used as
       the sole evidence for spontaneous weak surface association.
EOF
        ;;
    custom)
        cat << EOF
[INFO] PRODUCTION_MODE=custom uses manual restraint settings:
       protein=${PROTEIN_POSRES_IN_MD}, ligand_eq=${LIGAND_POSRES_IN_EQ}, ligand_md=${LIGAND_POSRES_IN_MD}.
EOF
        ;;
esac

mkdir -p prep ligands md_run

SEED_LOG="md_run/seeds.tsv"
echo -e "replicate\tseed\tmode" > "${SEED_LOG}"

echo "[1/9] Extracting standard protein residues from ${PROTEIN_PDB} (remove UNK/non-standard for pdb2gmx)"
awk '
BEGIN {
    split("ALA ARG ASN ASP CYS GLN GLU GLY HIS HIE HID HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASH GLH LYN CYM CYX ACE NME", a, " ")
    for (i in a) ok[a[i]] = 1
}
/^ATOM/ {
    res = substr($0,18,3)
    gsub(/ /, "", res)
    if (ok[res]) print $0
}
END { print "END" }
' "${PROTEIN_PDB}" > prep/protein_clean.pdb

ATOM_COUNT=$(grep -c '^ATOM' prep/protein_clean.pdb || true)
if [[ "${ATOM_COUNT}" -le 0 ]]; then
    echo "[ERROR] No standard protein ATOM records extracted from ${PROTEIN_PDB}" >&2
    exit 1
fi

echo "[2/9] Generating protein topology with ${FF}/${WATER}"
"${GMX_BIN}" pdb2gmx \
    -f prep/protein_clean.pdb \
    -o prep/protein_processed.gro \
    -p prep/topol_protein.top \
    -i prep/posre_protein.itp \
    -ff "${FF}" \
    -water "${WATER}" \
    -ignh

echo "[3/9] Parameterizing ligands by acpype (GAFF2 + AM1-BCC)"
run_acpype -i "${LIG1_SDF}" -b "${LIG1_NAME}" -a gaff2 -c bcc -n "${LIG1_CHARGE}" -o gmx
run_acpype -i "${LIG2_SDF}" -b "${LIG2_NAME}" -a gaff2 -c bcc -n "${LIG2_CHARGE}" -o gmx

[[ -d "${LIG1_NAME}.acpype" && -d "${LIG2_NAME}.acpype" ]] || {
    echo "[ERROR] acpype output folders missing." >&2
    exit 1
}

cp "${LIG1_NAME}.acpype/${LIG1_NAME}_GMX.itp" ligands/
cp "${LIG1_NAME}.acpype/posre_${LIG1_NAME}.itp" ligands/
cp "${LIG1_NAME}.acpype/${LIG1_NAME}_GMX.gro" ligands/
cp "${LIG2_NAME}.acpype/${LIG2_NAME}_GMX.itp" ligands/
cp "${LIG2_NAME}.acpype/posre_${LIG2_NAME}.itp" ligands/
cp "${LIG2_NAME}.acpype/${LIG2_NAME}_GMX.gro" ligands/

fix_ligand_names() {
    local lig="$1"
    local itp="ligands/${lig}_GMX.itp"
    local gro="ligands/${lig}_GMX.gro"
    local posre="ligands/posre_${lig}.itp"

    awk -v lig="${lig}" '
BEGIN { in_moleculetype = 0; in_atoms = 0; renamed = 0 }
/^\[ moleculetype \]/ { in_moleculetype = 1; in_atoms = 0; print; next }
/^\[ atoms \]/ { in_moleculetype = 0; in_atoms = 1; print; next }
/^\[/ {
    if ($0 !~ /^\[ moleculetype \]/) in_moleculetype = 0
    if ($0 !~ /^\[ atoms \]/) in_atoms = 0
}
in_moleculetype && renamed == 0 && $0 !~ /^;/ && NF >= 1 {
    printf "%-16s", lig
    for (i = 2; i <= NF; i++) printf " %s", $i
    printf "\n"
    renamed = 1
    next
}
in_atoms && $0 !~ /^;/ && NF >= 8 {
    $4 = lig
    print
    next
}
{ print }
' "${itp}" > "${itp}.tmp"
    mv "${itp}.tmp" "${itp}"

    awk -v lig="${lig}" '
NR <= 2 { print; next }
NF == 3 { print; next }
{
    resid = substr($0, 1, 5)
    rest = substr($0, 11)
    printf "%s%5s%s\n", resid, lig, rest
}
' "${gro}" > "${gro}.tmp"
    mv "${gro}.tmp" "${gro}"

    if [[ -f "${posre}" ]]; then
        sed -i.bak "s/\bMOL\b/${lig}/g" "${posre}"
        rm -f "${posre}.bak"
    fi
}

fix_ligand_names "${LIG1_NAME}"
fix_ligand_names "${LIG2_NAME}"

make_ligand_posres_variant() {
    local lig="$1"
    local fc="$2"
    local suffix="$3"
    local src="ligands/posre_${lig}.itp"
    local out="ligands/posre_${lig}_${suffix}.itp"

    [[ -f "${src}" ]] || {
        echo "[ERROR] Missing ligand position restraint file: ${src}" >&2
        exit 1
    }

    awk -v fc="${fc}" '
BEGIN { in_posres = 0 }
/^\[ position_restraints \]/ { in_posres = 1; print; next }
/^\[/ {
    if ($0 !~ /^\[ position_restraints \]/) in_posres = 0
}
in_posres && $0 !~ /^;/ && NF >= 5 {
    $3 = fc
    $4 = fc
    $5 = fc
    print
    next
}
{ print }
' "${src}" > "${out}"
}

make_ligand_posres_variant "${LIG1_NAME}" "${LIGAND_POSRES_EQ_FC}" "eq"
make_ligand_posres_variant "${LIG1_NAME}" "${LIGAND_POSRES_MD_FC}" "md"
make_ligand_posres_variant "${LIG2_NAME}" "${LIGAND_POSRES_EQ_FC}" "eq"
make_ligand_posres_variant "${LIG2_NAME}" "${LIGAND_POSRES_MD_FC}" "md"

# Merge atomtypes once, then strip [ atomtypes ] from each ligand topology.
awk '
BEGIN {
    print "; merged ligand atomtypes"
    print "[ atomtypes ]"
}
FNR == 1 { in_atomtypes = 0 }
/^\[ atomtypes \]/ { in_atomtypes = 1; next }
/^\[/ {
    if (in_atomtypes == 1) in_atomtypes = 0
}
in_atomtypes == 1 {
    if ($0 ~ /^;/ || NF == 0) next
    key = $1
    if (!(key in seen)) {
        seen[key] = 1
        print $0
    }
}
' "ligands/${LIG1_NAME}_GMX.itp" "ligands/${LIG2_NAME}_GMX.itp" > ligands/ligand_atomtypes.itp

awk '
BEGIN { skip = 0 }
/^\[ atomtypes \]/ { skip = 1; next }
/^\[/ {
    if (skip == 1) skip = 0
}
skip == 0 { print $0 }
' "ligands/${LIG1_NAME}_GMX.itp" > "ligands/${LIG1_NAME}_GMX_stripped.itp"

cat >> "ligands/${LIG1_NAME}_GMX_stripped.itp" << EOF

#ifdef POSRES_LIG_EQ
#include "posre_${LIG1_NAME}_eq.itp"
#endif

#ifdef POSRES_LIG_MD
#include "posre_${LIG1_NAME}_md.itp"
#endif
EOF

awk '
BEGIN { skip = 0 }
/^\[ atomtypes \]/ { skip = 1; next }
/^\[/ {
    if (skip == 1) skip = 0
}
skip == 0 { print $0 }
' "ligands/${LIG2_NAME}_GMX.itp" > "ligands/${LIG2_NAME}_GMX_stripped.itp"

cat >> "ligands/${LIG2_NAME}_GMX_stripped.itp" << EOF

#ifdef POSRES_LIG_EQ
#include "posre_${LIG2_NAME}_eq.itp"
#endif

#ifdef POSRES_LIG_MD
#include "posre_${LIG2_NAME}_md.itp"
#endif
EOF

echo "[4/9] Building initial complex coordinates"
"${GMX_BIN}" editconf -f prep/protein_processed.gro -o prep/protein_processed.pdb
"${GMX_BIN}" editconf -f "ligands/${LIG1_NAME}_GMX.gro" -o "ligands/${LIG1_NAME}_GMX.pdb"
"${GMX_BIN}" editconf -f "ligands/${LIG2_NAME}_GMX.gro" -o "ligands/${LIG2_NAME}_GMX.pdb"

{
    grep -E '^(ATOM|HETATM)' prep/protein_processed.pdb
    grep -E '^(ATOM|HETATM)' "ligands/${LIG1_NAME}_GMX.pdb"
    grep -E '^(ATOM|HETATM)' "ligands/${LIG2_NAME}_GMX.pdb"
    echo "TER"
    echo "END"
} > prep/complex_init.pdb

"${GMX_BIN}" editconf \
    -f prep/complex_init.pdb \
    -o prep/complex_box.gro \
    -bt "${BOX_TYPE}" \
    -d "${BOX_DIST}" \
    -c

echo "[5/9] Creating combined topology"
cp prep/topol_protein_Protein_chain_*.itp md_run/
cp prep/posre_protein_Protein_chain_*.itp md_run/
for chain_itp in md_run/topol_protein_Protein_chain_*.itp; do
    sed -i 's|"prep/posre_protein_|"posre_protein_|g' "${chain_itp}"
done

awk '
BEGIN { inserted = 0 }
/^#include "amber99sb-ildn\.ff\/forcefield\.itp"/ && inserted == 0 {
    print $0
    print "#include \"../ligands/ligand_atomtypes.itp\""
    print "#include \"../ligands/'"${LIG1_NAME}"'_GMX_stripped.itp\""
    print "#include \"../ligands/'"${LIG2_NAME}"'_GMX_stripped.itp\""
    print ""
    inserted = 1
    next
}
{ print }
' prep/topol_protein.top > md_run/topol.top

{
    printf "%-16s 1\n" "${LIG1_NAME}"
    printf "%-16s 1\n" "${LIG2_NAME}"
} >> md_run/topol.top

echo "[6/9] Solvation and ion addition"
"${GMX_BIN}" solvate \
    -cp prep/complex_box.gro \
    -cs spc216.gro \
    -o md_run/complex_solv.gro \
    -p md_run/topol.top

"${GMX_BIN}" grompp \
    -f "${MDP_DIR}/ions.mdp" \
    -c md_run/complex_solv.gro \
    -p md_run/topol.top \
    -o md_run/ions.tpr \
    -maxwarn 1

printf "SOL\n" | "${GMX_BIN}" genion \
    -s md_run/ions.tpr \
    -o md_run/complex_solv_ions.gro \
    -p md_run/topol.top \
    -pname NA \
    -nname CL \
    -neutral \
    -conc "${SALT_CONC}"

run_em_steep_then_cg

echo "[7b] Building protein-ligand index for temperature coupling"
awk -v lig1="${LIG1_NAME}" -v lig2="${LIG2_NAME}" -v tc_group="${TC_GROUP_NAME}" '
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
    system_atoms[++nsystem] = atomn
    gsub(/ /, "", resn)
    if (protein_res[resn] || resn == lig1 || resn == lig2) {
        pl[++npl] = atomn
    }
    if (resn == lig1) {
        l1[++nl1] = atomn
    }
    if (resn == lig2) {
        l2[++nl2] = atomn
    }
    if (resn == "SOL" || resn == "WAT" || resn == "HOH" ||
        resn == "NA" || resn == "CL" || resn == "K" ||
        resn == "MG" || resn == "ZN" || resn == "CA2") {
        wi[++nwi] = atomn
    }
}
END {
    if (npl == 0 || nwi == 0 || nl1 == 0 || nl2 == 0) {
        printf "[ERROR] Failed to build index from em1.gro: %s atoms=%d, Water_and_ions atoms=%d, %s atoms=%d, %s atoms=%d\n", tc_group, npl, nwi, lig1, nl1, lig2, nl2 > "/dev/stderr"
        exit 1
    }
    flush_group("System", system_atoms, nsystem)
    flush_group(lig1, l1, nl1)
    flush_group(lig2, l2, nl2)
    flush_group(tc_group, pl, npl)
    flush_group("Water_and_ions", wi, nwi)
}
' md_run/em1.gro > md_run/index.ndx

echo "[8/9] NVT + NPT equilibration"
eq_extra_define=""
if [[ "${LIGAND_POSRES_IN_EQ}" == "yes" ]]; then
    eq_extra_define="-DPOSRES_LIG_EQ"
fi

echo "[9/9] Replica equilibration + production MD (${N_REPLICAS} replicas)"
for rep in $(seq 1 "${N_REPLICAS}"); do
    rep_dir="md_run/rep${rep}"
    mkdir -p "${rep_dir}"

    if [[ "${REPLICA_SEED_MODE}" == "deterministic" ]]; then
        eq_seed="${BASE_SEED}"
    else
        eq_seed=$((BASE_SEED + rep))
    fi
    echo -e "rep${rep}\t${eq_seed}\t${REPLICA_SEED_MODE}" >> "${SEED_LOG}"

    echo "  [rep${rep}] NVT"
    make_nvt_mdp_for_seed "${eq_seed}" "${rep_dir}/nvt.mdp" "${eq_extra_define}"
    "${GMX_BIN}" grompp \
        -n md_run/index.ndx \
        -f "${rep_dir}/nvt.mdp" \
        -c md_run/em1.gro \
        -r md_run/em1.gro \
        -p md_run/topol.top \
        -o "${rep_dir}/nvt.tpr" \
        -maxwarn 1 \
        -v
    run_mdrun "${rep_dir}/nvt"

    echo "  [rep${rep}] NPT"
    make_mdp_with_extra_define "${MDP_DIR}/npt.mdp" "${rep_dir}/npt.mdp" "${eq_extra_define}"
    "${GMX_BIN}" grompp \
        -n md_run/index.ndx \
        -f "${rep_dir}/npt.mdp" \
        -c "${rep_dir}/nvt.gro" \
        -r "${rep_dir}/nvt.gro" \
        -t "${rep_dir}/nvt.cpt" \
        -p md_run/topol.top \
        -o "${rep_dir}/npt.tpr" \
        -maxwarn 1 \
        -v
    run_mdrun "${rep_dir}/npt"

    echo "  [rep${rep}] Production MD (${PRODUCTION_MODE}; protein restraints=${PROTEIN_POSRES_IN_MD}, ligand restraints=${LIGAND_POSRES_IN_MD})"
    md_defines=()
    if [[ "${PROTEIN_POSRES_IN_MD}" == "yes" ]]; then
        md_defines+=("-DPOSRES")
    fi
    if [[ "${LIGAND_POSRES_IN_MD}" == "yes" ]]; then
        md_defines+=("-DPOSRES_LIG_MD")
    fi
    md_extra_define="${md_defines[*]}"

    make_mdp_with_extra_define "${MDP_DIR}/md.mdp" "${rep_dir}/md.mdp" "${md_extra_define}"
    md_ref_opt=()
    if [[ "${PROTEIN_POSRES_IN_MD}" == "yes" || "${LIGAND_POSRES_IN_MD}" == "yes" ]]; then
        md_ref_opt=(-r "${rep_dir}/npt.gro")
    fi
    "${GMX_BIN}" grompp \
        -n md_run/index.ndx \
        -f "${rep_dir}/md.mdp" \
        -c "${rep_dir}/npt.gro" \
        "${md_ref_opt[@]}" \
        -t "${rep_dir}/npt.cpt" \
        -p md_run/topol.top \
        -o "${rep_dir}/md.tpr" \
        -maxwarn 1 \
        -v
    run_mdrun "${rep_dir}/md"
done

echo "Completed. Key outputs:"
echo "  md_run/rep1/md.xtc (and rep2/rep3)"
echo "  md_run/rep1/md.tpr (and rep2/rep3)"
echo "  md_run/rep1/md.edr (and rep2/rep3)"
echo "  md_run/rep1/npt.gro (and rep2/rep3)"
echo "  ${SEED_LOG}"
echo "  ${MODE_LOG}"

echo "NOTE: This script assumes ligand coordinates in ${LIG1_SDF}/${LIG2_SDF} are already in the correct binding pose frame relative to protein."
echo "NOTE: For this 3-chain/small-oligomer use case, use PRODUCTION_MODE=protein_hold for the main weak surface-contact analysis and compare against free/protein_hold_ligand_weak controls."
