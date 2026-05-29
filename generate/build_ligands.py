# author: Benjamin Tingle 8/17/2020
import sys
import os
import subprocess
import tarfile
import io
import time
import signal
import math

try:
    from openbabel import openbabel
except ImportError:
    import openbabel

from rdkit import Chem
from rdkit.Chem import rdMolAlign

DOCKBASE='/'.join(__file__.split("/")[0:-3])
sys.path.append(DOCKBASE + "/ligand/mol2db2_py3_strain")
sys.path.append(DOCKBASE + "/ligand/strain")

from mol2 import Mol2
from hydrogens import count_hydrogens
from hierarchy import Hierarchy
from mol2db2 import mol2db2_quick
from Torsion_Strain import calc_strain
from TL_Functions import Mol2MolSupplier_noF
from rdkit_conformers import (
    CONF_BACKEND_PROP,
    CONF_ENERGY_PROP,
    CONF_FORCEFIELD_PROP,
    generate_conformations,
    generate_seed_conformation,
)


def replace_mol2_name(mol2_text, name):
    lines = mol2_text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>MOLECULE") and i + 1 < len(lines):
            lines[i + 1] = name
            break
    return "\n".join(lines) + "\n"


def rdkit_mol_from_smiles(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", name)
    return mol


def rdkit_mol_from_mol2(mol2_text):
    def from_openbabel_roundtrip():
        obmol = openbabel.OBMol()
        converter = openbabel.OBConversion()
        converter.SetInAndOutFormats("mol2", "sdf")
        if not converter.ReadString(obmol, mol2_text):
            raise ValueError("OpenBabel could not read mol2")
        sdf_text = converter.WriteString(obmol)
        roundtrip_mol = Chem.MolFromMolBlock(sdf_text, sanitize=False, removeHs=False)
        if roundtrip_mol is None:
            raise ValueError("RDKit could not parse OpenBabel SDF")
        Chem.SanitizeMol(roundtrip_mol)
        return roundtrip_mol

    mol = Chem.MolFromMol2Block(mol2_text, sanitize=False, removeHs=False)
    if mol is None:
        return from_openbabel_roundtrip()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # RDKit's mol2 aromaticity parser is stricter than OpenBabel's writer.
        # A mol2->SDF round trip through OpenBabel preserves atom order for these
        # AMSOL mol2 files and gives RDKit a more reliable connection table.
        mol = from_openbabel_roundtrip()
    return mol


def rdkit_mol_to_mol2_text(mol, conf_id=0, name=None):
    if name is not None:
        mol = Chem.Mol(mol)
        mol.SetProp("_Name", name)
    sdf = Chem.MolToMolBlock(mol, confId=conf_id) + "\n$$$$\n"
    obmol = openbabel.OBMol()
    converter = openbabel.OBConversion()
    converter.SetInAndOutFormats("sdf", "mol2")
    if not converter.ReadString(obmol, sdf):
        raise ValueError("OpenBabel could not read RDKit SDF block")
    mol2_text = converter.WriteString(obmol)
    if name is not None:
        mol2_text = replace_mol2_name(mol2_text, name)
    return mol2_text


def mol2_text_from_template_and_rdkit(template, rdkit_mol):
    if rdkit_mol.GetNumAtoms() != len(template.atomNum):
        raise ValueError(
            "RDKit/template atom-count mismatch: {} != {}".format(
                rdkit_mol.GetNumAtoms(), len(template.atomNum)
            )
        )

    mol2 = template.copy()
    mol2.atomXyz = []
    mol2.inputEnergy = []
    mol2.inputTotalStrain = []
    mol2.inputMaxStrain = []
    mol2.inputHydrogens = []

    for conf in rdkit_mol.GetConformers():
        xyz = []
        for atom_idx in range(rdkit_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            xyz.append((pos.x, pos.y, pos.z))
        mol2.atomXyz.append(xyz)
        if conf.HasProp(CONF_ENERGY_PROP):
            mol2.inputEnergy.append(conf.GetDoubleProp(CONF_ENERGY_PROP))
        else:
            mol2.inputEnergy.append(9999.99)
        mol2.inputTotalStrain.append(9999.99)
        mol2.inputMaxStrain.append(9999.99)
        mol2.inputHydrogens.append(0)

    mol2.xyzCount = len(mol2.atomXyz)
    mol2.origXyzCount = mol2.xyzCount
    output = io.StringIO()
    mol2.writeMol2File(output)
    return output.getvalue()


def largest_rigid_atom_map(template):
    helper = object.__new__(Hierarchy)
    helper._getRigidStructures(len(template.atomNum), template.atomBonds, False)
    groups = {}
    for atom_idx, rigid_structure in enumerate(helper.rigidStructures):
        groups.setdefault(rigid_structure, []).append(atom_idx)
    largest = max(groups.values(), key=len)
    return [(atom_idx, atom_idx) for atom_idx in largest]


def pin_largest_rigid_fragment(rdkit_mol, template):
    if rdkit_mol.GetNumConformers() < 2:
        return rdkit_mol

    atom_map = largest_rigid_atom_map(template)
    if len(atom_map) < 2:
        return rdkit_mol

    rdkit_mol = Chem.Mol(rdkit_mol)
    ref_id = rdkit_mol.GetConformer(0).GetId()
    ref_conf = rdkit_mol.GetConformer(ref_id)
    ref_positions = {
        atom_idx: ref_conf.GetAtomPosition(atom_idx)
        for atom_idx, _ref_idx in atom_map
    }

    for conf in list(rdkit_mol.GetConformers())[1:]:
        rdMolAlign.AlignMol(
            rdkit_mol,
            rdkit_mol,
            prbCid=conf.GetId(),
            refCid=ref_id,
            atomMap=atom_map,
        )
        for atom_idx, pos in ref_positions.items():
            conf.SetAtomPosition(atom_idx, pos)

    return rdkit_mol


def angle_degrees(conf, atom_a, atom_b, atom_c):
    pa = conf.GetAtomPosition(atom_a)
    pb = conf.GetAtomPosition(atom_b)
    pc = conf.GetAtomPosition(atom_c)
    v1 = (pa.x - pb.x, pa.y - pb.y, pa.z - pb.z)
    v2 = (pc.x - pb.x, pc.y - pb.y, pc.z - pb.z)
    dot = sum(x * y for x, y in zip(v1, v2))
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return 0.0
    cosang = max(-1.0, min(1.0, dot / (n1 * n2)))
    return math.degrees(math.acos(cosang))


def reorder_seed_for_amsol(mol):
    """Put a non-collinear bonded triplet first for AMSOL internal coordinates."""
    if mol.GetNumAtoms() < 4:
        return mol

    conf = mol.GetConformer()
    heavy = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1}
    for center in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(center)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        heavy_neighbors = [idx for idx in neighbors if idx in heavy]
        candidates = heavy_neighbors if len(heavy_neighbors) >= 2 else neighbors
        for i, atom_a in enumerate(candidates):
            for atom_c in candidates[i + 1:]:
                angle = angle_degrees(conf, atom_a, center, atom_c)
                if 20.0 <= angle <= 160.0:
                    first = [atom_a, center, atom_c]
                    order = first + [idx for idx in range(mol.GetNumAtoms()) if idx not in first]
                    return Chem.RenumberAtoms(mol, order)
    return mol


def write_seed_mol2(smiles, name, outfile):
    mol = rdkit_mol_from_smiles(smiles, name)
    seed = generate_seed_conformation(mol)
    seed = reorder_seed_for_amsol(seed)
    mol2_text = rdkit_mol_to_mol2_text(seed, seed.GetConformer().GetId(), name)
    with open(outfile, "w") as out:
        out.write(mol2_text)
    return seed


# read in protomers
protomers = sys.argv[1]
protonated_flat = []
with open(protomers) as prot_f:
    for line in prot_f:
        smiles, name, prot_id = line.split()
        protonated_flat.append((name + "." + str(prot_id), smiles, int(prot_id)))

# sort the protonated list for the next step
protonated_flat = sorted(protonated_flat, key=lambda x:x[0])

# when protomers get expanded, expanded protomers will not be assigned a unique id
# e.g if a new protomer was expanded from ZINC123.1, the new protomer will still be named ZINC123.1
# we fix this here
protomers_aug = open(os.path.dirname(protomers) + "/protomers_in_rdkit", 'w')
name_prev = None
protonated_flat_aug = []
for protomer in protonated_flat:
    name = protomer[0]
    if name_prev == name:
        protid = int(name.split('.')[-1])
        name = '.'.join(name.split('.')[:-1] + [str(protid + 1)])
        protomers_aug.write(protomer[1] + " " + name + "\n")
        protonated_flat_aug.append((name, protomer[1], protid + 1))
        name_prev = name
    else:
        protomers_aug.write(protomer[1] + " " + name + "\n")
        protonated_flat_aug.append(protomer)
        name_prev = name
        dup_protid_cnt = 0
protomers_aug.close()
protonated_flat = protonated_flat_aug

# initialize working directories
os.makedirs("3d", exist_ok=True)
os.makedirs("solv", exist_ok=True)

# create one RDKit 3D seed mol2 per protomer for AMSOL. Ensemble generation
# happens after AMSOL so the downstream mol2 keeps AMSOL atom order/charges.
protonated_success = []
rdkit_seed_mols = {}
total_protomers = len(protonated_flat)
t_seed_embed = 0.0

for mol2_index, p in enumerate(protonated_flat, start=1):
    mol2fn = str(mol2_index)
    name, smiles, _prot_id = p
    try:
        start = time.time()
        rdkit_seed_mols[mol2fn] = write_seed_mol2(smiles, name, os.path.join("3d", mol2fn))
        t_seed_embed += time.time() - start
    except Exception as exc:
        print("RDKit seed embedding failed for {}: {}".format(name, exc), file=sys.stderr)
        continue
    protonated_success.append((mol2fn, p))

protonated_flat = protonated_success

print("{} / {} protomers built successfully".format(len(protonated_flat), total_protomers))

# class for handling the conversion between Mol2 classes needed by various stages of the pipeline
class MultiMol2:
    def __init__(self, data, rdkit_mol=None):
        self.data = data
        self.dockFormat = Mol2(mol2text=[line+'\n' for line in data.split('\n')])
        self.rdkitFormat = Chem.Mol(rdkit_mol) if rdkit_mol is not None else rdkit_mol_from_mol2(data)

    @staticmethod
    def rdkit2str(rdkit_mol, template):
        return mol2_text_from_template_and_rdkit(template, rdkit_mol)

    @staticmethod
    def rdkit2dock(rdkit_mol, template):
        s = MultiMol2.rdkit2str(rdkit_mol, template).split('\n')
        return Mol2(mol2text=[line+'\n' for line in s])

    @staticmethod
    def rdkit2supplier(rdkit_mol, template):
        s = MultiMol2.rdkit2str(rdkit_mol, template)
        return Mol2MolSupplier_noF(s) # new function added to TL_Functions.py, creates a mol supplier object without reading from a file

    @staticmethod
    def rdkit2dock_and_supplier(rdkit_mol, template):
        s = MultiMol2.rdkit2str(rdkit_mol, template)
        return (
            Mol2(mol2text=[line+'\n' for line in s.split('\n')]),
            Mol2MolSupplier_noF(s),
        )

start_total = time.time()
start = time.time()

# read mol2 data into memory
mol2_data = []
solv_data = []
for mol2, protomer in protonated_flat:
    name, smiles, prot_id = protomer

    # calculate solvation as we go along
    if not os.path.isdir("solv/" + str(mol2)):
        os.makedirs("solv/" + str(mol2), exist_ok=True)
        mol2_fullname = '.'.join([name, 'mol2'])
        subprocess.call(["ln", "-svfn", os.getcwd() + "/3d/" + str(mol2), "solv/" + str(mol2) + "/" + mol2_fullname])
        os.chdir("solv/" + str(mol2))
        csh_exe = os.environ.get("CSH", "/bin/csh")
        subprocess.call([csh_exe, "-f", DOCKBASE + "/ligand/amsol/calc_solvation.csh", mol2_fullname])
        os.chdir("../..")

    if os.path.isfile("solv/" + str(mol2) + "/output.solv"):

        with open("solv/" + str(mol2) + "/output.solv") as solv:
            solv_text = solv.read()

        with open("solv/" + str(mol2) + "/output.mol2") as mol2_f:
            try:
                mol = MultiMol2(mol2_f.read(), rdkit_seed_mols.get(mol2))
            except Exception as exc:
                print("RDKit mol2 parsing failed for {}: {}".format(name, exc), file=sys.stderr)
                continue
            solv_data.append(solv_text)
            mol.idx = mol2
            mol.name = name
            mol.smiles = smiles
            mol.prot_id = prot_id
            mol.charge = int(float(solv_text.split('\n')[0].split()[2]))
            mol2_data.append(mol)

t_solvation = (time.time() - start)
print("{} / {} solvation successful".format(len(mol2_data), len(protonated_flat)))

# hard-code this bad boy in here
mp_range = ["M500","M400","M300","M200","M100","M000","P000","P010","P020","P030","P040","P050","P060","P070","P080","P090","P100","P110","P120","P130","P140","P150","P160","P170","P180","P190","P200","P210","P220","P230","P240","P250","P260","P270","P280","P290","P300","P310","P320","P330","P340","P350","P360","P370","P380","P390","P400","P410","P420","P430","P440","P450","P460","P470","P480","P490","P500","P600","P700","P800","P900"]
def get_zinc_directory_hash(zinc_id):
    if not zinc_id.startswith("ZINC"):
        return "."
    digits  = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    hp_b62  = zinc_id[4:6]
    si_b62  = zinc_id[6:16]
    hash_l1 = zinc_id[14:16]
    hash_l2 = zinc_id[12:14]
    h, p    = digits.index(hp_b62[0]), digits.index(hp_b62[1])
    tranche = "H{0:>02d}".format(h) + mp_range[p]
    return '/'.join([tranche[0:3], tranche, hash_l2, hash_l1])

def convert(data, inf, otf):
    mol = openbabel.OBMol()
    converter = openbabel.OBConversion()
    converter.SetInAndOutFormats(inf, otf)
    if not converter.ReadString(mol, data):
        raise ValueError("OpenBabel could not read {} input".format(inf))
    return converter.WriteString(mol)

def write_to_tarball(ball, data, name):
    tar = tarfile.TarInfo(name=name)
    tar.size = len(data)
    ball.addfile(tar, io.BytesIO(data))

tar_name = lambda x, *y:((x+'/') if x else '')+'.'.join([*y])

if os.path.isfile("output.tar.gz"):
    print("found ")
    subprocess.call(["mv", "output.tar.gz", "restart.tar.gz"])

stop = False

with tarfile.open("output.tar.gz", mode='w:gz') as output:

    # make sure to handle the SIGUSR1 interrupt so we properly close the tarfile
    def timeup_handler(signum, frame):
        global stop
        print("received SIGUSR1 : build_ligands.py")
        stop = True
        #output.close()
        #sys.exit(1)

    signal.signal(10, timeup_handler)
    signal.signal(2, timeup_handler)

    # allow for restartability
    # really need to rewrite this at some point...
    if os.path.isfile("restart.tar.gz"):
        with tarfile.open("restart.tar.gz", mode='r:gz') as restart:
            found_mol2s = []
            dest_mol2_map = {
                tar_name(get_zinc_directory_hash(mol.name), '.'.join([mol.name, str(mol.prot_id), chr(mol.charge+78)]), 'mol2') : mol for mol in mol2_data
            }
            for name in restart.getnames():
                with restart.extractfile(name) as memberfile:
                    data = memberfile.read()
                    write_to_tarball(output, data, name)
                # remove any mol2s from the worklist that already exist in the restart tarball
                if name.endswith('.mol2'):
                    found_mol2s.append(name)
            # set() fixes a weird bug where distinct protomers are generated but given id 0, causing identical mol2s in the output
            # update from the future: this is fixed at the protomer generation step now, this fix is no longer necessary
            found_mol2s = list(set(found_mol2s))
            for mol2 in [dest_mol2_map[m] for m in found_mol2s]:
                mol2_data.remove(mol2)

    t_strain_tot = t_db2_tot = t_conformer_tot = t_convert_tot = 0
    conformer_failures = 0
    strain_failures = 0
    db2_failures = 0
    conversion_failures = 0
    built_outputs = 0
    
    skip_conformers = True if os.getenv("SKIP_OMEGA") or os.getenv("SKIP_RDKIT_CONFORMERS") else False

    for i, mol in enumerate(mol2_data):
        zinc_hash = get_zinc_directory_hash(mol.name)
        mol_fullname = '.'.join([mol.name, chr(mol.charge+78)])
        print(zinc_hash + '/' + mol_fullname)

        # new wrapper function added to hydrogens.py, count_hydrogens
        h = count_hydrogens(mol.dockFormat)
        if not skip_conformers:
            start = time.time()
            try:
                _seed, ensemble = generate_conformations(mol.rdkitFormat, h)
                ensemble = pin_largest_rigid_fragment(ensemble, mol.dockFormat)
                db2ins = [MultiMol2.rdkit2dock_and_supplier(ensemble, mol.dockFormat)]
                print(
                    "{} conformers: {} ({})".format(
                        ensemble.GetConformer(0).GetProp(CONF_BACKEND_PROP),
                        ensemble.GetNumConformers(),
                        ensemble.GetConformer(0).GetProp(CONF_FORCEFIELD_PROP)
                    )
                )
            except Exception as exc:
                conformer_failures += 1
                print("conformer generation failed for {}: {}".format(mol.name, exc), file=sys.stderr)
                continue
            t_conformer_tot += (time.time() - start)
        else:
            db2ins = [MultiMol2.rdkit2dock_and_supplier(mol.rdkitFormat, mol.dockFormat)]

        db2_all_data = ""
        mol_failed = False
        for j, db2in in enumerate(db2ins):
            if stop == True:
                output.close()
                sys.exit()

            db2in_standard, db2in_strain = db2in
            # new wrapper function added to Torsion_Strain.py, calc_strain
            start = time.time()
            try:
                tE, pE = calc_strain(*db2in_strain)
            except Exception as exc:
                strain_failures += 1
                strain_count = len(db2in_strain[0]) if db2in_strain and db2in_strain[0] else getattr(db2in_standard, "xyzCount", 1)
                print(
                    "strain calculation failed for {}: {}; using zero-strain fallback".format(mol.name, exc),
                    file=sys.stderr,
                )
                tE = [0.0] * strain_count
                pE = [0.0] * strain_count
            t_strain_tot += (time.time() - start)
            # new utility function added to mol2.Mol2 class, addStrainInfo
            db2in_standard.addStrainInfo(tE, pE)
            # new wrapper function added to mol2db2.py, mol2db2_quick
            start = time.time()
            try:
                db2_data = mol2db2_quick(db2in_standard, solvfile="solv/" + str(mol.idx) + "/output.solv", clashfile=DOCKBASE + "/ligand/mol2db2/clashfile.txt")
            except Exception as exc:
                db2_failures += 1
                print("db2 generation failed for {}: {}".format(mol.name, exc), file=sys.stderr)
                mol_failed = True
                break
            # Remove the protomer ID inside the DB2 file itself
            db2_data = db2_data[:2] + f"{mol.name.split('.')[0]:16}" + db2_data[18:]
            db2_all_data += db2_data
            t_db2_tot += (time.time() - start)

            print(j+1, '/', len(db2ins))

        if mol_failed:
            continue

        start = time.time()
        try:
            pdbqt_data = convert(mol2_data[i].data, 'mol2', 'pdbqt')
            sdf_data   = convert(mol2_data[i].data, 'mol2', 'sdf')
        except Exception as exc:
            conversion_failures += 1
            print("archive conversion failed for {}: {}".format(mol.name, exc), file=sys.stderr)
            continue
        t_convert_tot += (time.time() - start)

        write_to_tarball(output, solv_data[i].encode('utf-8'),      name=tar_name(zinc_hash, mol_fullname, 'solv'))
        write_to_tarball(output, mol2_data[i].data.encode('utf-8'), name=tar_name(zinc_hash, mol_fullname, 'mol2'))
        write_to_tarball(output, sdf_data.encode('utf-8'),          name=tar_name(zinc_hash, mol_fullname, 'sdf'))
        write_to_tarball(output, pdbqt_data.encode('utf-8'),        name=tar_name(zinc_hash, mol_fullname, 'pdbqt'))
        write_to_tarball(output, db2_all_data.encode('utf-8'),      name=tar_name(zinc_hash, mol_fullname, 'db2.gz'))
        built_outputs += 1

    write_to_tarball(output, 'version={}'.format(os.environ.get("DOCK_VERSION")).encode('utf-8'), name=".dock_version")
    print("elapsed times:")
    print("seed embed: {}".format(t_seed_embed))
    print("solvation: {}".format(t_solvation))
    print("archive conversions: {}".format(t_convert_tot))
    print("strain:    {}".format(t_strain_tot))
    print("db2:       {}".format(t_db2_tot))
    print("confs:     {}".format(t_conformer_tot))
    print("total:     {}".format(time.time() - start_total))
    print("outputs built: {}".format(built_outputs))
    print("conformer failures: {}".format(conformer_failures))
    print("strain failures: {}".format(strain_failures))
    print("db2 failures: {}".format(db2_failures))
    print("archive conversion failures: {}".format(conversion_failures))
# that's all!
