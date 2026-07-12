"""
parse_sepr2.py  Pure-Python .NET BinaryFormatter (MS-NRBF) parser
Reads ISS_WT83.sepr2 and Ground_WT83.sepr2 without PatternLab software.
Extracts MyProtein objects and their SpectralCount fields, then runs
the VEN panel permutation test (identical methodology to all other analyses).

MS-NRBF specification: https://docs.microsoft.com/en-us/openspecs/windows_protocols/ms-nrbf
"""

import struct, io, re, os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)

VEN_PANEL = {
    "Myelination": ["MBP","MOG","PLP1","MAG","CNP","MOBP","ERMN"],
    "FastSignalling": ["SCN1A","KCNQ2","ANK3","NEFH","NEFM","NEFL","SNCG"],
    "SocialCircuit": ["OXTR","AVPR1A","HTR2A","DRD1","CHRM1","GABRB2"],
    "LayerVProjection": ["FEZF2","BCL11B","TBR1","SATB2","CUX1"],
    "MetabolicSupport": ["VDAC1","ATP2B2","SLC17A7","SNAP25","SYP","NRXN1"],
}
N_PERM, RNG_SEED = 10_000, 42
rng = np.random.default_rng(RNG_SEED)

#  MS-NRBF parser  
# PrimitiveTypeEnum
PT = {1:'bool',2:'byte',3:'char',4:'decimal',5:'double',6:'int16',
      7:'int32',8:'int64',9:'sbyte',10:'single',11:'timespan',
      12:'datetime',13:'uint16',14:'uint32',15:'uint64',16:'null',17:'string'}
PT_SIZE = {1:1,2:1,3:2,5:8,6:2,7:4,8:8,9:1,10:4,13:2,14:4,15:8}

# BinaryTypeEnum
BT = {0:'Primitive',1:'String',2:'Object',3:'SystemClass',4:'Class',
      5:'ObjectArray',6:'StringArray',7:'PrimitiveArray'}

# RecordTypeEnum
RT_HEADER = 0
RT_CLASS_WITH_ID = 1
RT_SYS_CLASS_NO_TYPES  = 2
RT_CLASS_NO_TYPES = 3
RT_SYS_CLASS = 4
RT_CLASS = 5 # ClassWithMembersAndTypes
RT_STRING = 6
RT_ARRAY = 7
RT_PRIM_TYPED = 8
RT_MEMBER_REF = 9
RT_NULL = 10
RT_END = 11
RT_LIBRARY = 12
RT_NULL_MULTI_256 = 13
RT_NULL_MULTI = 14
RT_ARRAY_PRIM = 15
RT_ARRAY_OBJ = 16
RT_ARRAY_STR = 17

class Reader:
    def __init__(self, data):
        self.d = data
        self.p = 0
    def byte(self):
        v = self.d[self.p]; self.p += 1; return v
    def i32(self):
        v = struct.unpack_from('<i', self.d, self.p)[0]; self.p += 4; return v
    def u32(self):
        v = struct.unpack_from('<I', self.d, self.p)[0]; self.p += 4; return v
    def i16(self):
        v = struct.unpack_from('<h', self.d, self.p)[0]; self.p += 2; return v
    def i64(self):
        v = struct.unpack_from('<q', self.d, self.p)[0]; self.p += 8; return v
    def f32(self):
        v = struct.unpack_from('<f', self.d, self.p)[0]; self.p += 4; return v
    def f64(self):
        v = struct.unpack_from('<d', self.d, self.p)[0]; self.p += 8; return v
    def bool(self):
        v = self.d[self.p]; self.p += 1; return bool(v)
    def lenpfx_str(self):
        n = 0; shift = 0
        while True:
            b = self.d[self.p]; self.p += 1
            n |= (b & 0x7F) << shift
            if not (b & 0x80): break
            shift += 7
        s = self.d[self.p:self.p+n].decode('utf-8', errors='replace')
        self.p += n
        return s
    def primitive(self, pt):
        if pt == 1: return self.bool()
        if pt == 2: v = self.d[self.p]; self.p += 1; return v
        if pt == 3: v = self.d[self.p:self.p+2].decode('utf-16-le'); self.p += 2; return v
        if pt == 5: return self.f64()
        if pt == 6: return self.i16()
        if pt == 7: return self.i32()
        if pt == 8: return self.i64()
        if pt == 9: v = struct.unpack_from('<b', self.d, self.p)[0]; self.p += 1; return v
        if pt == 10: return self.f32()
        if pt == 13: v = struct.unpack_from('<H', self.d, self.p)[0]; self.p += 2; return v
        if pt == 14: v = struct.unpack_from('<I', self.d, self.p)[0]; self.p += 4; return v
        if pt == 15: v = struct.unpack_from('<Q', self.d, self.p)[0]; self.p += 8; return v
        if pt == 17: return self.lenpfx_str()
        if pt == 11: v = self.i64(); return v # TimeSpan as ticks
        if pt == 12: v = self.i64(); return v # DateTime as ticks
        return None

class NRBFParser:
    def __init__(self, data):
        self.r = Reader(data)
        self.classes = {} # objectId -> ClassDef dict
        self.objects = {} # objectId -> value (dict for class instances)
        self.libraries = {} # libId -> name

    def parse(self):
        r = self.r
        # Stream header
        rt = r.byte()
        if rt != RT_HEADER:
            raise ValueError(f"Expected stream header (0), got {rt}")
        root_id = r.i32(); header_id = r.i32(); major = r.i32(); minor = r.i32()

        while r.p < len(r.d) - 1:
            try:
                rt = r.byte()
            except IndexError:
                break
            if rt == RT_END:
                break
            try:
                self._dispatch(rt)
            except Exception as e:
                # Stop gracefully on parse error
                break
        return self.objects

    def _dispatch(self, rt):
        r = self.r
        if rt == RT_LIBRARY:
            lid = r.i32(); name = r.lenpfx_str()
            self.libraries[lid] = name

        elif rt == RT_STRING:
            oid = r.i32(); s = r.lenpfx_str()
            self.objects[oid] = s

        elif rt == RT_NULL:
            pass  # nothing to read

        elif rt == RT_NULL_MULTI_256:
            count = r.byte()

        elif rt == RT_NULL_MULTI:
            count = r.i32()

        elif rt == RT_MEMBER_REF:
            ref_id = r.i32()
            # returns a lazy reference; we'll resolve later
            return ('ref', ref_id)

        elif rt == RT_PRIM_TYPED:
            pt = r.byte(); val = r.primitive(pt)

        elif rt in (RT_CLASS, RT_SYS_CLASS):
            # ClassWithMembersAndTypes
            oid = r.i32()
            cls_name = r.lenpfx_str()
            n_members = r.i32()
            member_names = [r.lenpfx_str() for _ in range(n_members)]
            bin_types = [r.byte() for _ in range(n_members)]
            # Additional info per BinaryType
            add_info = []
            for bt in bin_types:
                if bt == 0: add_info.append(r.byte()) # PrimitiveTypeEnum
                elif bt == 3: add_info.append(r.lenpfx_str())  # SystemClass name
                elif bt == 4: add_info.append((r.lenpfx_str(), r.i32()))  # ClassName + LibId
                elif bt == 7: add_info.append(r.byte()) # PrimitiveTypeEnum for array
                else: add_info.append(None)
            if rt == RT_CLASS:
                lib_id = r.i32()
            else:
                lib_id = None
            cdef = {'name': cls_name, 'member_names': member_names,
                    'bin_types': bin_types, 'add_info': add_info}
            self.classes[oid] = cdef
            obj = self._read_members(cdef)
            obj['__class__'] = cls_name
            self.objects[oid] = obj

        elif rt == RT_CLASS_WITH_ID:
            # ClassWithId instance of a previously-defined class
            oid = r.i32()
            ref_id = r.i32()
            cdef = self.classes.get(ref_id)
            if cdef is None:
                # Unknown class, can't parse, skip to next record
                return
            obj = self._read_members(cdef)
            obj['__class__'] = cdef['name']
            self.classes[oid] = cdef   # same class def, new instance
            self.objects[oid] = obj

        elif rt in (RT_SYS_CLASS_NO_TYPES, RT_CLASS_NO_TYPES):
            oid = r.i32()
            cls_name = r.lenpfx_str()
            n_members = r.i32()
            member_names = [r.lenpfx_str() for _ in range(n_members)]
            if rt == RT_CLASS_NO_TYPES:
                lib_id = r.i32()
            # No type info, read nothing for members (they follow as subsequent records)

        elif rt == RT_ARRAY_PRIM:
            oid = r.i32()
            length = r.i32()
            pt = r.byte()
            arr = [r.primitive(pt) for _ in range(length)]
            self.objects[oid] = arr

        elif rt == RT_ARRAY_STR:
            oid = r.i32()
            length = r.i32()
            arr = []
            for _ in range(length):
                rt2 = r.byte()
                if rt2 == RT_STRING:
                    sid = r.i32(); s = r.lenpfx_str()
                    self.objects[sid] = s; arr.append(s)
                elif rt2 == RT_NULL:
                    arr.append(None)
                elif rt2 == RT_MEMBER_REF:
                    ref_id = r.i32(); arr.append(('ref', ref_id))
                else:
                    arr.append(None)
            self.objects[oid] = arr

        elif rt == RT_ARRAY_OBJ:
            oid = r.i32()
            length = r.i32()
            arr = []
            for _ in range(length):
                rt2 = r.byte()
                val = self._dispatch(rt2)
                arr.append(val)
            self.objects[oid] = arr

        elif rt == RT_ARRAY:
            # BinaryArray generic array
            oid = r.i32()
            arr_type = r.byte() # ArrayTypeEnum
            rank = r.i32()
            lengths = [r.i32() for _ in range(rank)]
            # lower bounds for some types
            if arr_type in (3, 4, 5):  # Jagged, Rectangular, JaggedOffset
                lower = [r.i32() for _ in range(rank)]
            bt = r.byte()
            add = None
            if bt == 0: add = r.byte()
            elif bt == 3: add = r.lenpfx_str()
            elif bt == 4: add = (r.lenpfx_str(), r.i32())
            elif bt == 7: add = r.byte()
            total = 1
            for l in lengths: total *= l
            arr = []
            for _ in range(total):
                rt2 = r.byte()
                val = self._dispatch(rt2)
                arr.append(val)
            self.objects[oid] = arr

    def _read_member_value(self, bt, add):
        """Read one member value given its BinaryType and additional info."""
        r = self.r
        if bt == 0:   # Primitive
            return r.primitive(add)
        elif bt == 1: # String
            rt2 = r.byte()
            if rt2 == RT_STRING:
                sid = r.i32(); s = r.lenpfx_str()
                self.objects[sid] = s; return s
            elif rt2 == RT_NULL: return None
            elif rt2 == RT_MEMBER_REF: return ('ref', r.i32())
            else: return None
        elif bt in (2, 3, 4, 5, 6, 7):  # Object / Class / Arrays
            rt2 = r.byte()
            if rt2 == RT_NULL: return None
            elif rt2 == RT_MEMBER_REF: return ('ref', r.i32())
            elif rt2 == RT_STRING:
                sid = r.i32(); s = r.lenpfx_str()
                self.objects[sid] = s; return s
            else:
                self._dispatch(rt2)
                return None
        return None

    def _read_members(self, cdef):
        obj = {}
        for name, bt, add in zip(cdef['member_names'], cdef['bin_types'], cdef['add_info']):
            val = self._read_member_value(bt, add)
            # Strip k__BackingField wrapper
            clean = name.replace('<','').replace('>k__BackingField','')
            obj[clean] = val
        return obj


def extract_proteins_from_sepr2(filepath):
    """
    Parse .sepr2 and return dict: gene_symbol -> spectral_count (int).
    Falls back to regex scan if NRBF parse fails.
    """
    print(f"\nParsing: {filepath}  ({os.path.getsize(filepath)//1024//1024} MB)")
    with open(filepath, 'rb') as f:
        data = f.read()

    #Primary: NRBF parse 
    proteins = {}
    try:
        parser = NRBFParser(data)
        objs = parser.parse()
        print(f"  Objects parsed: {len(objs):,}")

        # Find MyProtein objects
        for oid, obj in objs.items():
            if not isinstance(obj, dict): continue
            cls = obj.get('__class__','')
            if 'MyProtein' not in cls: continue

            # Extract locus (protein identifier)
            locus = (obj.get('Locus') or obj.get('locus') or
                     obj.get('ProteinName') or obj.get('proteinName') or '')
            if isinstance(locus, tuple): locus = ''  # unresolved ref

            # Parse gene from UniProt locus string  e.g. sp|Q9UQB8|NEFM_HUMAN
            gene = ''
            m = re.search(r'\|([A-Z][A-Z0-9]+)_HUMAN', str(locus))
            if m: gene = m.group(1)
            if not gene:
                gn = re.search(r'GN=(\S+)', str(locus))
                if gn: gene = gn.group(1)

            # Extract spectral count
            sc = (obj.get('SpectralCount') or obj.get('spectralCount') or
                  obj.get('NSpectralCount') or 0)
            if isinstance(sc, (int, float)) and gene:
                cur = proteins.get(gene, 0)
                proteins[gene] = max(cur, int(sc))

        if proteins:
            print(f" MyProtein objects with gene+count: {len(proteins):,}")
            return proteins
        else:
            print(" NRBF parse found no MyProtein objects with gene+count.")
    except Exception as e:
        print(f"  NRBF parse error: {e}")

    #  Fallback: count PSM-level gene occurrences 
    # In .sepr2, each PSM entry stores the protein's accession.
    # We count all UniProt accession occurrences (not just per-protein-list).
    print("  Falling back to PSM-level occurrence counting...")
    text = data.decode('latin-1')

    # Find the protein sequence blocks: each protein appears ONCE in the list,
    # but each PSM references the protein accession string. If we count all
    # sp|ACCESSION| occurrences including in PSM records, high-confidence
    # proteins will appear more often.
    gene_re = re.compile(r'GN=([A-Z][A-Z0-9]{1,15})(?:\s|PE=|\x00)')
    locus_re  = re.compile(r'[st][pr]\|[A-Z0-9]+\|([A-Z0-9]+)_HUMAN')

    from collections import Counter
    gene_occ  = Counter(gene_re.findall(text))
    locus_occ = Counter(m.group(1) for m in locus_re.finditer(text))

    # Merge: prefer GN= counts as spectral count proxy
    for gene, cnt in locus_occ.items():
        if gene not in gene_occ:
            gene_occ[gene] = cnt

    print(f"  Unique genes by occurrence: {len(gene_occ):,}")
    # Return raw occurrence counts as spectral count proxy
    return dict(gene_occ)

for fname in ["ISS_WT83.sepr2", "Ground_WT83.sepr2"]:
    if not os.path.exists(fname):
        print(f"\nMissing: {fname}")
        print("Place ISS_WT83.sepr2 and Ground_WT83.sepr2 in this directory.")
        sys.exit(0)

iss_counts = extract_proteins_from_sepr2("ISS_WT83.sepr2")
ground_counts = extract_proteins_from_sepr2("Ground_WT83.sepr2")

print("\nVEN panel raw counts (ISS / Ground):")
ALL_GENES = [g for gs in VEN_PANEL.values() for g in gs] + ["NOS1"]
for gene in ALL_GENES:
    ic = iss_counts.get(gene, 0)
    gc = ground_counts.get(gene, 0)
    print(f"  {gene:10s}: ISS={ic:5d}  Ground={gc:5d}")

PSEUDO = 1
all_genes = set(iss_counts) | set(ground_counts)
log2fc = {}
for gene in all_genes:
    ic = iss_counts.get(gene, 0) + PSEUDO
    gc = ground_counts.get(gene, 0) + PSEUDO
    log2fc[gene] = np.log2(ic / gc)

log2fc_s = pd.Series(log2fc, name="log2FC_spectral")
log2fc_s.to_csv("results/alysson_all_proteins_log2FC.csv", header=True)
print(f"\nTotal proteins with log2FC: {len(log2fc_s):,}")
print(f"Genome-wide mean={log2fc_s.mean():+.4f}  SD={log2fc_s.std():.4f}")

background = log2fc_s.values
rows = []
print(f"\n{'Category':20s} {'n':>4s} {'mean log2FC':>12s} {'perm p':>10s} {'SD above null':>14s}")
print("-"*64)

for category, gene_list in VEN_PANEL.items():
    present = [g for g in gene_list if g in log2fc_s.index and log2fc_s[g] != 0.0]
    if not present:
        print(f"{category:20s} 0 --- --- ---")
        rows.append({"Category":category,"n_present":0,"mean_log2FC":np.nan,
                     "SD_above_null":np.nan,"perm_p":np.nan,"sig":""})
        continue
    cat_fc   = log2fc_s[present].values
    cat_mean = cat_fc.mean()
    null_means = np.array([
        rng.choice(background, size=len(present), replace=False).mean()
        for _ in range(N_PERM)
    ])
    null_sd = null_means.std()
    sd_above = (cat_mean - null_means.mean()) / null_sd if null_sd > 0 else np.nan
    perm_p = np.mean(np.abs(null_means - null_means.mean()) >= abs(cat_mean - null_means.mean()))
    sig = "***" if perm_p<0.001 else ("**" if perm_p<0.01 else ("*" if perm_p<0.05 else ("+" if perm_p<0.10 else "")))
    rows.append({"Category":category,"n_present":len(present),"mean_log2FC":round(cat_mean,4),
                 "SD_above_null":round(sd_above,3),"perm_p":round(perm_p,4),"sig":sig,
                 "genes_present":", ".join(present)})
    print(f"{category:20s} {len(present):>4d} {cat_mean:>+12.4f} {perm_p:>10.4f}{sig} {sd_above:>14.2f}")

results_df = pd.DataFrame(rows)
results_df.to_csv("results/alysson_VEN_panel_WT83_ISS_vs_Ground.csv", index=False)
print("\nSaved: results/alysson_VEN_panel_WT83_ISS_vs_Ground.csv")

CAT_ORDER  = ["Myelination","FastSignalling","SocialCircuit","LayerVProjection","MetabolicSupport"]
CAT_COLORS = {"Myelination":"#4472C4","FastSignalling":"#70AD47",
              "SocialCircuit":"#ED7D31","LayerVProjection":"#FFC000",
              "MetabolicSupport":"#7030A0"}
CAT_LABELS = ["Myelin-\nation","Fast\nSignal.","Social\nCircuit","Layer V\nProj.","Metabolic\nSupport"]

res_idx = results_df.set_index("Category")
sds = [float(res_idx.loc[c,"SD_above_null"]) if c in res_idx.index and not np.isnan(float(res_idx.loc[c,"SD_above_null"])) else np.nan for c in CAT_ORDER]
pps = [float(res_idx.loc[c,"perm_p"]) if c in res_idx.index else np.nan for c in CAT_ORDER]

fig, ax = plt.subplots(figsize=(9,5))
valid = [v for v in sds if not np.isnan(v)]
if valid:
    ax.bar(range(len(CAT_ORDER)), [v if not np.isnan(v) else 0 for v in sds],
           color=[CAT_COLORS[c] for c in CAT_ORDER], edgecolor="black", lw=0.8, alpha=0.88)
    ax.axhline(0, color="black", lw=0.9)
    ax.axhline( 1.96, color="#888888", lw=0.9, ls="--", alpha=0.65, label="p=0.05 (±1.96 SD)")
    ax.axhline(-1.96,  color="#888888", lw=0.9, ls="--", alpha=0.65)
    ylo, yhi = min(valid)-0.5, max(valid)+1.2
    ax.set_ylim(ylo, yhi)
    for i,(sd,pp) in enumerate(zip(sds,pps)):
        if np.isnan(sd) or np.isnan(pp): continue
        if pp < 0.10:
            sym = "***" if pp<0.001 else ("**" if pp<0.01 else ("*" if pp<0.05 else "+"))
            ax.text(i, sd+(yhi-ylo)*0.025, sym, ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(CAT_ORDER)))
ax.set_xticklabels(CAT_LABELS, fontsize=9, linespacing=1.3)
ax.set_ylabel("SD above genome-wide permutation null\n(N=10,000 random proteins, seed 42)", fontsize=9)
ax.set_title("VEN Gene Panel — Alysson Lab Proteomics (PXD069807)\n"
             "WT83 ISS vs Ground (30 days) · Jourdon et al. 2026", fontsize=10, fontweight="bold")
ax.legend(loc="upper right", fontsize=8.5)
ax.yaxis.grid(True, alpha=0.22, lw=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("results/alysson_VEN_panel_permutation.png", dpi=200, bbox_inches="tight")
plt.close()
print("Saved: results/alysson_VEN_panel_permutation.png\nDone.")
