"""
Inspect .sepr2 binary around a known protein to find the spectral count field.
"""
import re
import struct

FILE = "ISS_WT83.sepr2"

with open(FILE, 'rb') as f:
    data = f.read()

text = data.decode('latin-1')

# Find NEFM (neurofilament medium) - a large structural protein, likely high spectral count
# We look for the full FASTA header and then inspect the bytes around it
target_patterns = ["GN=NEFM", "GN=ANK3", "GN=SNAP25", "GN=SYP", "GN=NRXN1"]

for pat in target_patterns:
    positions = [m.start() for m in re.finditer(re.escape(pat), text)]
    print(f"\n{pat}: found at {len(positions)} positions")
    if positions:
        # Look at first 3 occurrences
        for pos in positions[:3]:
            # Context: 50 chars before, 150 after
            ctx_start = max(0, pos - 50)
            ctx_end = min(len(data), pos + 150)
            snippet = data[ctx_start:ctx_end]
            # Print hex dump of region after the GN= tag
            after_pos = pos + len(pat)
            after_bytes = data[after_pos:after_pos+40]
            hex_str = ' '.join(f'{b:02x}' for b in after_bytes)
            # Try to read integers and floats from those bytes
            vals = []
            for offset in range(0, 36, 4):
                if offset + 4 <= len(after_bytes):
                    i32 = struct.unpack_from('<i', after_bytes, offset)[0]
                    vals.append(f"int32@{offset}={i32}")
            for offset in range(0, 32, 8):
                if offset + 8 <= len(after_bytes):
                    f64 = struct.unpack_from('<d', after_bytes, offset)[0]
                    if 0 < abs(f64) < 1e9:
                        vals.append(f"float64@{offset}={f64:.4f}")
            # Text after
            after_text = data[after_pos:after_pos+80].decode('latin-1')
            printable = ''.join(c if 32 <= ord(c) < 127 else '.' for c in after_text)
            print(f" pos={pos}: hex={hex_str[:60]}")
            print(f" text: {printable[:80]}")
            print(f" numeric candidates: {vals[:6]}")
