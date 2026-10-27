"""Replace the simulation tables of the manuscript by freshly generated ones (e.g. after the full run of 500 replications).

usage:  python3 refresh_manuscript_tables.py <manuscript.tex> <tex_dir>

Every table environment of the manuscript whose \\label{...} matches the label of a table file in <tex_dir> (the output
of summarize_sim.R) is replaced by the content of that file; commented-out tables and tables without a generated
counterpart are left as they are. The figures (simcase<c>.pdf, simS_bylevel.pdf) only have to be copied next to the
manuscript. The numbers quoted in the text are not changed and have to be checked against key_numbers.txt.
"""
import glob, os, re, sys

tex_file, tex_dir = sys.argv[1], sys.argv[2]
main = open(tex_file).read()

generated = {}
for fn in glob.glob(os.path.join(tex_dir, "*.tex")):
    body = open(fn).read().rstrip() + "\n"
    for lab in re.findall(r"\\label\{([^}]*)\}", body):
        generated[lab] = (fn, body)

table = re.compile(r"^\\begin\{table\}.*?^\\end\{table\}\n", re.S | re.M)   # uncommented tables only (start of line)
done = []

def swap(m):
    labs = re.findall(r"\\label\{([^}]*)\}", m.group(0))
    hit = [l for l in labs if l in generated]
    if not hit:
        return m.group(0)
    done.append(hit[0])
    return generated[hit[0]][1]

main = table.sub(swap, main)
open(tex_file, "w").write(main)
print("replaced %d tables: %s" % (len(done), ", ".join(done)))
missing = sorted(set(generated) - set(done))
if missing:
    print("generated but not found in the manuscript:", ", ".join(missing))
