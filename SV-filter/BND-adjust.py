import pandas as pd

df = pd.read_csv("INV.breakpoint.txt", sep="\t")

bnd_df = df[df['ID'].str.contains('BND')].copy()

tolerance = 100

pairs = []
used = set()

for i, row1 in bnd_df.iterrows():
    for j, row2 in bnd_df.iterrows():
        if i == j or (j, i) in used:
            continue
        if (row1['CHR'] == row2['CHR2'] and row1['CHR2'] == row2['CHR'] and
            abs(row1['POS'] - row2['END']) <= tolerance and
            abs(row1['END'] - row2['POS']) <= tolerance):
            pairs.append({
                'id1': row1['ID'],
                'id2': row2['ID'],
                'idx1': i,
                'idx2': j,
                'start': min(row1['POS'], row1['END']),
                'end': max(row1['POS'], row1['END']),
                'svlen': max(int(row1['SVLEN']), int(row2['SVLEN']))
            })
            used.add((i, j))

to_remove = set()
import itertools
for a, b in itertools.combinations(pairs, 2):
    if max(a['start'], b['start']) <= min(a['end'], b['end']):  # overlap
        if a['svlen'] >= b['svlen']:
            to_remove.update([b['id1'], b['id2']])
        else:
            to_remove.update([a['id1'], a['id2']])

paired_ids = set([p['id1'] for p in pairs] + [p['id2'] for p in pairs])
all_bnd_ids = set(bnd_df['ID'])
unpaired_ids = all_bnd_ids - paired_ids

final_remove = to_remove.union(unpaired_ids)

with open("deleted_ids.txt", "w") as f:
    for id in sorted(final_remove):
        f.write(f"{id}\n")

filtered_df = df[~df['ID'].isin(final_remove)]
filtered_df.to_csv("INV.filtered.txt", sep="\t", index=False)
print("筛选完成，保留写入 INV.filtered.txt，剔除ID写入 deleted_ids.txt")
