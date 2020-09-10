# coding: utf-8
import sys
import json
data = json.load(open('smirnoff_parameter_assignments.json'))
mset = []
a = sys.argv[1]
for m, terms in data.items():
    for t in terms:
        if t['id'] == a:
            mset.append([m, t])
            
print(f'Angles matching {a}')
for m, t in mset:
    print(m, t)

print('Torsions containing above angles')            
for m, t_angle in mset:
    for t in data[m]:
        if t['type'] == 'ProperTorsions' and set(t['atoms']).issuperset(set(t_angle['atoms'])):
            print(m, t)
            
            
