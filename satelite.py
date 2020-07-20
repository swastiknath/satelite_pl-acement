import itertools
import random
import dimod
import neal
import networkx as nx


num_satelites  = 12
num_constellations = 4

constellation_size = num_satelites // num_constellations


coverage_map ={0: 0.90,
            1: 0.36,
            2: 0.79,
            3: 0.78,
            4: 0.46,
            5: 0.27,
            6: 0.86,
            7: 0.52,
            8: 0.78,
            9: 0.99,
            10: 0.25,
            11: 0.91}

score_threshold = 0.4

bqm = dimod.BinaryQuadraticModel.empty(dimod.BINARY)

for constellation in itertools.combinations(range(num_satelites), constellation_size):
    score = sum(coverage_map[v] for v in constellation) / constellation_size
    if score < score_threshold :
        continue
    bqm.add_variable(frozenset(constellation), -score)


for c0, c1 in itertools.combinations(bqm.variables, 2):
    if c0.isdisjoint(c1):
        continue
    bqm.add_interaction(c0, c1, 2)

bqm.update(dimod.generators.combinations(bqm, num_constellations, strength=1))

sampleset = neal.Neal().sample(bqm, num_reads=100).aggregate()
constellations = [constellation for constellation, chosen in sampleset.first.sample.items() if chosen]
print(constellations)