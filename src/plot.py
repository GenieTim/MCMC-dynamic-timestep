import glob
import json
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# manually decided:
# NOTE: either 4 & 5 + 6 & 7 are CO & NH and 14 & 15 + 16 & 17
# OR: 4 & 5 + 6 & 7 are NH & CO and 14 & 15 + 16 & 17
atom_types = {
    'H': [0, 2, 3, 7, 9, 11, 12, 13, 17, 19, 20, 21],
    'C': [1, 4, 8, 10, 14, 18],
    'N': [6, 16],
    'O': [5, 15]
}

# NOTE: switch if above they were wrong
dihedral_angle_atoms = {
    'phi': [
        6, 7, 8, 10
    ],
    'psi': [
        14, 15, 8, 10
    ]
}


def get_atom_type(idx):
    for atom_type, idxs in atom_types.items():
        if (idx in idxs):
            return atom_type


def points_to_vector(a1, a2):
  # utility to get a directional vector from two points
    return [
        a1[0] - a2[0],
        a1[1] - a2[1],
        a1[2] - a2[2]
    ]


def angle_between_vectors(v1, v2):
    # utility to calculate angle between two vectors
    # v1 is your first vector
    # v2 is your second vector
    return np.arccos(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))


def torsional_angle_from_atoms(a1, a2, b1, b2):
    # utility to calculate angle between two lines given by points
    return angle_between_vectors(points_to_vector(a1, a2), points_to_vector(b1, b2))


def torsional_angle(name, positions):
  # utility to calculate the torsional angle from the whole positions list, by name of the torsional angle
  # see: dihedral_angle_atoms
    idxs = dihedral_angle_atoms[name]
    a1 = positions[idxs[0]]
    a2 = positions[idxs[1]]
    b1 = positions[idxs[2]]
    b2 = positions[idxs[3]]
    return torsional_angle_from_atoms(a1, a2, b1, b2)


# read data
os.chdir(os.path.dirname(os.path.realpath(__file__)) + "/../")
# TODO: replace with your experiment
# "./out/results_TorsionNeglectingMove.json"))
data = json.load(open("./out/test_results_False_1.json"))

# check whether we know which positions correspond to which atom
# test_pos = np.array(data[0]['positions']['value'])

# fig = plt.figure()
# ax3D = fig.add_subplot(111, projection='3d')
# ax3D.scatter(test_pos[:, 0], test_pos[:, 1], test_pos[:, 2], marker='o')
# indices = range(0, len(test_pos))
# for i in range(len(test_pos)):
#     ax3D.text(test_pos[i, 0], test_pos[i, 1],
#               test_pos[i, 2], '{}: {}'.format(str(i), str(get_atom_type(i))))

# plt.show()

#
energies = []
phi = []
psi = []

for i in range(len(data)):
    state = data[i]
    energies.append(state['unmodified_energy'])
    positions = np.array(state['positions']['value'])
    phi.append(torsional_angle('phi', positions) / 1e-1)
    psi.append(torsional_angle('psi', positions) / 1e-1)

print("Energy range: between {} and {}".format(max(energies), min(energies)))
norm = mpl.colors.Normalize(vmin=min(energies), vmax=max(energies))
color_map = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.hot)
colors = []
for i in range(len(energies)):
    colors.append(color_map.to_rgba(energies[i]))

fig, ax = plt.subplots()
ax.scatter(phi, psi, c=colors, marker="+")
# emulate fix for https://github.com/matplotlib/matplotlib/issues/6015
# see also https://stackoverflow.com/questions/45757260/matplotlibs-autoscale-doesnt-seem-to-work-on-y-axis-for-small-values
ax.set_ylabel("Psi $\Psi$")
ax.set_xlabel("Phi $\Phi$")

dx = (max(phi) - min(phi))*0.1
dy = (max(psi) - min(psi))*0.1

ax.set_xlim(min(phi)-dx, max(phi)+dx)
ax.set_ylim(min(psi)-dy, max(psi)+dy)

ax.yaxis.set_ticks(np.linspace(min(psi), max(psi), num=5))
ax.xaxis.set_ticks(np.linspace(min(phi), max(phi), num=5))

fig.savefig("./img/psi_vs_phi.png", bbox_inches='tight')
# plt.show()

# plot runtime gain
perf_test_res = pd.read_csv('./out/test_results_perf_True.tsv', sep="\t")

fig, ax = plt.subplots()
ax.plot(perf_test_res['Timestepsize'], perf_test_res['Time [s]'])
ax.set_title("Runtime of the MC simulation, accepting all steps")
ax.set_xlabel("Nr. of timesteps before updating torsion cache")
ax.set_ylabel("Runtime [s]")
ax.set_ylim(0, 1.05 * max(perf_test_res['Time [s]']))
fig.savefig("./img/runtime.png", bbox_inches='tight')

# plot deviations between current energy
acceptrate_test_res = pd.read_csv(
    './out/test_results_perf_False.tsv', sep="\t")

fig, ax = plt.subplots()
ax.plot(acceptrate_test_res['Timestepsize'],
        acceptrate_test_res['Accepted'] / acceptrate_test_res['Proposed'])
ax.set_title("Acceptance rate of the MC simulation")
ax.set_xlabel("Nr. of timesteps before updating torsion cache")
ax.set_ylabel("Acceptance rate")
ax.set_ylim(0, 1)
fig.savefig("./img/acceptrate.png", bbox_inches='tight')

# TODO: plot deviations between current energy
mc_res_files = sorted(glob.glob("./out/test_results_False_*.json"))
timesteps = []
max_energy_deviations = []

for f in mc_res_files:
    timestep_re = re.search('([0-9]*)\.json', f)
    timestep = timestep_re.group(1)
    timesteps.append(timestep)
    data = pd.read_json(f)
    
    energy_deviations_this_time = data['proposed_energy'] - data['unmodified_energy']
    max_energy_deviations.append(max(max(energy_deviations_this_time), abs(min(energy_deviations_this_time))))
    # TODO: may want to plot these per timestep here
    # instead of just the average afterwards

fig, ax = plt.subplots()
ax.plot(timesteps,
        max_energy_deviations)
ax.set_title("Energy deviation from correct value while iterating")
ax.set_xlabel("Nr. of timesteps before updating torsion cache")
ax.set_ylabel("Maximum absolute energy deviation [kJ/mol]")
fig.savefig("./img/energy_deviations.png", bbox_inches='tight')

