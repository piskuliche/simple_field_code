import MDAnalysis as mda
import sys

grofile = str(sys.argv[1])
xtcfile = str(sys.argv[2])

u = mda.Universe(grofile, xtcfile)
aa = u.select_atoms("all")

with mda.Writer("traj.xyz", aa.n_atoms) as W:
    f = open("L.dat", 'w')
    for ts in u.trajectory:
        W.write(aa)
        f.write("%10.5f %10.5f %10.5f\n" %
                (ts.dimensions[0], ts.dimensions[1], ts.dimensions[2]))
    f.close()
