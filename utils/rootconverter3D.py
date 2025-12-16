import csv
import itertools
import ROOT as root

def read_field(filename):
    with open(filename,"r") as f:
        lines = itertools.islice(f,9,None)
        reader = csv.reader(lines)
        data = []
        for line in reader:
            point = []
            for col in line:
                point.append(float(col))
            data.append(tuple(point))
    return data

data = read_field('2DGGemE-field.csv')
nt = root.TNtupleD("weighting","Weighting Field","x:y:z:wx:wy:wz")
n = len(data)
print('read: ',n)

ff = root.TFile('2DGGEM_5_2_5-1V.root','recreate')
for (x,y,z,ex,ey,ez) in data:
    nt.Fill(x,y,z,ex,ey,ez)

nt.Write()
ff.Close()

#nt = root.TNtupleD("drift","Drift Field","x:y:z:ex:ey:ez")
#
