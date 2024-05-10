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

data = read_field('2DWField.csv')
nt = root.TNtupleD("weighting","Weighting Field","x:y:wx:wy")
n = len(data)
print 'read: ',n

ff = root.TFile('2DMgem1Volt.root','recreate')
for (x,y,ex,ey) in data:
    nt.Fill(x,y,ex,ey)

nt.Write()
ff.Close()
