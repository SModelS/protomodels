f = open('Obs_UL.csv', 'r')

a,b = [],[]
for line in f:
    if line[0] != '#' and line[0] != "" and len(line) > 1:
        if line[1] != '$':
            l = line[1:-2].split()[0]
            for i in range(len(l)):
                if l[i] == ",":
                    a.append(float(l[:i-1]))
                    b.append(float(l[i+2:]))

f.close()

x,y,z = [],[],[]
lh = int(len(a)*0.5)                
for i in range(lh):
    x.append(a[i])
    y.append(b[i])
    z.append(b[i+lh])

f = open('Obs_UL_fixed.csv', 'w')

for i in range(len(x)):
    f.write(str(x[i]) + ',' + str(y[i]) + ',' + str(z[i]) + '\n')
    
f.close()
