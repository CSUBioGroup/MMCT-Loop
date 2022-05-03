import sys
with open(sys.argv[1],'r') as fin:
    line = fin.readline()
    line = fin.readline()
    while line != "":
        line = line.strip().split()
        if int(float(line[-1])) == 1:
            sa = line[-3]
            sb = line[-2]
            namea,alr = sa.split(":")
            nameb,blr = sb.split(":")
            al, ar = alr.split("-")
            bl, br = blr.split("-")
            print("\t".join([namea,al,ar,nameb,bl,br]))
        line = fin.readline()
