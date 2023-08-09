import sys
with open(sys.argv[2],'r') as fin:
    line = fin.readline()
    line = fin.readline()
    if sys.argv[1] == '-b':
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
    else:
        while line != "":
            line = line.strip().split()
            if int(float(line[-1])) == 1:
                sa = line[-3]
                sb = line[-2]
                namea,alr = sa.split(":")
                nameb,blr = sb.split(":")
                al, ar = alr.split("-")
                bl, br = blr.split("-")
                es, fdr, pvalue = line[5], line[6], line[7]
                print("\t".join([namea,al,ar,nameb,bl,br,es,fdr,pvalue]))
            line = fin.readline()