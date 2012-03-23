#!/usr/bin/python

def getid(line):
    parts=line.split("\"")
    return eval(parts[1])

def getiso(line):
    parts=line.split("\"")
    return eval(parts[7])

def main():
    f = open('weights/unbiased_OptimizedCuts_chaIsoOnly_weights_eta14442_Pt10To20.xml', 'r')

    wplines = [251,281,293,311,323]
    wp = [70,80,85,90,95]
    
    files = ['weights/unbiased_OptimizedCuts_weights_eta10_Pt10To20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta14442_Pt10To20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta1566_Pt10To20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta20_Pt10To20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta25_Pt10To20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta10_Pt20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta14442_Pt20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta1566_Pt20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta20_Pt20.xml',
             'weights/unbiased_OptimizedCuts_weights_eta25_Pt20.xml']


    etabins = [0,1,2,3,4,0,1,2,3,4]
    ptbins =  [0,0,0,0,0,1,1,1,1,1]

    etalabels = ['/&eta;/<1','1</&eta;/<1.4442','1.4442</&eta;/<1.566','1.566</&eta;/<2.0','2.0</&eta;/<2.5']

    # 2 pt bins x 5 eta bins x 5 WPs
    idtable = [ [ [ 999 for i in range(5) ] for j in range(5) ] for k in range(2) ]
    isotable = [ [ [ 999 for i in range(5) ] for j in range(5) ] for k in range(2) ]

    bin=0
    for ifile in files:
        etabin = etabins[bin]
        ptbin = ptbins[bin]

        eleid = []
        eleiso = []
        f = open(ifile)
        counter=0
        for line in f:
            counter = counter+1
            for i in range(len(wp)):
                if(counter==wplines[i]):
                    eleid.append(getid(line))
                    eleiso.append(getiso(line))

        for i in range(len(eleid)):
            idtable[ptbin][etabin][i] = eleid[i]
            isotable[ptbin][etabin][i] = eleiso[i]
        bin = bin+1


    # write the table now
    for ipt in range (2):
        print "pt bin = "+str(ipt)
        print "| *WP* | ",
        for ieta in range(5):
            print etalabels[ieta]+" | ",
            if(ieta==4):
                print ""
        for iwp in range(5):
            print "| *WP "+str(wp[iwp])+"* | ",
            for ieta in range(5):
                print "%.3f" % idtable[ipt][ieta][iwp]+" / %.3f" % isotable[ipt][ieta][iwp]+" | ",
                if(ieta==4):
                    print ""


if __name__ == "__main__":
    main()
