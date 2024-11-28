import numpy as np
import math

def Dictionalyofq(reads):
    '''Encoder for quality score'''
    qscoredict = {"!":0, '"':1,"#":2,'$':3,'%':4,
                  '&':5,"'":6,"(":7,")":8,"*":9,"+":10,
                  ",":11,"-":12,".":13,"/":14,"0":15,
                  "1":16,"2":17,"3":18,"4":19,"5":20,
                  "6":21,"7":22,"8":23,"9":24,":":25,
                  ";":26,"<":27,"=":28,">":29,"?":30,
                  "@":31,"A":32,"B":33,"C":34,"D":35,
                  "E":36,"F":37,"G":38,"H":39,"I":40,
                  "J":41,"K":42,"L":43,"M":44,"N":45,
                  "O":46,"P":47,"Q":48,"R":49,"S":50,
                  "T":51,"U":52,"V":53,"W":54,"X":55,
                  "Y":56,"Z":57,"[":58,"\\":59,"]":60,
                  "^":61,"_":62,"`":63,"a":64,"b":65,
                  "c":66,"d":67,"e":68,"f":69,"g":70,
                  "h":71,"i":72,"j":73,"k":74,"l":75,
                  "m":76,"n":77,"o":78,"p":79,"q":80,
                  "r":81,"s":82,"t":83,"u":84,"v":85,
                  "w":86,"x":87,"y":88,"z":89,"{":90,
                  "|":91,"}":92,"~":93}
    return qscoredict

def valueofquality(reads):
    '''Total quality score per read'''
    total = 0
    qscoredict = Dictionalyofq(reads)
    for i in reads:
        if i in qscoredict:
            value = qscoredict[i]
            total+= value
    return total
               
def calmean(reads):
    '''Calculate mean quality score per read'''
    lengthofreads=len(reads)
    qualityscores=valueofquality(reads)/lenofreads(reads)
    # print(qualityscores)
    return qualityscores

def lenofreads(reads):
    '''Count length'''
    return len(reads)

def SDquality(reads):
    '''Calculate standard deviation of quality score per read'''
    sumxmutiply = 0
    summean = (valueofquality(reads))**2   
    qscoredict = Dictionalyofq(reads)
    for i in reads:
        if i in qscoredict:
            valusesq = qscoredict[i]
            sumxmutiply += valusesq**2
    if lenofreads(reads) == 1:
        return 0
    return math.sqrt(((lenofreads(reads)*sumxmutiply) - summean)/(lenofreads(reads)*(lenofreads(reads)-1)))

def calquantile(reads, IQR = False):
    '''Calculate median quality score per read'''
    qscoredict = Dictionalyofq(reads)
    value = []
    for i in reads:
        if i in qscoredict:
            value.append(qscoredict[i])
    
    median = np.quantile(value, .50)
    if IQR == True:
        Q1 = np.quantile(value, .25)
        Q3 = np.quantile(value, .75)
        IQR = Q3 - Q1
        return median, IQR

    return median
      



def callenandbar(dictcal,key):
    '''Calculation for each barcode'''
    #  print(dictcal)
    total = 0
    sumxmutiply = 0
    for i in dictcal[key]:
        # print(i)
        total += i
        sumxmutiply += i**2

    summean = (total)**2
    sd = math.sqrt(((len(dictcal[key])*sumxmutiply) - summean)/(len(dictcal[key])*(len(dictcal[key])-1)))
    minvalue = min(dictcal[key])
    maxvalue = max(dictcal[key])
    mean_q = total/len(dictcal[key])
    return mean_q, sd, minvalue, maxvalue

def callenandbar_median(dictcal,key):
    '''Calculation for each barcode by median'''

    Q1 = np.quantile(dictcal[key], .25)
    median = np.quantile(dictcal[key], .50)
    Q3 = np.quantile(dictcal[key], .75)
    IQR = Q3 - Q1
    minvalue = min(dictcal[key])
    maxvalue = max(dictcal[key])

    return median, IQR, minvalue, maxvalue

if __name__ == '__main__':
   reads = "ASWDJIDJ(&"
#    quality_values = valueofquality(reads)
#    print("Total of q score : ", quality_values)

#    lengthofreads=len(reads)
#    print(lengthofreads)

#    qualityscores= quality_values/lengthofreads
#    print(qualityscores)

   print(calquantile(reads))