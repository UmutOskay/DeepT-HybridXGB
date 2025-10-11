#ITTCA-RF- Feature extraction
#!/usr/bin/env python
#_*_coding:utf-8_*_

from collections import Counter
import numpy as np
import re
import math


def Count_1(seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum
def Count_2(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 0 for i in cutoffNums]
    code = []
    for cutoff in cutoffNums:
        myCount = 0
        if cutoff == 0:
            code.append(0)
        else:
            for i in range(len(sequence)):
                if sequence[i] in aaSet:
                    myCount += 1
                    if myCount == cutoff:
                        code.append((i + 1) / len(sequence) * 100)
                        break
            if myCount == 0:
                code.append(0)
    return code

AA = 'ACDEFGHIKLMNPQRSTVWY'
group = {
        'Alphaticr': 'GAVLMI',
        'Aromatic': 'FYW',
        'Positive Charger': 'KRH',
        'Negative Charger': 'DE',
        'Uncharger': 'STCPNQ'}

group1 = {
        'Hydrophobicity': 'RKEDQN',  # 疏水性特征
        'Normwaals Volume': 'GASCTPD',  # 范德华力
        'Polarity': 'LIFWCMVY',  # 极性
        'Polarizability': 'GASDT',  # 极化性质
        'Charge': 'KR',  # 电荷性质
        'Surface Tension': 'GQDNAHR',  # 表面张力
        'Secondary Struct': 'EALMQKRH',  # 二级结构
        'Solvent Access': 'ALFCGIVW'}  # 溶剂可及性

group2 = {
        'Hydrophobicity': 'GASTPHY',
        'Normwaals Volume': 'NVEQIL',
        'Polarity': 'PATGS',
        'Polarizability': 'CPNVEQIL',
        'Charge': 'ANCQGHILMFPSTWYV',
        'Surface Tension': 'KTSEC',
        'Secondary Struct': 'VIYCWFT',
        'Solvent Access': 'RKQEND'}
group3 = {
        'Hydrophobicity': 'CVLIMFW',
        'Normwaals Volume': 'MHKFRYW',
        'Polarity': 'HQRKNED',
        'Polarizability': 'KMHFRYW',
        'Charge': 'DE',
        'Surface Tension': 'ILMFPWYV',
        'Secondary Struct': 'GNPSD',
        'Solvent Access': 'MPSTHY'}
property = ('Hydrophobicity', 'Normwaals Volume',
                'Polarity', 'Polarizability', 'Charge', 'Surface Tension', 'Secondary Struct', 'Solvent Access')

def get_feature_names():
    AA = 'ACDEFGHIKLMNPQRSTVWY'  # AA değişkenini burada tanımlıyoruz.
    feature_names = []

    # AAC features
    feature_names.extend([f"AAC_{aa}" for aa in AA])

    # CTDC, CTDT, CTDD features
    for p in property:
        feature_names.extend([f"CTDC {p} c1", f"CTDC {p} c2", f"CTDC {p} c3"])
        feature_names.extend([f"CTDT {p} c1221", f"CTDT {p} c1331", f"CTDT {p} c2332"])
        feature_names.extend([f"CTDD {p} Cutoff {i}" for i in range(1, 16)])

    # GAAC features
    feature_names.extend([f"GAAC {g}" for g in group.keys()])

    # GDPC features
    groupKey = group.keys()
    dipeptide = [f"{g1}.{g2}" for g1 in groupKey for g2 in groupKey]
    feature_names.extend([f"GDPC {d}" for d in dipeptide])

    # GTPC features
    triple = [f"{g1}.{g2}.{g3}" for g1 in groupKey for g2 in groupKey for g3 in groupKey]
    feature_names.extend([f"GTPC {t}" for t in triple])

    # PAAC features
    dataFile = '/content/PAAC.txt'
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    feature_names.extend([f"PAAC {aa}" for aa in AA])
    for n in range(1, min(3, len(AA))):
        feature_names.extend([f"PAAC Theta {n}"])

    return feature_names



def get_features(fastas):
    def AAC():
        encoding = []
        for sequence in fastas:
            count = Counter(sequence)
            for key in count:
                count[key] = count[key]/len(sequence) * 100
            code = []
            for aa in AA:
                code.append(count[aa])
            encoding.append(code)
        return encoding

    def CTDC(p):
        encodings = []
        for sequence in fastas:
            code = []
            c1 = Count_1(group1[p], sequence)/len(sequence)*100
            c2 = Count_1(group2[p], sequence)/len(sequence)*100
            c3 = 100 - c1 - c2
            code = code + [c1, c2, c3]
            encodings.append(code)
        return encodings


    def CTDT(p):
        encodings = []
        for sequence in fastas:
            code = []
            aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
            if not aaPair:
              code = [0,0,0]
            else:
              c1221, c1331, c2332 = 0, 0, 0
              for pair in aaPair:
                  if (pair[0] in group1[p] and pair[1] in group2[p]) or (
                          pair[0] in group2[p] and pair[1] in group1[p]):
                      c1221 = c1221 + 1
                      continue
                  if (pair[0] in group1[p] and pair[1] in group3[p]) or (
                          pair[0] in group3[p] and pair[1] in group1[p]):
                      c1331 = c1331 + 1
                      continue
                  if (pair[0] in group2[p] and pair[1] in group3[p]) or (
                          pair[0] in group3[p] and pair[1] in group2[p]):
                      c2332 = c2332 + 1
              code = code + [c1221/len(aaPair)*100, c1331/len(aaPair)*100, c2332/len(aaPair)*100]
            encodings.append(code)
        return encodings

    def CTDD(p):
        encodings = []
        for sequence in fastas:
            code = []
            code = code + Count_2(group1[p], sequence) + Count_2(group2[p], sequence) + Count_2(group3[p], sequence)
            encodings.append(code)
        return encodings

    def GAAC():
        encoding = []
        groupKey = group.keys()
        for sequence in fastas:
            code = []
            count = Counter(sequence)
            myDict = {}
            for key in groupKey:
                for aa in group[key]:
                    myDict[key] = myDict.get(key, 0) + count[aa]
            for key in groupKey:
                code.append(myDict[key] / len(sequence))
            encoding.append(code)
        return encoding

    def GDPC():
        groupKey = group.keys()
        #baseNum = len(groupKey)
        dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]
        index = {}
        for key in groupKey:
            for aa in group[key]:
                index[aa] = key
        encodings = []
        for sequence in fastas:
            code = []
            myDict = {}
            for t in dipeptide:
                myDict[t] = 0

            sum = 0
            for j in range(len(sequence) - 2 + 1):
                myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                    sequence[j + 1]]] + 1
                sum = sum + 1

            if sum == 0:
                for t in dipeptide:
                    code.append(0)
            else:
                for t in dipeptide:
                    code.append(myDict[t] / sum)
            encodings.append(code)
        return encodings

    def GTPC():
        groupKey = group.keys()
        baseNum = len(groupKey)
        triple = [g1 + '.' + g2 + '.' + g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]
        index = {}
        for key in groupKey:
            for aa in group[key]:
                index[aa] = key
        encodings = []

        for sequence in fastas:
            code = []
            myDict = {}
            for t in triple:
                myDict[t] = 0

            sum = 0
            for j in range(len(sequence) - 3 + 1):
                myDict[index[sequence[j]] + '.' + index[sequence[j + 1]] + '.' + index[sequence[j + 2]]] = myDict[index[sequence[j]] + '.' +index[sequence[j + 1]] + '.' +index[sequence[j + 2]]] + 1
                sum = sum + 1
            if sum == 0:
                for t in triple:
                    code.append(0)
            else:
                for t in triple:
                    code.append(myDict[t] / sum)
            encodings.append(code)
        return encodings

    def Rvalue(aa1, aa2, AADict, Matrix):
        return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)


    def PAAC():
        dataFile = '/content/PAAC.txt' #r'PAAC.TXT'
        with open(dataFile) as f:
            records = f.readlines()
        AA = ''.join(records[0].rstrip().split()[1:])
        # AA = "ARNDCQEGHILKMFPSTWYV"
        AADict = {}
        for i in range(len(AA)):  # 20
            AADict[AA[i]] = i
        AAProperty = []
        AAPropertyNames = []
        for i in range(1, len(records)):  # llen(records) 4
            array = records[i].rstrip().split() if records[i].rstrip() != '' else None
            AAProperty.append([float(j) for j in array[1:]])
            AAPropertyNames.append(array[0])

        AAProperty1 = []
        for i in AAProperty:
            meanI = sum(i) / 20
            fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
            AAProperty1.append([(j - meanI) / fenmu for j in i])

        encodings = []

        for sequence in fastas:
            code = []
            theta = []
            for n in range(1, min(3, len(sequence))):
            #for n in range(1, 3):
                theta.append(
                    sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in
                         range(len(sequence) - n)]) / (
                            len(sequence) - n))
            myDict = {}
            for aa in AA:
                myDict[aa] = sequence.count(aa)
            code = code + [myDict[aa] / (1 + 0.05 * sum(theta)) for aa in AA]
            code = code + [(0.05 * j) / (1 + 0.05 * sum(theta)) for j in theta]
            encodings.append(code)
        return encodings

    print('Feature extraction...')
    encoding = []
    encoding.append(AAC())
    for p in property:
        encoding.append(CTDC(p))
        encoding.append(CTDT(p))
        encoding.append(CTDD(p))
    encoding.append(GAAC())
    encoding.append(GDPC())
    encoding.append(GTPC())
    encoding.append(PAAC())
    return np.column_stack(encoding)


