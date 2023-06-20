"""****************************************************************************
**************** RLP prediction by machine learning aproaches *****************
* Autores: Silva et al. 2019
* Date:14/05/2019
* Version: 0.9 alpha
****************************************************************************"""
""" **************************************************************************************************************************** """
""" ************************************************************************************************************************** """
##################################################### Declarations #############################################################

import re
import sys
import math
import pickle
import pymc3 as pm
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from collections import defaultdict
from FeacturesExtraction import FeacturesExtraction
AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

AADipeptide = {}
for i in AA:
    for j in AA:
        Dipeptide = i+j
        AADipeptide[Dipeptide] = int(1)

AATripeptide = {}
for i in AA:
    for j in AA:
        for k in AA:
            kmer = (i+j+k)
            AATripeptide[kmer] = int(1)

###################################################### AA composition ###########################################################
# AAComposition
# Dipepitide
# Tripepitide


def CalculateAAComposition(aaseq):
    results = {}
    for i in AA:
        results[i] = 0 if aaseq.count(i) == 0 else round(
            float(aaseq.count(i)) / float(len(aaseq)), 5)
    return results


def CalculateDipeptide(aaseq):
    results = {}
    for i in AA:
        for j in AA:
            Dipeptide = i+j
            results[Dipeptide] = 0 if aaseq.count(Dipeptide) == 0 else round(
                float(aaseq.count(Dipeptide)) / float(len(aaseq)), 5)
    return results


def CalculateTripeptide(aaseq):
    results = {}
    for aa in AATripeptide:
        results[aa] = 0 if aaseq.count(aa) == 0 else round(
            float(aaseq.count(aa)) / float(len(aaseq)), 5)
    return results


""" *********************************************************************** """


def Run_Signalp(filefasta):
    dic = defaultdict(list)
    cmd = 'signalp/signalp ' + filefasta + ' > ' + filefasta + '_signalp.out'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()

    out = '' + filefasta + '_signalp.out'

    f = open(out, 'r')
    for l in f:
        if '#' not in l:
            l = re.sub(r"\s+", ' ', l)
            l = l.split(' ')
            dic[l[0]] = l[1::]
    return dic


""" *********************************************************************** """


def Run_Tmhmm(filefasta):
    dic = defaultdict(list)
    cmd = './tmhmm/bin/tmhmm ' + filefasta + ' > ' + filefasta + '_tmhmm.out'

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()

    out = '' + filefasta + '_tmhmm.out'

    f = open(out, 'r')
    for l in f:
        if 'Total prob' in l:
            l = re.sub(r"\s+", ' ', l)
            l = l.split(' ')
            dic[l[1]].append(l[6])
        if 'TMhelix' in l:
            l = re.sub(r"\s+", ' ', l)
            l = l.split(' ')
            if int(l[4]) > 60:
                dic[l[0]].append(l[3])
                dic[l[0]].append(l[4])

    return dic


def Run_Phobius(filefasta):
    dicSignal = defaultdict(list)
    dicTransm = defaultdict(list)
    cmd = './phobius/phobius ' + filefasta + ' > ' + filefasta + '_phobius.out'

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()

    out = '' + filefasta + '_phobius.out'

    f = open(out, 'r')
    for l in f:
        if 'ID' in l:
            l = re.sub(r"\s+", ' ', l)
            l = l.split(' ')
            key = l[1]
        elif 'SIGNAL' in l:
            l = re.sub(r"\s+", ' ', l)
            l = l.split(' ')
            dicSignal[key].append(l[2:4])
        elif 'TRANSMEM' in l:
            l = re.sub(r"\s+", ' ', l)
            l = l.split(' ')
            dicTransm[key].append(l[2:4])
    return dicSignal, dicTransm


""" ***************************************************************************************************************** """
fasta = {}
with open(sys.argv[1], "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        fasta[record.id] = record.seq

""" ****************************** Signal peptide and transmembrane sement ***************************************** """
signalp = Run_Signalp(sys.argv[1])
tmhmm = Run_Tmhmm(sys.argv[1])
phobiusS, phobiusT = Run_Phobius(sys.argv[1])

""" ************************** Feature extraction ************************* """
datasetFeatureCS = defaultdict(dict)
datasetFeatureSP = defaultdict(dict)
datasetFeatureTM = defaultdict(dict)

for key in fasta:
    """ ************************** Complete sequence ************************* """
    # Features for complete sequence
    featureCS = defaultdict(list)
    featureSP = defaultdict(list)
    featureTM = defaultdict(list)

    # Tripeptide
    aux = list()
    res = CalculateTripeptide(fasta[key])
    for r in res:
        aux.append(res[r])
    featureCS['Tripeptide'].append(aux)

    # Dipeptide
    aux = list()
    res = CalculateDipeptide(fasta[key])
    for r in res:
        aux.append(res[r])
    featureCS['Dipeptide'].append(aux)

    # AAComposition
    aux = list()
    res = CalculateAAComposition(fasta[key])
    for r in res:
        aux.append(res[r])
    featureCS['AAComposition'].append(aux)

    # AAComposition_N_C
    medium = int(len(fasta[key]) / 2)
    aux = list()
    res = CalculateAAComposition(fasta[key][0:medium])
    for r in res:
        aux.append(res[r])

    res = CalculateAAComposition(fasta[key][medium::])
    for r in res:
        aux.append(res[r])

    featureCS['AA_Composition_N_C'].append(aux)

    # CPAASC (chemical properties of amino acid side chains)
    ef = FeacturesExtraction(fasta[key])
    ef.calculate_proportion()
    feat = ef.get_all_features()
    featureCS['CPAASC'].append(feat)

    # CPAASC_N_C
    medium = int(len(fasta[key]) / 2)
    ef = FeacturesExtraction(fasta[key][0:medium])
    ef.calculate_proportion()
    feat = ef.get_all_features()
    aux = list()
    for f in feat:
        aux.append(f)

    ef = FeacturesExtraction(fasta[key][medium::])
    ef.calculate_proportion()
    feat = ef.get_all_features()
    for f in feat:
        aux.append(f)

    featureCS['CP_AASC_N_C'].append(aux)

    datasetFeatureCS[key] = featureCS

""" **************** Load the model of the disk ******************************* """
models = defaultdict(list)
model_f = open(sys.argv[2], 'r')

for f in sorted(model_f):
    f = f.strip().split('\t')
    models[f[1]].append(f[0])

""" ******************************* First layer ******************************* """

fl = defaultdict(list)
for k1 in sorted(datasetFeatureCS):
    for k2 in sorted(datasetFeatureCS[k1]):
        for k3 in models:
            if 'RLPs_vs_noRLPs' in k3:
                for m in models[k3]:
                    if k2 in m:
                        p = list()
                        loaded_model = pickle.load(open(str(m), 'rb'))
                        X = np.asarray(datasetFeatureCS[k1][k2])
                        p = loaded_model.predict_proba(X)
                        pred = loaded_model.predict(X)
                        fl[k1].append([round(p[0][pred[0]], 3), pred[0]])
                        # print (k1 +'\t'+ k2 +'\t'+ k3 +'\t'+ str(m) +'\t'+ str(p[0][0])+'\t'+ str(p[0][1]) +'\t'+ str(pred[0]) )

""" ******************************* Second layer ******************************* """

sl = defaultdict(list)
for k1 in sorted(datasetFeatureCS):
    for k2 in sorted(datasetFeatureCS[k1]):
        for k3 in models:
            if 'RLPs_vs_RLKs' in k3:
                for m in models[k3]:
                    if k2 in m:
                        loaded_model = pickle.load(open(str(m), 'rb'))
                        X = np.asarray(datasetFeatureCS[k1][k2])
                        p = loaded_model.predict_proba(X)
                        pred = loaded_model.predict(X)
                        sl[k1].append([round(p[0][pred[0]], 3), pred[0]])
                        # print (k1 +'\t'+ k2 +'\t'+ k3 +'\t'+ str(m) +'\t'+ str(p[0][0])+'\t'+ str(p[0][1]) +'\t'+ str(pred[0]) )

""" ******************************* Third layer ******************************* """

tl = defaultdict(list)
for k1 in sorted(datasetFeatureCS):
    for k2 in sorted(datasetFeatureCS[k1]):
        for k3 in models:
            if 'RLPs_Subfamily' in k3:
                for m in models[k3]:
                    if k2 in m:
                        loaded_model = pickle.load(open(str(m), 'rb'))
                        X = np.asarray(datasetFeatureCS[k1][k2])
                        pred = loaded_model.predict(np.array(X))
                        p = loaded_model.predict_proba(np.array(X))
                        tl[k1].append([round(max(p[0]), 3), pred[0]])

""" ******************************************  Bayesian evolutional classification ************************************************************** """
compilated_layer = {}  # Armazena todo resultado
decision_vector = dict()
output_file = '' + sys.argv[1] + '_predictions.tsv'
output_file = open(output_file, 'w')
output_file.write('Accession\tSignal peptide\tTransmembrane segment\tRLP-noRLP\tRLP-noRLP Probability\tRLP-RLK\tRLP-RLK Probability\tRLP-Subfamily\tRLP-Subfamily Probability\tClassification\tDecision probability\tComplete\n')

for k in sorted(fl):
    # phobiusS,phobiusT
    print('************************************************************************ ' +
          k + ' ************************************************************************')
    decision = [0, 0, 0, 0, 0]
    res = list()
    index_name = ''
    if signalp[k][8] == 'Y':
        res.append(signalp[k][4])
        decision[0] = 1
    elif len(phobiusS[k]) > 0:
        res.append(phobiusS[k][0])
        decision[0] = 1
    else:
        res.append('No signal peptide')
        decision[0] = 0

    # TMHMM
    if len(tmhmm[k]) > 1:
        res.append(tmhmm[k][0])
        decision[1] = 1
    elif len(phobiusT[k]) > 0:
        res.append(phobiusT[k][0])
        decision[1] = 1
    else:
        res.append('No TM')
        decision[1] = 0
    evol1 = [0, 0]
    """ ******************************************  Transmembrane Topology and Signal Peptide positions ******************************************* """
    if signalp[k][8] == 'Y':
        index_name = '__SignalP: 1-' + signalp[k][1] + ' p:' + res[0]
        evol1[0] = 1
    if len(phobiusS[k]) > 0:
        index_name += '__SignalPhobius: ' + \
            phobiusS[k][0][0] + '-' + phobiusS[k][0][1]
        evol1[0] = 1
    if len(tmhmm[k]) > 1:
        index_name += '__TMHMM: '
        index_name += tmhmm[k][1] + '-' + tmhmm[k][2]
        evol1[1] = 1
    if len(phobiusT[k]) > 0:
        index_name += '__phobiusT: '
        for t in phobiusT[k]:
            index_name += t[0] + '-' + t[1]
            evol1[1] = 1

    """ *************************************************** Integration of the first layer ******************************************************** """
    print('\nIntegrating the first layer.\n')
    n = len(fl[k])
    h = 0
    for l in fl[k]:
        if l[1] == 0:
            h += 1

    # Update of information about SP and TM
    if evol1[1] == 1 or sum(evol1) == 2:
        n += 1
        h += 1
    else:
        n += 1

    # Alpha e beta
    alpha = 0.01
    beta = 0.01

    with pm.Model() as model:
        p = pm.Beta('RLP-noRLK', alpha=alpha, beta=beta)
        y = pm.Binomial('y', n=n, p=p, observed=h)
        trace = pm.sample(2000, tune=1000, cores=4, init='advi')

    RLP_noRLP = pm.summary(trace)
    RLP_noRLP.index = ['RLP-noRLP']

    if RLP_noRLP['mean'][0] > 0.6 and decision[1] != 0:
        decision[2] = 1
    else:
        decision[2] = 0
        decision[3] = 0
        decision[4] = 0

    """ ********************************************* Integration of the second layer ************************************************************** """
    print('\nIntegrating the second layer.\n')
    n = len(sl[k])
    h = 0
    for l in sl[k]:
        if l[1] == 0:
            h += 1
    # Update of information of the first layer
    if decision[2] == 1:
        n += 1
        h += 1
    else:
        n += 1

    # Alpha e beta
    alpha = 0.01
    beta = 0.01

    with pm.Model() as model:
        p = pm.Beta('RLPs-RLKs', alpha=alpha, beta=beta)
        y = pm.Binomial('y', n=n, p=p, observed=h)
        trace = pm.sample(2000, tune=1000, cores=4, init='advi')

    RLP_RLK = pm.summary(trace)
    RLP_RLK.index = ['RLP-RLK']

    if RLP_RLK['mean'][0] > 0.6 and decision[2] != 0:
        decision[3] = 1
    else:
        decision[3] = 0
        decision[4] = 0
    """ ****************************************** Integration of the third layer ******************************************************************* """
    print('\nIntegrating the third layer.\n')

    RLPclass = {'L-Lectin-RLK': 0, 'LRR-RLK': 1, 'S-domain-RLK': 2, 'Malectin-RLK': 3, 'Salt-stress-response/antifungal-RLK': 4,
                'WAK-RLK': 5, 'B-Lectin-RLK': 6, 'Unknown-RLK': 7, 'Other-RLK': 8, 'LysM-RLK': 9, 'PAN-RLK': 10, 'Ethylene-responsive-RLK': 11,
                'Thaumatin-RLK': 12, 'RCC1-RLK': 13, 'Glycosyl-hydrolases-RLK': 14, 'C-Lectin-RLK': 15, 'Unknown': 16}

    pred_class = {}
    for t in tl[k]:
        if t[0] > 0.3:  # Alterado
            t[1] = t[1].replace('RLK', 'RLP')
            pred_class[t[1]] = 0
        else:
            pred_class['Unknown'] = 0

    for t in tl[k]:
        if t[0] > 0.3:  # Alterado
            t[1] = t[1].replace('RLK', 'RLP')
            pred_class[t[1]] += 1
        else:
            pred_class['Unknown'] += 1

    lv = len(pred_class.keys())
    vector = np.full((1, lv), 0)
    vector_class = list()
    index = 0
    for class_p in pred_class:
        vector[0][index] = pred_class[class_p]
        vector_class.append(class_p)
        index += 1

    with pm.Model() as model:
        theta = pm.Dirichlet('theta', a=vector, shape=vector.shape)
        post = pm.Multinomial('post', n=vector.sum(), p=theta, observed=vector)
        trace = pm.sample(2000, tune=1000, cores=4, init='advi')

    RLP_aux = pm.summary(trace)
    RLP_aux.index = vector_class
    RLP_Sub = RLP_aux

    maxid = RLP_Sub.idxmax()['mean']
    for index, row in RLP_aux.iterrows():
        if index == maxid:
            continue
        else:
            RLP_Sub = RLP_Sub.drop(index, axis=0)

    if RLP_Sub['mean'][0] > 0.30  and decision[3] != 0:
        decision[4] = 1
    else:
        decision[4] = 0

    """ ********************************************* Bayesian integration decision ************************************************************** """
    print('\nBayesian integration decision:\n')
    n = len(decision)
    h = sum(decision)

    alpha = 0.01
    beta = 0.01

    with pm.Model() as model:
        p = pm.Beta('Decision', alpha=alpha, beta=beta)
        y = pm.Binomial('y', n=n, p=p, observed=h)
        trace = pm.sample(2000, tune=1000, cores=4, init='advi')

    Desc = pm.summary(trace)
    Desc.index = ['(' + RLP_Sub.index[0] + ')']
    RLPs_Pred = pd.concat([RLP_noRLP, RLP_RLK, RLP_Sub, Desc])
    "-----------------------------------------------------------------------------"
    RLPs_Pred.index.name = index_name

    aux = RLPs_Pred.index.name.split('__')
    signal_peptide = '-'
    Transm_segment = '-'
    RLP_noRLP = ''
    RLP_RLK = ''
    RLP = ''
    descaux = ''
    namedec = ''
    name = ''
    hpd1 = ''
    hpd2 = ''
    sd = ''

    for ax in aux:
        if 'SignalP' in ax or 'SignalPhobius' in ax:
            signal_peptide += ' ' + ax
        elif 'TMHMM' in ax or 'phobiusT' in ax:
            Transm_segment += ax + ' '

    for index, row in RLPs_Pred.iterrows():
        if 'RLP-noRLP' in str(index):
            RLP_noRLP = str(round(row['mean'], 4))
            if row['mean'] > 0.6:
                RLP_noRLP = 'RLP\t' + RLP_noRLP
            else:
                RLP_noRLP = 'noRLP\t' + RLP_noRLP

        elif 'RLP-RLK' in str(index):
            RLP_RLK = str(round(row['mean'], 4))
            if row['mean'] > 0.6:
                RLP_RLK = 'RLP\t' + RLP_RLK
            else:
                RLP_RLK = 'RLK-like\t' + RLP_RLK
        elif '(' in str(index):
            decaux = str(round(row['mean'], 4))

            if row['mean'] > 0.70:
                namedec = index
            else:
                namedec = 'No RLP'
        else:
            RLP = str(round(row['mean'], 4))
            if row['mean'] > 0.3:
                RLP = index +'\t'+ RLP
            else:
                RLP = 'Unknown\t' + RLP

    if '-' == signal_peptide:
        signal_peptide = 'No SP'
    if '-' == Transm_segment:
        Transm_segment = 'No TM'

    out = str(k + '\t' + str(signal_peptide) + '\t'+str(Transm_segment) + '\t'+str(RLP_noRLP) +
              '\t'+str(RLP_RLK)+'\t'+str(RLP)+'\t'+namedec+'\t'+str(decaux)+'\tOk\n')
    output_file.write(out)

output_file.close()

# https://ericmjl.github.io/bayesian-stats-talk/
# https://juanitorduz.github.io/intro_pymc3/
# https://towardsdatascience.com/hands-on-bayesian-statistics-with-python-pymc3-arviz-499db9a59501if
