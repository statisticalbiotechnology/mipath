import pandas as pd
import numpy as np
import warnings
from sklearn.metrics.cluster import expected_mutual_information, mutual_info_score, adjusted_mutual_info_score

def expected_conditional_mutual_information(X,Y,Z):
    joint = pd.crosstab(Z,[X,Y], rownames = ['Z'], colnames = ['X','Y'], dropna=False)
    c = joint.sum(axis=1)
    N = sum(c)
    ecmi=0
    for k in range(len(c)):
        ct = joint.iloc[k].unstack().values
        term1 = c[k]/N
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            term2 = expected_mutual_information(ct, c[k])
        ecmi += term1*term2
    return ecmi

def conditional_mutual_information(X,Y,Z):
    joint = pd.crosstab(Z,[X,Y], rownames = ['Z'], colnames = ['X','Y'], dropna=False)
    c = joint.sum(axis=1)
    N = sum(c)
    cmi=0
    for k in range(len(c)):
        ct = joint.iloc[k].unstack().values
        term1 = c[k]/N
        term2 = mutual_info_score([],[],contingency=ct)
        cmi += term1*term2
    return cmi


def max_conditional_mutual_information(X,Y,Z):
    joint = pd.crosstab(Z,[X,Y], rownames = ['Z'], colnames = ['X','Y'], dropna=False)
    c = joint.sum(axis=1)
    N = sum(c)
    ceX=0
    ceY=0
    for k in range(len(c)):
        ct = joint.iloc[k].unstack().values
        a = np.sum(ct, axis=1)
        b = np.sum(ct, axis=0)
        a = a[a>0]
        b = b[b>0]
        a_sum = np.sum(a)
        b_sum = np.sum(b)
        eA = -np.sum((a / a_sum) * (np.log(a) - np.log(a_sum)))
        eB = -np.sum((b / b_sum) * (np.log(b) - np.log(b_sum)))
        term1 = c[k]/N
        ceX += term1*eA
        ceY += term1*eB
    max_cmi = (ceX+ceY)/2
    return max_cmi

    
def adjusted_contitional_mutual_information(X,Y,Z):
    cmi = conditional_mutual_information(X,Y,Z)
    max_cmi = max_conditional_mutual_information(X,Y,Z)
    ecmi = expected_conditional_mutual_information(X,Y,Z)
    num = cmi - ecmi
    den = max_cmi - ecmi
    acmi = num/den
    return acmi

#Wrapers for now
def adjusted_mutual_information(X,Y):
    ami = adjusted_mutual_info_score(X,Y)
    return ami



