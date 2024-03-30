import numpy as np

def jaccard(set1, set2):
    intersection = len(set(set1).intersection(set2))
    union = (len(set1) + len(set2)) - intersection
    if union == 0: return 1
    return float(intersection) / union