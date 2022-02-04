import numpy as np
import pandas as pd
from scipy.stats import rankdata

def wcov(x,y,mx,my,w):
    return np.sum(w * (x - mx) * (y - my))

def wpearson(x,y,w):
    mx, my = (np.sum(i * w) / np.sum(w) for i in [x, y])
    return wcov(x,y,mx,my,w)/np.sqrt(wcov(x,x,mx,mx,w) * wcov(y,y,my,my,w))

