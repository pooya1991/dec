x = [0] * len(features)
for i in range(len(features)):
    x[i] = features[i][3] - features[i][2]

min_mz = [feature[2] for feature in features]
import matplotlib.pyplot as plt