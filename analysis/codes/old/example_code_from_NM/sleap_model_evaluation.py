# get model evaluation
# from https://sleap.ai/notebooks/Model_evaluation.html

import sleap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.style.use("seaborn-deep")
sleap.versions()

help(sleap.load_metrics)

model = "F:/node17/Rspe-tandem/models/221219_125502.centered_instance.n=353"
metrics = sleap.load_metrics(model, split="val")
print("\n".join(metrics.keys()))

print("Error distance (50%):", metrics["dist.p50"])
print("Error distance (90%):", metrics["dist.p90"])
print("Error distance (95%):", metrics["dist.p95"])

plt.figure(figsize=(6, 3), dpi=150, facecolor="w")
sns.histplot(metrics["dist.dists"].flatten(), binrange=(0, 20), kde=True, kde_kws={"clip": (0, 20)}, stat="probability")
plt.xlabel("Localization error (px)");
plt.show()


plt.figure(figsize=(6, 3), dpi=150, facecolor="w")
sns.histplot(metrics["oks_voc.match_scores"].flatten(), binrange=(0, 1), kde=True, kde_kws={"clip": (0, 1)}, stat="probability")
plt.xlabel("Object Keypoint Similarity");
plt.show()

#An easy way to s ummarize this analysis is to take the average over all of these thresholds to compute the mean Average Precision (mAP) and mean Average Recall (mAR) which are widely used in the pose estimation literature.
print("mAP:", metrics["oks_voc.mAP"])
print("mAR:", metrics["oks_voc.mAR"])