import numpy as np
import matplotlib.pyplot as plt

# Load data (replace filename if needed)
data = np.loadtxt('coulcorrcompare.betatest.txt')

# Extract columns
beta_T = data[:, 0]
avg_abs_diff = data[:, 1]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(beta_T, avg_abs_diff, 'o-', markersize=8)

plt.xlabel(r'$\beta_T$', fontsize=18)
plt.ylabel(r'$\langle\,|\mathrm{rel.\ diff.}|\,\rangle$', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.grid(True)

# Add parameter labels
param_text = (
    r'$\alpha = 1.2$' + '\n' +
    r'$R_{\mathrm{out}} = 5.7$ fm' + '\n' +
    r'$R_{\mathrm{side}} = 5.0$ fm' + '\n' +
    r'$R_{\mathrm{long}} = 6.8$ fm' + '\n' +
    r'$\lambda = 0.8$'
)
plt.text(0.3, 0.9, param_text, transform=plt.gca().transAxes,
         fontsize=18, verticalalignment='top')

plt.tight_layout()
plt.savefig("coulcorrcompare.betatest.png")
