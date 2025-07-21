import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("coulcorrtest.out", comments="#")
q = data[:, 0]
C_qout = data[:, 1]
C_qside = data[:, 2]
C_qlong = data[:, 3]

# Plotting setup
fig, axs = plt.subplots(1, 3, figsize=(12, 4), sharey=True)

C_list = [C_qout, C_qside, C_qlong]
labels = [r'$q_{\mathrm{out}}$ [GeV/$c$]', r'$q_{\mathrm{side}}$ [GeV/$c$]', r'$q_{\mathrm{long}}$ [GeV/$c$]']
param_text = (
    r'$\alpha = 1.2$' + '\n' +
    r'$R_{\mathrm{out}} = 5.3$ fm' + '\n' +
    r'$R_{\mathrm{side}} = 5.1$ fm' + '\n' +
    r'$R_{\mathrm{long}} = 5.8$ fm' + '\n' +
    r'$\lambda = 0.8$'
)

for i, ax in enumerate(axs):
    ax.plot(q, C_list[i], 'o-', markersize=3)
    ax.set_xlim(0, 0.18)
    ax.set_ylim(bottom=0.95)
    ax.set_xticks([0.0, 0.05, 0.10, 0.15])
    ax.set_xlabel(labels[i], fontsize=18)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(False)

    if i == 0:
        ax.set_ylabel(r'$C(\mathbf{q})$', fontsize=18)
    else:
        ax.set_yticklabels([])

    if i == 1:
        ax.text(0.4, 0.9, param_text, transform=ax.transAxes,
                fontsize=16, va='top')

# Remove space between subplots
fig.subplots_adjust(left=0.08, right=1.00, top=0.95, bottom=0.15, wspace=0.00)

plt.savefig("coulcorrtest.png", dpi=300)
plt.show()
