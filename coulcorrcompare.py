import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("coulcorrcompare.out", comments="#")
q = data[:, 0]
full = [data[:, 1], data[:, 3], data[:, 5]]
approx = [data[:, 2], data[:, 4], data[:, 6]]

labels = [r'$q_{\mathrm{out}}$ [GeV/$c$]', r'$q_{\mathrm{side}}$ [GeV/$c$]', r'$q_{\mathrm{long}}$ [GeV/$c$]']
param_text = (
    r'$\alpha = 1.2$' + '\n' +
    r'$R_{\mathrm{out}} = 5.7$ fm' + '\n' +
    r'$R_{\mathrm{side}} = 5.0$ fm' + '\n' +
    r'$R_{\mathrm{long}} = 6.8$ fm' + '\n' +
    r'$\lambda = 0.8$, $\beta_T = 0.9$'
)

fig, axs = plt.subplots(1, 3, figsize=(12, 4), sharey=True)

for i, ax in enumerate(axs):
    ax.plot(q, full[i], 'o-', label='Full', markersize=3)
    ax.plot(q, approx[i], '-', label='Approx', linewidth=1.5)
    ax.set_xlim(0, 0.1)
    ax.set_xticks([0.0, 0.02, 0.04, 0.06, 0.08])
    ax.set_xlabel(labels[i], fontsize=18)
    ax.tick_params(axis='both', labelsize=14)

    if i == 0:
        ax.set_ylabel(r'$C(\mathbf{q})$', fontsize=18)
    #else:
    #    ax.set_yticklabels([])

    if i == 1:
        ax.text(0.20, 0.60, param_text, transform=ax.transAxes,
                fontsize=14, va='top')

fig.subplots_adjust(left=0.08, right=1.00, top=0.95, bottom=0.15, wspace=0.00)
axs[0].legend(fontsize=14, loc='lower right')
plt.savefig("coulcorrcompare.png", dpi=300)
plt.show()
