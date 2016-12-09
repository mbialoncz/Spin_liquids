import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import mean_field as mf
import minimization as mm
import scipy.linalg as lin
import scipy.optimize as opt
import matplotlib.pyplot as plt

def main():
    """ Check everything """
    ansatzs = ['T1', 'K01']
    Nss = [8 + it for it in range(4)]
    kappas = [1.0 + 0.2 * it for it in range(4)]
    for ansatz in ansatzs:
        for Ns in Nss:
            for kappa in kappas:
                print 'Current parameters:', ansatz, Ns, kappa
                calculate(ansatz, kappa, Ns)

def calculate(ansatz, kappa, Ns):
    """ MAIN """
    # Prepare parameters
    Q_minimal, mu_minimal = mm.find_q_mu(kappa, Ns, ansatz)

    H_values = np.arange(0, 2, 0.01)
    magnetization_values = []

    for H in H_values:
        #find the wave vector corresponding to zero mode of energy
        g = lambda k : np.abs(mf.dispersion(k[0],
                                            k[1],
                                            [Q_minimal,0],
                                            [0,0],
                                            0,
                                            0,
                                            mu_minimal,
                                            kappa,
                                            Ns,
                                            ansatz)[0] - H/2)
        sol = opt.minimize(g, [-2.,2.], method = 'Nelder-Mead')
        k_minim = sol.x
        n1 = Ns * k_minim[0]/(2 * np.pi)
        n2 = Ns * k_minim[1]/(2 * np.pi)

        M, dim = mf.HamiltonianMatrix(n1,
                                      n2,
                                      [Q_minimal,0],
                                      [0,0],
                                      0,
                                      0,
                                      mu_minimal,
                                      kappa,
                                      Ns,
                                      ansatz)
        B = np.identity(dim)
        B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)

        w, v = lin.eig(np.dot(B, M))
        # print 'condensate:', v[0]

        # normalization factor
        r = np.abs(v[0][0])**2 + np.abs(v[0][1])**2  

        result = (np.abs(v[0][0])**2 - np.abs(v[0][1])**2)/r
        magnetization_values.append(result)

    title = 'ansats-{}-kappa-{}-Ns-{}'.format(ansatz, kappa, Ns)
    savepath = 'tmp/' + title + '.png'
    plt.plot(H_values, magnetization_values, 'm--o', lw=3, alpha=0.666)
    plt.xlabel('H')
    plt.ylabel('S_z')
    plt.title(title)
    plt.savefig(savepath)
    plt.show()
    plt.clf()

if __name__ == '__main__':
    main()
