import numpy as np
import mean_field as mf
import scipy.linalg as lin
import scipy.optimize as opt
import minimzation_T1 as mt1
import matplotlib.pyplot as plt

def main():
    """ MAIN """
    # Prepare parameters
    Ns = 4
    kappa = 1.0
    Q_minimal, mu_minimal = mt1.find_q_mu(kappa, Ns)

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
                                            'T1')[0] - H/2)
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
                                      'T1')
        B = np.identity(dim)
        B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)

        w, v = lin.eig(np.dot(B, M))
        # print 'condensate:', v[0]

        # normalization factor
        r = np.abs(v[0][0])**2 + np.abs(v[0][1])**2      

        result = (np.abs(v[0][0])**2 - np.abs(v[0][1])**2)/r
        magnetization_values.append(result)

        # only one condensate 
        # computing z component of magnetization?
        print 'Result:', result 
        
    plt.plot(H_values, magnetization_values)
    plt.xlabel('H')
    plt.ylabel('S_z')
    plt.savefig('magnetization_triangular')
    plt.show()

if __name__ == '__main__':
    main()
