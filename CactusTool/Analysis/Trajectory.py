from re import A
from .DiscreteFunction import VectorSeries
from ..Lib.EccRed import tCoal0PN, m_omega_r35SPN, make_ansatz1, make_ansatz2a, make_ansatz2, lambda_r1PNResSum, lambda_t1PNResSum, Pr35PN, delta_rmas
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

class Trajectory:

    def __init__(self, t, x1, y1, z1, x2, y2, z2):
        self.t = t
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.x2 = x2
        self.y2 = y2
        self.z2 = z2

    @property
    def r(self):
        return VectorSeries(self.t, np.stack([self.x1-self.x2, self.y1-self.y2, self.z1-self.z2], axis=1))

    @property
    def init_separation(self):
        return self.r.y[0][0]
    
    def tmerger(self, eps=1e-2):
        index = self.r.norm.NearestValueIndex(eps)
        return self.t[index]

    @property
    def omega(self):
        r = self.r.clip(self.t[0], self.tmerger(1e-5))
        v = r.gradient()
        r_norm = r.norm
        n = len(r)
        omega = np.zeros((n, 3))
        for i in range(n):
            omega[i] = np.cross(r[i], v[i])/(r_norm[i]**2)
        return VectorSeries(r.t, omega)

    def e_Omega(self, ADM_Mass, eta, chi_1, chi_2, MinTime=100, MaxTime=800):
        q = (1 + np.sqrt(1 - 4*eta) - 2*eta)/(2*eta)
        D0 = self.init_separation
        s1z = chi_1[2]
        s2z = chi_2[2]

        Omega_0 = m_omega_r35SPN(D0,q,chi_1,chi_2,1.) if m_omega_r35SPN(D0,q,chi_1,chi_2,1.) > 10e-6 else 0
        omega_data = self.omega.norm.clip(MinTime, max(MaxTime, 2/3 * self.tmerger()))

        MAXTime = tCoal0PN(eta=eta, total_mass=ADM_Mass, separation=D0)
        popt1, pcov1 = curve_fit(make_ansatz1(eta, s1z, s2z, MAXTime),omega_data.t, omega_data.y)
        a, t_0 = popt1

        # first fit only for e and phase offset t_1
        popt2a, pcov2a = curve_fit(make_ansatz2a(eta, s1z, s2z, MAXTime, Omega_0, a, t_0, omega_1=1., ke=0., k_omega=0.), omega_data.t, omega_data.y, p0=[1e-2,  -np.pi], bounds=([1e-5,-2*np.pi], [5e-2,+2*np.pi]))
        e, t_1 = popt2a

        # # Fitting with ansatz1 params not fixed
        popt2, pcov2 = curve_fit(make_ansatz2(eta, s1z, s2z, MAXTime, Omega_0), omega_data.t, omega_data.y, p0=[e, a, t_0, 1.0, t_1, 0., 0.], bounds=([1e-5, 0.5, -5.0, -5.0, -2*np.pi, -10., -10.], [5e-2,3.0,+5.0,+5.0,+2*np.pi,+10.,+10.]))
        e, a, t_0, omega_1, t_1, ke, k_omega = popt2
        perr = np.sqrt(np.diag(pcov2))
        delta_e = perr[0]
        delta_ke = perr[5]
        delta_t_1 = perr[4]
        t = 0.
        amp = e*(1-.0001*ke*t)
        delta_amp = delta_e*np.abs(1 - .0001*ke*t) + np.abs(e*-.0001*t)*delta_ke
        e_Omega = .5*amp/Omega_0
        delta_e_Omega = .5*delta_amp/Omega_0

        # plot of fitted function
        plt.plot(omega_data.t, omega_data.y, label="data")
        plt.plot(omega_data.t, make_ansatz1(eta, s1z, s2z, MAXTime)(omega_data.t, *popt1), linestyle="dashed", label="first ansatz")
        plt.plot(omega_data.t, make_ansatz2(eta, s1z, s2z, MAXTime, Omega_0)(omega_data.t, *popt2), linestyle="dashed", label="second ansatz")
        plt.legend()
        plt.xlabel('t [M]')
        plt.ylabel(r'orbital frequency $\Omega$')
        msg = "ecc estimate: %g +/- %g" % (e_Omega, delta_e_Omega)
        plt.title(msg)

        self.lambda_R = 1/lambda_r1PNResSum(amp, t_1, Pr35PN(D0,q,chi_1,chi_2,1.), m_omega_r35SPN(D0,q,chi_1,chi_2,1.),D0,eta,1.) #radial_factor
        self.lambda_T = 1/lambda_t1PNResSum(amp, t_1, m_omega_r35SPN(D0,q,chi_1,chi_2,1.),D0,eta,1.) #tangential factor
        self.delta_R = delta_rmas(amp,D0,m_omega_r35SPN(D0,q,chi_1,chi_2,1.),q,1.) #radius correction

    # def eccentricity(self, p0=[0, 0, 0.01, 0.002, 0]):
    #     """
    #     $$
    #     \dot{D}(t)=A_{0}+A_{1} t-e D_{0} \omega_{e} \sin \left(\omega_{e} t+\phi_{e}\right)
    #     $$
    #     where $e$ is the eccentricity and $D_{0}$ the initial coordinate interbinary distance.

    #     arXiv:1605.03424
    #     """
    #     def orbital_evolution(t, A0, A1, e, w, phi):
    #         return A0 + A1*t - e*self.init_separation*w*np.sin(w*t + phi)
        
    #     # The fit is performed in the time interval between $t_{\mathrm{ret}}=50 \mathrm{M}$ and $t_{\mathrm{ret}}=\frac{2}{3} t_{\mathrm{merger}}$ to avoid the initial spurious radiation and the plunge phase but having at least one eccentricity cycle included.
    #     ddot = self.r.norm.gradient().clip(50, 2/3 * self.tmerger())
    #     t = ddot.t
    #     y = ddot.y
        
    #     params, params_covariance = curve_fit(orbital_evolution, t, y, p0=p0)
    #     print('Orbital Eccentricity =', params[2])

    #     fig, ax = plt.subplots()
    #     ax.plot(t, y)
    #     ax.set_xlabel('time [M]')
    #     ax.set_ylabel(r'$\dot{D}$')
    #     ax.plot(t, orbital_evolution(t, params[0], params[1], params[2], params[3], params[4]), '--', label='Fitted function')
    #     plt.legend(loc='best')
    #     plt.show()

    def preview(self):
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")

        fig = go.Figure(data=go.Scatter3d(
            x=self.x1-self.x2, 
            y=self.y1-self.y2, 
            z=self.z1-self.z2,
            marker=dict(
                size=3,
                colorscale='Viridis',
            ),
            line=dict(
                width=2
            )
        ))
        fig.show()
