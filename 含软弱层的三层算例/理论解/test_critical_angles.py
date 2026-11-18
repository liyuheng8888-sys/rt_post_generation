"""Critical-angle regressions against an independent displacement/traction ODE.

Run: python -m unittest discover -s <this directory> -p test_critical_angles.py
The ODE reference is used for thin layers, where its transfer matrix is stable.
"""
import copy
import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.linalg import expm

_spec = importlib.util.spec_from_file_location(
    'critical_theory', Path(__file__).with_name('rt_coefficients_theorical_solution.py'))
theory = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(theory)
KEYS = ('R_P', 'R_SV', 'T_P', 'T_SV')


def elastic_matrix(m, p):
    rho = m['rho']; mu = rho*m['vs']**2
    lam = rho*m['vp']**2 - 2*mu; longitudinal = lam + 2*mu
    return np.array([
        [0, -p, 1/mu, 0],
        [-lam*p/longitudinal, 0, 0, 1/longitudinal],
        [rho-4*mu*(lam+mu)*p*p/longitudinal, 0, 0, -lam*p/longitudinal],
        [0, rho, -p, 0]], complex)


def elastic_columns(m, p, sign):
    rho, vp, vs = m['rho'], m['vp'], m['vs']; mu = rho*vs*vs
    # Match the input slowness rounding at exterior grazing modes; changing
    # (1/v)**2 to v**-2 can perturb a zero radicand across its square-root branch.
    qp = sign*np.sqrt(complex((1/vp)**2-p*p))
    qs = sign*np.sqrt(complex((1/vs)**2-p*p))
    return np.array([[vp*p, -vs*qs], [vp*qp, vs*p],
                     [2*mu*vp*p*qp, mu*vs*(p*p-qs*qs)],
                     [vp*(rho-2*mu*p*p), 2*mu*vs*p*qs]], complex)


def reference(model, incident):
    side, wave = incident['side'], incident['wave']
    p = math.sin(math.radians(incident['angle_deg'])) / model[side]['vp' if wave=='P' else 'vs']
    omega = 2*np.pi*incident['frequency_hz']
    impedance = model['upper']['rho']*model['upper']['vp']
    scale = np.array([1, 1, impedance, impedance])
    generator = elastic_matrix(model['middle'], p)*scale[None, :]/scale[:, None]
    transfer = expm(1j*omega*model['middle_thickness_m']*generator)
    upper = elastic_columns(model['upper'], p, -1)/scale[:, None]
    lower = elastic_columns(model['lower'], p, 1)/scale[:, None]
    j = 0 if wave=='P' else 1
    incoming = elastic_columns(model[side], p, 1 if side=='upper' else -1)[:, j]/scale
    if side=='upper':
        matrix, rhs = np.column_stack((transfer@upper, -lower)), -transfer@incoming
    else:
        matrix, rhs = np.column_stack((-lower, transfer@upper)), incoming
    return np.linalg.solve(matrix, rhs)


class CriticalAngles(unittest.TestCase):
    def check_elastic(self, model, incident):
        result = theory.solve(model, incident)
        actual = np.array([result['coefficients'][key] for key in KEYS])
        np.testing.assert_allclose(actual, reference(model, incident), rtol=2e-10, atol=2e-10)
        self.assertLess(abs(result['R_energy']+result['T_energy']-1), 2e-10)
        self.assertLess(abs(result['middle_flux_loss']), 2e-10)
        self.assertLess(result['residual'], 1e-10)
        return result

    def test_soft_layer_p_critical_from_both_sides(self):
        model = copy.deepcopy(theory.MODEL); model['middle']['vp'] = 2400.
        for side in ('upper', 'lower'):
            angle = math.degrees(math.asin(model[side]['vs']/2400))
            for offset in (-1e-3, -1e-7, -1e-11, 0., 1e-11, 1e-7, 1e-3):
                for frequency in (20., 95., 300.):
                    for h in (.01, 1., 10.):
                        with self.subTest(side=side, offset=offset, frequency=frequency, h=h):
                            model['middle_thickness_m'] = h
                            self.check_elastic(model, dict(side=side, wave='SV',
                                angle_deg=angle+offset, frequency_hz=frequency))

    def test_internal_sv_and_no_weak_baseline_critical(self):
        stiff = copy.deepcopy(theory.MODEL); stiff['middle'].update(vp=6000., vs=3600.)
        baseline = copy.deepcopy(theory.MODEL); baseline['middle'] = copy.deepcopy(baseline['lower'])
        for model in (stiff, baseline):
            for side in ('upper', 'lower'):
                for wave, speed_key in (('P', 'vp'), ('SV', 'vs')):
                    speed = model[side][speed_key]
                    for middle_key in ('vp', 'vs'):
                        target = model['middle'][middle_key]
                        if target <= speed:continue
                        angle = math.degrees(math.asin(speed/target))
                        for offset in (-1e-7, 0., 1e-7):
                            with self.subTest(side=side, wave=wave, target=target, offset=offset):
                                self.check_elastic(model, dict(side=side, wave=wave,
                                    angle_deg=angle+offset, frequency_hz=95.))

    def test_basis_transition_and_internal_field(self):
        model = copy.deepcopy(theory.MODEL); model['middle']['vp'] = 2400.
        top = model['upper_thickness_m']; h = model['middle_thickness_m']
        omega = 2*np.pi*95
        for sign in (-1, 1):
            for eta in (0., .999e-3, 1.001e-3):
                p = math.sqrt(1+sign*eta*eta)/2400
                angle = math.degrees(math.asin(p*model['upper']['vs']))
                result = self.check_elastic(model, dict(side='upper', wave='SV',
                    angle_deg=angle, frequency_hz=95.))
                # Independently propagate the actual internal state to interior depths.
                scale = np.array([1, 1, 6e6, 6e6]); indices = [0, 1, 3, 4]
                def field(depth):
                    return sum(entry['coefficient']*theory.basis_state(entry, omega, depth)
                               for entry in result['internal'])[indices]/scale
                generator = elastic_matrix(model['middle'], result['p'])*scale[None, :]/scale[:, None]
                for fraction in (0., .25, .5, .75, 1.):
                    expected = expm(1j*omega*fraction*h*generator)@field(top)
                    np.testing.assert_allclose(field(top+fraction*h), expected, atol=2e-10, rtol=2e-10)
                if eta==0:
                    self.assertEqual(sum(entry.get('basis') in ('even', 'odd')
                                         for entry in result['internal']), 2)
                    self.assertLess(result['condition_scaled'], 100.)

    def test_acoustic_internal_critical(self):
        model = copy.deepcopy(theory.MODEL)
        for layer, vp, rho in (('upper', 1500., 1000.), ('middle', 2400., 1200.), ('lower', 1700., 1100.)):
            model[layer].update(type='fluid', vp=vp, vs=0., rho=rho, poroelastic_enabled=False)
        critical = math.degrees(math.asin(1500/2400))
        for offset in (-1e-7, 0., 1e-7):
            incident = dict(side='upper', wave='P', angle_deg=critical+offset, frequency_hz=95.)
            result = theory.solve(model, incident); p = result['p']; omega = 2*np.pi*95
            def column(layer, sign):
                m = model[layer]; q = sign*np.sqrt(complex((1/m['vp'])**2-p*p))
                return np.array([m['vp']*q, m['rho']*m['vp']/1e6])
            m = model['middle']
            generator = np.array([[0, ((1/m['vp'])**2-p*p)*1e6/m['rho']], [m['rho']/1e6, 0]])
            transfer = expm(1j*omega*model['middle_thickness_m']*generator)
            expected = np.linalg.solve(np.column_stack((transfer@column('upper', -1), -column('lower', 1))),
                                       -transfer@column('upper', 1))
            np.testing.assert_allclose([result['coefficients']['R_P'], result['coefficients']['T_P']],
                                       expected, atol=1e-10, rtol=1e-10)
            self.assertLess(abs(result['R_energy']+result['T_energy']-1), 1e-10)

    def test_thick_evanescent_layer_remains_bounded(self):
        model = copy.deepcopy(theory.MODEL)
        model['middle'].update(vp=6000., vs=3600.)
        model['middle_thickness_m'] = 10000.
        incident = dict(side='upper', wave='SV', angle_deg=60., frequency_hz=95.)
        result = theory.solve(model, incident); p = result['p']
        scale = np.array([1, 1, 6e6, 6e6])
        lhs = np.column_stack((elastic_columns(model['upper'], p, -1),
                               -elastic_columns(model['middle'], p, 1)))/scale[:, None]
        rhs = -elastic_columns(model['upper'], p, 1)[:, 1]/scale
        reflected = np.linalg.solve(lhs, rhs)[:2]
        np.testing.assert_allclose([result['coefficients']['R_P'], result['coefficients']['R_SV']],
                                   reflected, atol=1e-10, rtol=1e-10)
        self.assertLess(result['T_energy'], 1e-14)
        self.assertLess(abs(result['R_energy']-1), 1e-10)


if __name__=='__main__':
    unittest.main()
