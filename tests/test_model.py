import unittest, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'content'))
import numpy as np
from teaching_model import Parameters,simulate,preflight,transport_slope_sensitivity,hydraulic_transport,classroom_estimate
class Tests(unittest.TestCase):
 def test_equilibrium(self):
  r=simulate(Parameters(duration_s=300)); self.assertLess(np.max(np.abs(r['bed_m'][-1]-r['bed_m'][0])),1e-11)
 def test_aggradation(self):
  r=simulate(Parameters(feed_ratio=1.3,duration_s=600)); self.assertGreater(r['storage_change_m3'],0);self.assertLess(abs(r['mass_balance_error_m3']),1e-8)
 def test_degradation(self):
  r=simulate(Parameters(feed_ratio=.7,duration_s=600)); self.assertLess(r['storage_change_m3'],0);self.assertLess(abs(r['mass_balance_error_m3']),1e-8)
 def test_zero_transport(self):
  r=simulate(Parameters(discharge_m3s=.1,duration_s=100));self.assertTrue(np.all(r['transport_faces_m2s']==0))
 def test_analytic_derivative(self):
  p=Parameters();s=np.array([.003,.008,.02]);eps=1e-6
  _,_,qp=hydraulic_transport(s+eps,p);_,_,qm=hydraulic_transport(s-eps,p)
  np.testing.assert_allclose(transport_slope_sensitivity(s,p),(qp-qm)/(2*eps),rtol=2e-4,atol=1e-10)
 def test_preflight(self):
  d=classroom_estimate();self.assertGreater(d['cells'],0);self.assertGreater(d['planned_dt_s'],0)
 def test_limit(self):
  with self.assertRaises(ValueError):simulate(Parameters(dx_m=.1,duration_s=0))
if __name__=='__main__':unittest.main(verbosity=2)
