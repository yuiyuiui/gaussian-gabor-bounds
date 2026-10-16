import json
from gaussian_gabor import GaussianPlan, choose_resolution

resolution = choose_resolution(1.5, x=0.2, y=1.1, analytic_error=0.05)
result = GaussianPlan(1.5, x=0.2, y=1.1, resolution=resolution).bounds()
assert sum(result["errors"].values()) <= 0.05
assert result["configuration"]["roundoff_certified"] is False
print(json.dumps(result, indent=2))
