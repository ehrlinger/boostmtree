# Spirometry Data

Data consists of 9471 longitudinal evaluations of forced 1-second
expiratory volume (FEV1-percentage of predicted) after lung transplant
from 509 patients who underwent lung transplant (LTx) at the Cleveland
Clinic. Twenty three patient/procedure variables were collected at the
time of the transplant. The major objectives are to evaluate the
temporal trend of FEV1 after LTx, and to identify factors associated
with post-LTx FEV1 and assessing the differences in the trends after
Single LTx versus Double LTx.

## Format

A list containing four elements:

1.  The 23 patient variables (features).

2.  Time points (time).

3.  Unique patient identifier (id).

4.  FEV1-outcomes (y).

## References

Mason D.P., Rajeswaran J., Li L., Murthy S.C., Su J.W., Pettersson G.B.,
Blackstone E.H. Effect of changes in postoperative spirometry on
survival after lung transplantation. *J. Thorac. Cardiovasc. Surg.*,
144:197-203, 2012.

## Examples

``` r
data(spirometry, package = "boostmtree")
```
