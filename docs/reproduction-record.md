# Reproduction record

22 November 2026

The source-only release was tested in an isolated Python 3.12 environment with the dependency versions in `requirements.txt`: matplotlib 3.11.0, numpy 2.4.5, pandas 3.0.5 and scipy 1.18.1.

| Execution | Outcome |
| --- | --- |
| `python -m unittest -v test_depth_ghg_model.py` | All nine tests passed |
| `python depth_ghg_model.py --no-plot` | All 23 scenarios completed and produced summary and time-series CSV files |
| `python figures/generate_current_figures.py` | All three figures exported successfully |
| Figure layout checks | Every label remained inside its canvas; all three figures passed the text-collision checks |
| Figure 1 geometry checks | Depth fields and transects retained their reference values; area-fraction histograms were normalized |
| Figure 3 curve checks | All eight stored arrays matched the reference definitions; panel d had non-negative excess storage and equal cumulative release within numerical tolerance |

The numerical program, tests, parameter and scenario tables, pinned dependencies and licence retain the bytes of the reviewed model archive. The public release adds the adopted figure sources, repository documentation and citation metadata. Generated execution outputs can be reproduced with the commands above.

The model checks establish numerical and prescribed directional behaviour. The experiments and independent budgets in `evaluation-plan.md` provide the design for evaluating field responses.
