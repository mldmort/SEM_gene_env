# 02_analysis

Analysis scripts. They read only `data/assembled/gxe_<cohort>.parquet` from
stage 1 and write to `results/` (not tracked by git).

| Folder | Contents |
|---|---|
| `exploratory/` | Quick-look plots (e.g. DNMs vs parental age and parental EA PGS) |
| `sem/` | Structural equation model of the parental age → ASD association (lavaan) |

Load a table with the shared helpers:

```python
import pandas as pd
from gxe import assembled_path

df = pd.read_parquet(assembled_path("gxe_spark_wgs.parquet"))
children = df[df["role"] == "child"]
```
