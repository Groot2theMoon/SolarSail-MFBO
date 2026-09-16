"""
mfbo.py 코드용 mf - currin exponential 2d function 데이터 생성 코드.
currin function에 대해서는 아래 코드 참조.
https://www.sfu.ca/~ssurjano/curretal88exp.html
"""
import numpy as np

def evaluate_currin_mf(x, s):
    x = np.atleast_2d(x)

    x1 = x[:, 0]
    x2 = x[:, 1]

    def f_high(x1_val, x2_val):
        x2_val = np.clip(x2_val, 1e-12, 1.0)
        term1 = 1 - np.exp(-1 / (2 * x2_val))
        num = 2300 * x1_val**3 + 1900 * x1_val**2 + 2092 * x1_val + 60
        den = 100 * x1_val**3 + 500 * x1_val**2 + 4 * x1_val + 20
        return term1 * (num / den)

    if s > 0.5:
        # High Fidelity (s=1)
        return f_high(x1, x2)
    else:
        # Low Fidelity (s=0)
        val_lf = 0.25 * (
            f_high(x1 + 0.05, x2 + 0.05) +
            f_high(x1 + 0.05, x2 - 0.05) +
            f_high(x1 - 0.05, x2 + 0.05) +
            f_high(x1 - 0.05, x2 - 0.05)
        )
        return val_lf