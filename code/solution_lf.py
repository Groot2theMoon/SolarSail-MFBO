"""
지정한 x1, x2 범위 내에서 전수조사를 실행하는 코드. 
hf 도 실행 가능하나, 보통 lf 용으로 사용했음.
"""
import subprocess
import os
import time
import numpy as np
import torch
import warnings
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
import wandb
from datetime import datetime
import plotly.graph_objects as go
import pandas as pd

warnings.filterwarnings("ignore")

def get_abaqus(new_x, new_s):
    """
    new_x : [x1, x2] (물리적 좌표)
    new_s : 0.0(LF) or 1.0(HF)
    """
    mode_str = "HF" if new_s > 0.5 else "LF"
    x1, x2 = new_x[0], new_x[1]

    command = f"abaqus cae noGUI=run_abaqus.py -- {mode_str} {x1} {x2}"
    print(f"--- Running Abaqus [{mode_str}] x1: {x1:.6f} x2: {x2:.6e} ---")
    
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True 
        )
        
        parsed_data = None
        for line in result.stdout.splitlines():
            if line.strip().startswith("RESULTS:"):
                data_str = line.strip().split("RESULTS:")[1]
                parsed_data = [float(val) for val in data_str.split(',')]
                break
        
        if parsed_data is None:
            print("!!! Error: Could not find 'RESULTS:' tag.")
            print(result.stderr)
            return 1e6, -1e6, -1e6, 1e6

        if mode_str == "HF":
            if len(parsed_data) >= 2:
                print(f"  >> [Done] LF: {parsed_data[0]:.6e}, {parsed_data[1]:.6e}, {parsed_data[2]:.6e} | HF: {parsed_data[3]:.6e}")
                return parsed_data[0], parsed_data[1], parsed_data[2], parsed_data[3]
            else:
                 print("!!! Error: HF result should have at least 2 values.")
                 return 1e6, -1e6, -1e6, 1e6
        else:
            print(f"  >> [Done] LF: {parsed_data[0]:.6e}, {parsed_data[1]:.6e}, {parsed_data[2]:.6e}")
            return parsed_data[0], parsed_data[1], parsed_data[2], None

    except subprocess.CalledProcessError as e:
        print(f"!!! Abaqus execution failed (Return Code {e.returncode})")
        print(e.stderr)
        return 1e6, -1e6, -1e6, 1e6
    except Exception as e:
        print(f"!!! System Error: {e}")
        return 1e6, -1e6, -1e6, 1e6

# 원하는 x1, x2 범위, 간격 지정
X1_RANGE = np.linspace(0.65, 0.95, 4)   # Aspect Ratio: 1.5 ~ 2.5
X2_RANGE = np.exp(np.linspace(np.log(1e-4), np.log(1), 10)) # Thickness: 5e-6 ~ 2e-5

wandb.init(
    entity="dongdong0615-postech",
    # Set the wandb project where this run will be logged.
    project="solarsail-solution-space",
    name=f"grid-lf-solar-sail_{datetime.now().strftime('%Y%m%d_%H%M')}",
)

print(f"--- Starting Grid Search ---")
print(f"x_c (coord) points: {X1_RANGE}")
print(f"d_c (disp) points: {X2_RANGE}")
print(f"Total Simulations: {len(X1_RANGE) * len(X2_RANGE)}")

results_lf_1 = []
results_lf_2 = []
results_lf_3 = []
results_hf = []
design_points = [] # [x1, x2]
loop_time = []

history_x1 = []
history_x2 = []
history_lf_1 = []
history_lf_2 = []
history_lf_3 = []

# --- Grid Loop ---
start_time = time.time()

for i, x1 in enumerate(X1_RANGE):
    for j, x2 in enumerate(X2_RANGE):
        loop_start = time.perf_counter() 
        print(f"\n[Grid Step] Processing x1={x1:.4f}, x2={x2:.4e} ({len(design_points)+1}/{len(X1_RANGE)*len(X2_RANGE)})")
        
        lf_val_1, lf_val_2, lf_val_3, _ = get_abaqus([x1, x2], 0.0)
        
        if not np.isnan(lf_val_1):
            history_x1.append(x1)
            history_x2.append(x2)
            history_lf_1.append(lf_val_1)
            history_lf_2.append(lf_val_2)
            history_lf_3.append(lf_val_3)
        
        loop_duration = time.perf_counter() - loop_start
        loop_time.append(loop_duration)

        fig_3d = go.Figure(data=[go.Scatter3d(
            x=history_x1, y=history_x2, z=history_lf_2,
            mode='markers',
            marker=dict(
                size=5, 
                color=loop_time, 
                colorscale='Viridis',
                colorbar=dict(title="loop duration")
            )
        )])

        fig_3d.update_layout(
            title=f"Solution Space (Step {len(history_x1)})",
            scene=dict(
                xaxis_title='x_c',
                yaxis_title='d_c',
                zaxis_title='y'
            ),
            margin=dict(l=0, r=0, b=0, t=40)
        )

        metrics = {
            "x_c": x1,
            "d_c": x2,
            "LF1": lf_val_1 if lf_val_1 != 1e6 else np.nan,
            "LF2": lf_val_2 if lf_val_2 != -1e6 else np.nan,
            "LF3": lf_val_3 if lf_val_3 != -1e6 else np.nan,
            "sim_duration_sec": loop_duration,
            "status": "success" if lf_val_1 != 1e6 else "failed"
        }
        
        wandb.log({**metrics, "solution_space_3d": fig_3d})
        results_lf_1.append(lf_val_1)
        results_lf_2.append(lf_val_2)
        results_lf_3.append(lf_val_3)
        design_points.append([x1, x2])
        print(f"\n   Logged to WandB: x1={x1:.2f}, x2={x2:.2e}, time duration : {loop_duration:.2f} sec")
        loop_time.append(loop_duration)

end_time = time.time()
print(f"\n--- Grid Search Finished in {end_time - start_time:.2f} sec ---")

wandb.finish()