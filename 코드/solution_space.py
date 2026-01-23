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
            return -1e6, -1e6

        if mode_str == "HF":
            if len(parsed_data) >= 2:
                print(f"  >> [Done] LF: {parsed_data[0]:.6e}, HF: {parsed_data[1]:.6e}")
                return parsed_data[0], parsed_data[1]
            else:
                 print("!!! Error: HF result should have at least 2 values.")
                 return -1e6, -1e6
        else:
            print(f"  >> [Done] LF: {parsed_data[0]:.6e}")
            return parsed_data[0], None

    except subprocess.CalledProcessError as e:
        print(f"!!! Abaqus execution failed (Return Code {e.returncode})")
        print(e.stderr)
        return -1e6, -1e6
    except Exception as e:
        print(f"!!! System Error: {e}")
        return -1e6, -1e6

X1_RANGE = np.linspace(1.5, 2.5, 8)   # Aspect Ratio: 1.5 ~ 2.5 
X2_RANGE = np.linspace(5e-6, 2.0e-5, 8) # Thickness: 5e-6 ~ 2e-5
TARGET_STRAIN = 0.02

wandb.init( 
    entity="dongdong0615-postech",
    # Set the wandb project where this run will be logged.
    project="solarsail-solution-space",
    name=f"grid-search-strain-{TARGET_STRAIN}_{datetime.now().strftime('%Y%m%d_%H%M')}",
    config={
        "target_strain": 0.05,
        "mesh_size": 0.0005,
        "imp_scale_factor": 2.5,
        "x1_range": [1.5, 2.5],
        "x2_range": [7.5e-6, 2.5e-5]
    }
)

print(f"--- Starting Grid Search ---")
print(f"X1 (AR) points: {X1_RANGE}")
print(f"X2 (Thick) points: {X2_RANGE}")
print(f"Total Simulations: {len(X1_RANGE) * len(X2_RANGE)}")

results_lf = []
results_hf = []
design_points = [] # [x1, x2]
loop_time = []

history_x1 = []
history_x2 = []
history_lf = []
history_hf = []

# --- Grid Loop ---
start_time = time.time()

for i, x1 in enumerate(X1_RANGE):
    for j, x2 in enumerate(X2_RANGE):
        loop_start = time.perf_counter() 
        print(f"\n[Grid Step] Processing x1={x1:.4f}, x2={x2:.4e} ({len(design_points)+1}/{len(X1_RANGE)*len(X2_RANGE)})")
        
        lf_val, hf_val = get_abaqus([x1, x2], 1.0)
        
        if lf_val == -1e6 or hf_val == 1e6:
            print("  -> Simulation Failed. Skipping...")
            lf_val, hf_val = np.nan, np.nan # 그래프 그릴 때 제외
        
        loop_duration = time.perf_counter() - loop_start

        history_x1.append(x1)
        history_x2.append(x2)
        history_lf.append(lf_val)
        history_hf.append(hf_val)
        
        fig_3d = go.Figure(data=[go.Scatter3d(
            x=history_x1,
            y=history_x2,
            z=history_hf,
            mode='markers',
            marker=dict(
                size=5,
                color=history_lf,
                colorscale='Viridis',
                colorbar=dict(title="LF"),
                opacity=0.8
            ),
            text=[f"HF: {h:.4f}" for h in history_hf] # 마우스 오버 시 HF 값 표시
        )])

        fig_3d.update_layout(
            title=f"Solution Space (Step {len(history_x1)})",
            scene=dict(
                xaxis_title='Aspect Ratio (x1)',
                yaxis_title='Thickness (x2)',
                zaxis_title='HF'
            ),
            margin=dict(l=0, r=0, b=0, t=40)
        )

        metrics = {
            "x1_AR": x1,
            "x2_Thick": x2,
            "LF": lf_val if lf_val != -1e6 else np.nan,
            "HF": hf_val if hf_val != -1e6 else np.nan,
            "sim_duration_sec": loop_duration,
            "status": "success" if hf_val != -1e6 else "failed"
        }
        
        # WandB에 데이터 전송
        wandb.log({**metrics, "solution_space_3d": fig_3d})
        results_lf.append(lf_val)
        results_hf.append(hf_val)
        design_points.append([x1, x2])
        print(f"\n   Logged to WandB: x1={x1:.2f}, x2={x2:.2e}, time duration : {loop_duration:.2f} sec")
        loop_time.append(loop_duration)

end_time = time.time()
print(f"\n--- Grid Search Finished in {end_time - start_time:.2f} sec ---")

valid_indices = [k for k, v in enumerate(results_hf) if not np.isnan(v)]
clean_lf = np.array([results_lf[k] for k in valid_indices])
clean_hf = np.array([results_hf[k] for k in valid_indices])

valid_indices = [k for k, v in enumerate(results_hf) if not np.isnan(v)]
clean_lf = np.array([results_lf[k] for k in valid_indices])
clean_hf = np.array([results_hf[k] for k in valid_indices])
clean_times = np.array([loop_time[k] for k in valid_indices])

if len(clean_lf) > 2:

    slope, intercept, r_value, p_value, std_err = linregress(clean_lf, clean_hf)
    r_squared = r_value**2
    
    spearman_rho, _ = spearmanr(clean_lf, clean_hf)
    
    best_idx = np.argmin(clean_hf) # HF(Loss)가 가장 낮은 인덱스 (Minimize)
    best_hf_val = clean_hf[best_idx]
    best_x1 = design_points[valid_indices[best_idx]][0]
    best_x2 = design_points[valid_indices[best_idx]][1]
    
    wandb.run.summary["correlation/R2"] = r_squared
    wandb.run.summary["correlation/pearson"] = r_value
    wandb.run.summary["correlation/spearman"] = spearman_rho
    wandb.run.summary["optimum/best_HF_loss"] = best_hf_val
    wandb.run.summary["optimum/best_AR"] = best_x1
    wandb.run.summary["optimum/best_Thickness"] = best_x2
    wandb.run.summary["efficiency/avg_time_sec"] = np.mean(clean_times)
    wandb.run.summary["efficiency/total_time_min"] = (time.time() - start_time) / 60.0
    wandb.run.summary["efficiency/success_rate"] = len(clean_hf) / (len(X1_RANGE) * len(X2_RANGE))

    plt.figure(figsize=(6, 6))
    plt.scatter(clean_lf, clean_hf, c='blue', alpha=0.6)
    plt.plot(clean_lf, slope * clean_lf + intercept, 'r--', label=f'R2={r_squared:.2f}')
    plt.xlabel('LF (Energy Density)')
    plt.ylabel('HF (Thrust Loss Log)')
    plt.title(f'Correlation Analysis (Strain={TARGET_STRAIN})')
    plt.legend()
    plt.grid(True)
    
    wandb.log({"plots/correlation_scatter": wandb.Image(plt)})
    plt.close()

    table = wandb.Table(columns=["x1", "x2", "LF", "HF", "Time"])
    for k in valid_indices:
        table.add_data(
            design_points[k][0], 
            design_points[k][1], 
            results_lf[k], 
            results_hf[k], 
            loop_time[k]
        )
    wandb.log({"raw_data_table": table})

    print(f"--- Final Stats Logged to WandB ---")
    print(f"  R2: {r_squared:.4f}, Best HF: {best_hf_val:.4f}")

else:
    print("Not enough data for statistics.")
    wandb.run.summary["status"] = "failed_insufficient_data"

wandb.finish()