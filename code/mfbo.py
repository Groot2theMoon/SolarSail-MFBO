"""
Solar-Sail의 형상 변수(x1, x2)를 최적화하여 추력 성능 저하를 최소화하는 스크립트.

스크립트(MFBO)와 구조해석 솔버(Abaqus)는 의도적으로 분리됨
1. MFBO Script (Current File): 
    - PyTorch/BoTorch 기반으로 GPU를 활용하여 Gaussian Process 학습 및 Acquisition Function 계산.
    - Python 3.x 환경.
2. Abaqus Script (run_abaqus.py): 
    - Abaqus 내장 Python 인터프리터(보통 Python 2.7)에서 실행됨.
    - 라이센스 체크 및 FEM 해석 수행.

-> 두 환경의 충돌(Python 버전, GPU 점유 등)을 막기 위해 subprocess를 통해 
    CLI 명령어로 Abaqus를 호출하고, 표준 출력(stdout) 텍스트를 파싱하여 결과값을 통신함.

[최적화 문제 정의: Minimization via Negation]
- 목적(Objective): 주름(Wrinkle)에 의한 추력 감소량 & 압축 영역 최소화.
- BoTorch 특성: 기본적으로 최대화(Maximization) 문제를 풂.
- 해결 전략: Abaqus에서 나온 결과값(양수, Penalty)에 (-1)을 곱해 음수로 변환.
    (예: 실제값 10 -> 최적화 대상 -10). -10이 -100보다 크므로, 
    BoTorch는 0에 가까운(실제 값이 작은) 해를(실제로는 최솟값) 찾게 됨.
"""

import subprocess
import os
import torch 
import warnings
from scipy.stats import qmc
from botorch.models.gp_regression_fidelity import SingleTaskMultiFidelityGP
from botorch.models.transforms.outcome import Standardize
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition.max_value_entropy_search import  qMultiFidelityLowerBoundMaxValueEntropy
from botorch.acquisition.cost_aware import InverseCostWeightedUtility
from botorch.models.cost import AffineFidelityCostModel
from botorch.optim import optimize_acqf_mixed
import wandb

warnings.filterwarnings("ignore")

#wandb 에 대한 부분은 사용자의 실제 id와 사용하고자 하는 프로젝트명 / id에 따라 수정 가능

CHECKPOINT_FILE = "mfbo_checkpoint.pt"
WANDB_PROJECT = "solar-sail-mfbo"
WANDB_RUN_ID_FILE = "wandb_run_id.txt"
FAIL = float("nan")

def save_checkpoint(train_x, train_y, iteration, wandb_run_id):
    """
    현재 상태를 파일로 저장
    시스템 크래시나 중단 시, 이 파일을 통해 이어서 학습가능
    """
    torch.save({
        'train_x': train_x,
        'train_y': train_y,
        'iteration': iteration, # 몇 번째 루프까지 돌았는지
        'wandb_run_id': wandb_run_id
    }, CHECKPOINT_FILE)
    print(f"  [Checkpoint] Saved at iteration {iteration}")

def load_checkpoint():
    """
    체크포인트 파일이 있으면 로드
    Returns:
        tuple: (train_x, train_y, iteration, wandb_run_id) or (None...)
    """
    if os.path.exists(CHECKPOINT_FILE):
        print("  [Checkpoint] Found existing checkpoint. Loading...")
        ckpt = torch.load(CHECKPOINT_FILE)
        return ckpt['train_x'], ckpt['train_y'], ckpt['iteration'], ckpt['wandb_run_id']
    else:
        return None, None, 0, None

# run_abaqus.py 기반으로 작성
def get_abaqus(new_x, new_s):
    """
    new_x : [x1, x2]
    new_s : 0.0(LF) or 1.0(HF)

    Returns:
        tuple or float: 
            - HF 모드일 때: (lf_val, hf_val) 반환. HF 해석 시 LF 결과도 부산물로 얻어짐.
            - LF 모드일 때: (lf_val, None) 반환.
            - 에러 발생 시: 1e6 (Penalty) 또는 -1e6 반환.
            
    Note:
        Abaqus는 표준 출력(stdout)에 "RESULTS:val1,val2..." 형식으로 결과를 찍어야 함.
    """
    mode_str = "HF" if new_s > 0.5 else "LF"
    x1, x2 = new_x[0], new_x[1]
    # M-15: Abaqus 없이 배관을 검증하는 mock oracle (MFBO_MOCK=1 일 때만)
    if os.environ.get("MFBO_MOCK") == "1":
        import numpy as _np
        from eval_currin_mf import evaluate_currin_mf
        _lf = float(_np.ravel(evaluate_currin_mf([[x1, x2]], 0.0))[0])
        _hf = float(_np.ravel(evaluate_currin_mf([[x1, x2]], 1.0))[0])
        print(f"--- [MOCK] {mode_str} currin lf={_lf:.6f} hf={_hf:.6f} ---")
        return (_lf, _hf) if mode_str == "HF" else (_lf, None)

    command = f"abaqus cae noGUI=run_abaqus.py -- {mode_str} {x1} {x2}"
    # [P0-C] HF 는 같은 설계점의 LF 결과와 짝지어야 한다. HF 모델에는 Step-HighTension 이 없어
    #        eval_abaqus.py 의 get_lf(args[-3]) 가 실패한다 -> 먼저 LF 해석을 수행해 둔다.
    #        (stale 산출물 오염 방지를 위해 기존 파일 삭제 후 실행)
    if mode_str == "HF":
        _stale = ("LF_Analysis.odb", "LF_Analysis.lck", "LF_Analysis.msg", "LF_Analysis.sta",
                  "LF_Analysis.dat", "LF_Analysis.fil", "LF_Analysis.prt")
        for _f in _stale:
            if os.path.exists(_f):
                try:
                    os.remove(_f)
                except OSError:
                    pass
        print(f"--- Running Abaqus [LF prerequisite] x1: {x1:.6f} x2: {x2:.6f} ---")
        subprocess.run(f"abaqus cae noGUI=run_abaqus.py -- LF {x1} {x2}",
                       shell=True, check=False, capture_output=True, text=True)
        if not os.path.exists("LF_Analysis.odb"):
            print("!!! HF prerequisite LF run produced no LF_Analysis.odb - aborting this HF point.")
            return FAIL, FAIL

    print(f"--- Running Abaqus [{mode_str}] x1: {x1:.6f} x2: {x2:.6f} ---")
    
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            capture_output=True, 
            text=True  # 문자열로 디코딩
        )
        
        output_lines = result.stdout.splitlines()
        parsed_data = None

        
        for line in output_lines:
            if line.strip().startswith("RESULTS:"):
                # "RESULTS:" 뒷부분 파싱 (예: "RESULTS:0.001,0.02")
                data_str = line.strip().split("RESULTS:")[1].strip()
                if data_str.upper() == "FAIL":        # eval_abaqus.py 의 명시적 실패 신호
                    print("!!! Simulation reported FAIL.")
                    return FAIL, FAIL
                try:
                    parsed_data = [float(val) for val in data_str.split(',')]
                except ValueError:
                    print("!!! Unparsable RESULTS payload: %r" % data_str)
                    return FAIL, FAIL
                break
        
        if parsed_data is None:
            print("!!! Error: Could not find 'RESULTS:' tag in output.")
            print("--- Stderr Log ---")
            print(result.stderr) # 에러 로그 출력
            print("--- Stdout (tail 40) ---")
            print("\n".join(result.stdout.splitlines()[-40:]))
            return FAIL, FAIL

        # 데이터 반환 로직
        if mode_str == "HF":
            # HF는 lf1, lf2, lf3, hf 4개를 반환
            if len(parsed_data) >= 4:
                print(f"  >> [Done] LF: {parsed_data[0]:.6e}, HF: {parsed_data[1]:.6e}")
                return parsed_data[0], parsed_data[3]
            else:
                 print("!!! Error: HF result should have at least 2 values.")
                 return FAIL, FAIL
        else:
            # LF는 lf1, lf2, lf3 총 3개
            print(f"  >> [Done] LF: {parsed_data[0]:.6e}")
            return parsed_data[0], None

    except subprocess.CalledProcessError as e:
        print(f"!!! Abaqus execution failed (Return Code {e.returncode})")
        print(e.stderr)
        return FAIL, FAIL
    except Exception as e:
        print(f"!!! System Error: {e}")
        return FAIL, FAIL

def unnormalize_params(norm_x):
    """
    Unit Hypercube [0, 1] 상의 정규화된 파라미터를 실제 물리적 범위로 변환.
    Params:
        x1: [0, 1] -> [0, 1] (클램프 위치)
        x2: [0, 1] -> [1e-4, 1] (꼭짓점 케이블에 대한 클램프 케이블의 변위비율)
    """
    
    x1_min, x1_max = 0.05, 0.95
    x2_min, x2_max = 1e-4, 1.0
    
    real_x1 = x1_min + norm_x[0] * (x1_max - x1_min)
    real_x2 = x2_min * (x2_max / x2_min) ** norm_x[1] # 지수적 스케일링
    
    return [real_x1, real_x2]

N_ITERATIONS = 5
LF_INIT = 5
HF_INIT = 1

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dtype = torch.double
print(f"Using device: {device}")

train_x, train_y, start_iter, run_id = load_checkpoint()

if run_id is None:
    # 처음 시작하는 경우
    run = wandb.init(project=WANDB_PROJECT, resume="allow")
    run_id = run.id
else:
    print(f"  [WandB] Resuming run {run_id}...")
    run = wandb.init(project=WANDB_PROJECT, id=run_id, resume="must")

if train_x is None:
    print("--- Initializing New Experiment (Independent LHS Strategy) ---")
    
    # 총 LF 개수: LF_INIT
    # HF 수행 개수: HF_INIT
    # 순수 LF 수행 개수: LF_INIT - HF_INIT
    
    n_pure_lf = LF_INIT - HF_INIT
    if n_pure_lf < 0:
        raise ValueError("HF_INIT cannot be larger than LF_INIT")

    # HF용 LHS 샘플링
    sampler_hf = qmc.LatinHypercube(d=2, rng=123)
    samples_hf_norm = torch.tensor(sampler_hf.random(n=HF_INIT), device=device, dtype=dtype)
    
    # 2. 순수 LF용 LHS 샘플링
    sampler_lf = qmc.LatinHypercube(d=2, rng=456) # 다른 시드 사용
    samples_lf_norm = torch.tensor(sampler_lf.random(n=n_pure_lf), device=device, dtype=dtype)
    
    train_x_list = []
    train_y_list = []

    print(f"--- Initializing: {HF_INIT} HF points + {n_pure_lf} LF points ---")

    # HF 샘플 (LF도 같이 저장)
    for i, x_norm in enumerate(samples_hf_norm):
        real_params = unnormalize_params(x_norm.cpu().numpy())
        print(f"  [Init HF {i+1}/{HF_INIT}] Running HF (and LF)...")
        
        lf_val, hf_val = get_abaqus(real_params, 1.0)
        if lf_val != lf_val or hf_val != hf_val:      # NaN = 시뮬레이션 실패 -> 관측 미추가
            print("  !! HF run failed - observation skipped (no data added).")
            continue

        # 최소화를 위해 값을 음수로 변환
        lf_val = -lf_val 
        hf_val = -hf_val
        
        x_lf = torch.cat([x_norm, torch.tensor([0.0], device=device, dtype=dtype)])
        train_x_list.append(x_lf)
        train_y_list.append(lf_val)
        
        x_hf = torch.cat([x_norm, torch.tensor([1.0], device=device, dtype=dtype)])
        train_x_list.append(x_hf)
        train_y_list.append(hf_val)

    # 순수 LF 샘플
    for i, x_norm in enumerate(samples_lf_norm):
        real_params = unnormalize_params(x_norm.cpu().numpy())
        print(f"  [Init LF {i+1}/{n_pure_lf}] Running LF Only...")
        
        lf_val, _ = get_abaqus(real_params, 0.0)
        if lf_val != lf_val:      # NaN = 시뮬레이션 실패 -> 관측 미추가
            print("  !! LF run failed - observation skipped (no data added).")
            continue
        lf_val = -lf_val
        
        x_lf = torch.cat([x_norm, torch.tensor([0.0], device=device, dtype=dtype)])
        train_x_list.append(x_lf)
        train_y_list.append(lf_val)

    # 텐서 변환
    train_x = torch.stack(train_x_list)
    train_y = torch.tensor(train_y_list, device=device, dtype=dtype).view(-1, 1)

    print("--- Initialization Complete ---")
    n_lf_total = (train_x[:, -1] == 0.0).sum().item()
    n_hf_total = (train_x[:, -1] == 1.0).sum().item()
    print(f"Total Data Points: {len(train_y)} (LF: {n_lf_total}, HF: {n_hf_total})")
    
    save_checkpoint(train_x, train_y, 0, run_id)
    start_iter = 0

# LF로 선형해석, HF로 포스트버클링을 수행하므로 두 충실도간 노이즈 차이를 다르게 설정
LF_NOISE = 1e-6
HF_NOISE = 1e-3

# MFBO Loop
# Cost Model 설정: HF 해석이 LF 해석보다 약 10배 비싸다고 가정 (9.0 + 1.0)
cost_model = AffineFidelityCostModel(fidelity_weights={2: 9.0}, fixed_cost=1.0)
cost_utility = InverseCostWeightedUtility(cost_model=cost_model)
new_hf = new_lf = log_hf = log_lf = real_best_val = None
for i in range(N_ITERATIONS):
    try:
        print(f"\nMFBO Iteration {i+1}/{N_ITERATIONS}")

        train_yvar = None   # M-11: 노이즈를 GP 가 학습 (고정하려면 아래 두 줄 주석 해제)
        # train_yvar = torch.full_like(train_y, LF_NOISE)
        # train_yvar[(train_x[:,-1] > 0.5)] = HF_NOISE
        # Surrogate model로 singletask - MF - GP
        model = SingleTaskMultiFidelityGP(
            train_x, 
            train_y,
            train_Yvar=train_yvar,
            data_fidelities=[2],
            outcome_transform=Standardize(m=1)
        )

        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)

        n_candidates = 2000
        x_dim = train_x.shape[-1]

        candidate_set = torch.rand(n_candidates, x_dim, device=device, dtype=dtype)
        candidate_set[..., -1] = 1.0

        #print(f"Candidate Shape {candidate_set.shape}")
        
        # acquisition function 으로 GIBBON 사용 (MES의 variation)
        mf_gibbon = qMultiFidelityLowerBoundMaxValueEntropy(
            model=model,
            candidate_set=candidate_set,
            cost_aware_utility=cost_utility,
        )

        bounds = torch.tensor([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], device=device, dtype=dtype)

        candidates, acq_value = optimize_acqf_mixed(
            acq_function=mf_gibbon.to(device),
            bounds=bounds,
            fixed_features_list=[{2: 0.0}, {2: 1.0}],
            q=1,
            num_restarts=20,
            raw_samples=1024,
        )

        new_x = candidates[0, :-1] 
        
        new_x_np = new_x.cpu().numpy()
        new_s = candidates[0, -1].item()

        real_params = unnormalize_params(new_x_np)

        if new_s >= 0.5: # HF (Trace-Aware Suggestion)
            print(f"-> Running Trace-Aware Simulation (HF runs, LF inferred)")
            
            lf_val, hf_val = get_abaqus(real_params, 1.0)
            if lf_val != lf_val or hf_val != hf_val:      # NaN = 시뮬레이션 실패 -> 관측 미추가
                print("  !! HF run failed - observation skipped (no data added).")
                continue
            
            log_lf, log_hf = lf_val, hf_val 
            lf_val = -lf_val
            hf_val = -hf_val

            # train_x에 두 점 추가
            new_hf = torch.cat([new_x, torch.tensor([1.0], device=device, dtype=dtype)])
            new_lf = torch.cat([new_x, torch.tensor([0.0], device=device, dtype=dtype)])
            train_x = torch.cat([train_x, new_hf.unsqueeze(0), new_lf.unsqueeze(0)], dim=0)
            
            # 반전된 값(hf_val, lf_val) 추가
            new_y = torch.tensor([[hf_val], [lf_val]], device=device, dtype=dtype)
            train_y = torch.cat([train_y, new_y], dim=0)
            
            print(f"   Captured Trace: LF={log_lf:.4f}, HF={log_hf:.4f} (Optim Target: {hf_val:.4f})")

        else: # LF Suggestion
            print(f"-> Running Single Fidelity Simulation (LF only)")
            
            lf_val, _ = get_abaqus(real_params, 0.0)
            if lf_val != lf_val:      # NaN = 시뮬레이션 실패 -> 관측 미추가
                print("  !! LF run failed - observation skipped (no data added).")
                continue
            
            log_lf = lf_val
            lf_val = -lf_val
            
            new_lf = torch.cat([new_x, torch.tensor([0.0], device=device, dtype=dtype)])
            train_x = torch.cat([train_x, new_lf.unsqueeze(0)], dim=0)
            train_y = torch.cat([train_y, torch.tensor([[lf_val]], device=device, dtype=dtype)], dim=0)
            
            print(f"   Captured LF: {log_lf:.4f} (Optim Target: {lf_val:.4f})")

        print(f"-> Suggestion Processed: x= {unnormalize_params(new_x_np.tolist())}, s = {new_s}")

        # 고충실도(HF) 데이터만 필터링
        is_hf_mask = (train_x[:, 2] > 0.5) 

        if is_hf_mask.any():
            # HF에 해당하는 결과값과 입력값 분리
            hf_y_raw = train_y[is_hf_mask]
            hf_x_raw = train_x[is_hf_mask, :2] # x1, x2 값만 추출
            best_opt_val, best_idx = torch.max(hf_y_raw, dim=0)
            real_best_val = -best_opt_val.item() 
            best_x = unnormalize_params(hf_x_raw[best_idx].squeeze())
            
            print("-" * 30)
            print(f"  [Iteration {i+1} Summary]")
            print(f"  - Current Suggested s: {new_s}")
            print(f"  - Best HF Value (Minimized): {real_best_val:.6f}") # [수정] 멘트 변경
            print(f"  - Best Design: x1={best_x[0].item():.4f}, x2={best_x[1].item():.4f}")
            print("-" * 30)
        else:
            print("  -> Status: No High-Fidelity data has been collected yet.")

        
        # 매 반복 체크포인트 저장 (크래시 시 진행 보존)
        save_checkpoint(train_x, train_y, i + 1, run_id)

        # WandB 로깅
        wandb.log({
            "iteration": i + 1,
            "current_x1": new_x[0],
            "current_x2": new_x[1],
            "current_HF": new_hf,
            "current_LF": new_lf, 
            "best_hf_value_minimized": real_best_val,
            "current_fidelity": new_s
        })
    except Exception as e:
        print(f"!!! CRASH at iteration {i}: {e}")
        save_checkpoint(train_x, train_y, i, run_id)
        raise e

print("\n--- Optimization Finished ---")
wandb.finish()
