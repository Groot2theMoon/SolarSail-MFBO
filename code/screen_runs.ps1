<#
  screen_runs.ps1 — 포스트버클링 스크리닝: 여러 설정을 100 증분씩 돌리고 지표로 판정한다.

  왜 --inc 인가
    --tstop 은 스텝의 timePeriod 를 줄여 로그의 증분값을 압축된 시간축의 값으로 만든다.
    실측: 증분당 람다가 완주 런의 1/6.2 가 되어 외삽 4.6배 느려 보였다. --inc 는 스텝핑 논리를
    전혀 건드리지 않고 증분 수만 자르므로, 완주 런과 같은 조건에서 비교할 수 있다.

  왜 한 번에 하나인가
    --tstop 과 --speed-discont 를 동시에 걸면 어느 쪽이 원인인지 귀속할 수 없다(실측 교란).
    이 스크립트는 덱 A(기준)와 덱 B(--speed-discont 단독 차이)를 분리해 돌린다.

  돌리는 런 (모두 --inc=100, 완주 런과 같은 라이브러리 사용)
    A_c1  기준              cpus=1
    A_c4  병렬 배수          cpus=4
    A_c8  병렬 배수          cpus=8
    B_c1  speed-discont 효과 cpus=1
    B_c4  조합               cpus=4

  사용법
    .\screen_runs.ps1                 # 덱 생성 + 5개 런 + 판정표
    .\screen_runs.ps1 -SkipDeck       # 이미 만든 덱으로 실행만
    .\screen_runs.ps1 -ParseOnly      # 실행 없이 이미 돌린 잡만 판정

  사전 조건
    - 이 파일을 code\ 폴더에 둔다 (aba\HF_Postbuckle.inp 가 있는 곳)
    - git pull 로 riks_patch_input.py (build 2026-10-01c 이상) 와 screen_report.py 를 받아둔다
#>
param(
  [switch]$SkipDeck,
  [switch]$ParseOnly,
  [int[]]$Cpus = @(1, 4, 8)
)

$ErrorActionPreference = 'Continue'
$CODE = $PSScriptRoot
if (-not $CODE) { $CODE = (Get-Location).Path }
Set-Location $CODE

$LIB = @('--line-search', '--relax-corr')      # 완주 런과 같은 완화 조합
$INC = '--inc=100'                              # 시간축을 건드리지 않는 스크리닝 컷

$DECKS = @{
  A = ($LIB + $INC)
  B = ($LIB + $INC + '--speed-discont')
}

$RUNS = @()
foreach ($c in $Cpus) { $RUNS += @{ n = "A_c$c"; d = 'A'; c = $c } }
$RUNS += @{ n = 'B_c1'; d = 'B'; c = 1 }
$RUNS += @{ n = 'B_c4'; d = 'B'; c = 4 }
$RUNS += @{ n = 'B_c8'; d = 'B'; c = 8 }

if (-not $ParseOnly -and -not $SkipDeck) {
  foreach ($k in @('A', 'B')) {
    $dst = "aba\scr_$k.inp"
    if (Test-Path $dst) { Remove-Item $dst -Force }
    Write-Host "[덱] $dst" -ForegroundColor Cyan
    $f = $DECKS[$k]
    python riks_patch_input.py aba\HF_Postbuckle.inp $dst Step-Postbuckle @f
    if (-not (Test-Path $dst)) {
      Write-Host "[중단] 덱 생성 실패: $dst  (패처 build 확인: 첫 줄에 [패처] build ... 가 찍힌다)" -ForegroundColor Red
      exit 1
    }
  }
}

if (-not $ParseOnly) {
  foreach ($r in $RUNS) {
    $log = "$($r.n).log"
    if (Test-Path $log) { Remove-Item $log -Force }
    Write-Host ("[런] {0}  cpus={1}  시작 {2:HH:mm:ss}" -f $r.n, $r.c, (Get-Date)) -ForegroundColor Cyan
    # ask_delete=OFF: 새 잡 이름만 쓰므로 안전. 2>&1 로 라이선스/스레드 줄까지 로그에 남긴다.
    abaqus job=$($r.n) input="aba\scr_$($r.d).inp" cpus=$($r.c) ask_delete=OFF interactive 2>&1 |
      Tee-Object -FilePath $log | Out-Null
    Write-Host ("      끝   {0:HH:mm:ss}" -f (Get-Date)) -ForegroundColor DarkGray
  }
}

Write-Host ''
python screen_report.py A_c1 A_c4 A_c8 B_c1 B_c4 B_c8
