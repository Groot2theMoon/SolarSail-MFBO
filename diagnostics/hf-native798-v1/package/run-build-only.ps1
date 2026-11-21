param([Parameter(Mandatory=$true)][string]$ApprovedPackageSha256)

function Get-PackagePaths {
    @(
        'inputs/baseline-native-memberships.json', 'inputs/boundary-policy-manifest.json',
        'inputs/coating-sides.json', 'inputs/derived-offsets.csv', 'inputs/derived-seeded-coordinates.csv',
        'inputs/hf-build-config.json', 'inputs/hf-q2.json', 'inputs/reference-identity.json',
        'inputs/reference-mesh.inp', 'inputs/verified-transferred-table.txt', 'run-build-only.ps1',
        'staged-code/aba_grid_mesh.py', 'staged-code/aba_imperfection.py', 'staged-code/analysis_config.py',
        'staged-code/build_only_guard.py', 'staged-code/build_planar_hf.py',
        'staged-code/core/__init__.py', 'staged-code/core/constants.py', 'staged-code/ctx.py',
        'staged-code/eval_abaqus.py', 'staged-code/post.py', 'staged-code/riks_patch_input.py',
        'staged-code/run_abaqus.py', 'staged-code/run_abaqus_buckle.py', 'staged-code/solver.py',
        'staged-code/state2_run_evidence.py'
    )
}

function Get-LocalIdentity([string]$Path, [switch]$Exclusive) {
    $share = if ($Exclusive) { [IO.FileShare]::None } else { [IO.FileShare]::Read }
    $stream = [IO.File]::Open($Path, [IO.FileMode]::Open, [IO.FileAccess]::Read, $share)
    try {
        $sha = [Security.Cryptography.SHA256]::Create()
        try {
            @{sha256=([BitConverter]::ToString($sha.ComputeHash($stream)).Replace('-','').ToLowerInvariant()); bytes=$stream.Length}
        } finally { $sha.Dispose() }
    } finally { $stream.Dispose() }
}

function Get-AttemptTree([string]$Root) {
    $result = @{files=@(); directories=@(); issues=@()}
    try {
        $rootItem = Get-Item -LiteralPath $Root -Force -ErrorAction Stop
        if ($rootItem -isnot [IO.DirectoryInfo]) { throw 'Attempt root is not a directory' }
        for ($ancestor = $rootItem; $null -ne $ancestor; $ancestor = $ancestor.Parent) {
            if ($ancestor.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw "Reparse point: $($ancestor.FullName)" }
        }
    } catch { $result.issues += $_.Exception.Message; return $result }
    $queue = New-Object 'System.Collections.Generic.Queue[string]'
    $queue.Enqueue('')
    while ($queue.Count -gt 0) {
        $prefix = $queue.Dequeue()
        try {
            $dir = if ($prefix) { Join-Path $Root $prefix } else { $Root }
            foreach ($item in Get-ChildItem -LiteralPath $dir -Force -ErrorAction Stop) {
                $relative = if ($prefix) { "$prefix/$($item.Name)" } else { $item.Name }
                if ($item.Attributes -band [IO.FileAttributes]::ReparsePoint) {
                    $result.issues += "Reparse point: $relative"
                } elseif ($item -is [IO.DirectoryInfo]) {
                    $result.directories += $relative
                    $queue.Enqueue($relative)
                } elseif ($item -is [IO.FileInfo] -and -not ($item.Attributes -band [IO.FileAttributes]::Device)) {
                    $result.files += $relative
                    foreach ($stream in Get-Item -LiteralPath $item.FullName -Force -Stream '*' -ErrorAction Stop) {
                        if ($stream.Stream -ne ':$DATA') { $result.issues += "Alternate stream: $relative" }
                    }
                } else { $result.issues += "Unexpected item type: $relative" }
            }
        } catch { $result.issues += $_.Exception.Message }
    }
    $result
}

function Assert-ExactPaths([string[]]$Actual, [string[]]$Required, [string[]]$Optional=@()) {
    $legal = [Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
    $seen = [Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
    foreach ($name in ($Required + $Optional)) { [void]$legal.Add($name) }
    foreach ($name in $Actual) {
        if (-not $legal.Contains($name)) { throw "Unexpected path: $name" }
        if (-not $seen.Add($name)) { throw "Missing or duplicate path: $name" }
    }
    foreach ($name in $Required) {
        if (-not $seen.Contains($name)) { throw "Missing or duplicate path: $name" }
    }
}

function Assert-Inventory([string]$Root, [string]$ApprovedDigest) {
    if ($ApprovedDigest -cnotmatch '\A[a-f0-9]{64}\z') { throw 'Malformed approved package digest' }
    $inventory = Join-Path $Root 'build-files.sha256'
    $bytes = [IO.File]::ReadAllBytes($inventory)
    $sha = [Security.Cryptography.SHA256]::Create()
    try { $actual = [BitConverter]::ToString($sha.ComputeHash($bytes)).Replace('-','').ToLowerInvariant() }
    finally { $sha.Dispose() }
    if ($actual -cne $ApprovedDigest) { throw 'Unapproved inventory digest' }
    $ascii = [Text.Encoding]::ASCII.GetString($bytes)
    if (-not $ascii.EndsWith("`n")) { throw 'Inventory is not LF-terminated' }
    $rows = $ascii.Substring(0,$ascii.Length - 1).Split("`n")
    [string[]]$expected = @(Get-PackagePaths)
    [Array]::Sort($expected,[StringComparer]::Ordinal)
    if ($rows.Count -ne $expected.Count) { throw 'Inventory row count differs from the package allowlist' }
    for ($i=0; $i -lt $rows.Count; $i++) {
        if ($rows[$i] -cnotmatch '\A([a-f0-9]{64})  ([a-zA-Z0-9._/-]+)\z' -or $Matches[2] -cne $expected[$i]) {
            throw "Unexpected inventory row at index $i"
        }
        $digest = $Matches[1]
        if ((Get-LocalIdentity (Join-Path $Root $expected[$i])).sha256 -cne $digest) {
            throw "Package file mismatch: $($expected[$i])"
        }
    }
}

function Assert-PreparedPackage([string]$Root, [string]$ApprovedDigest) {
    $tree = Get-AttemptTree $Root
    if ($tree.issues) { throw ($tree.issues -join '; ') }
    Assert-ExactPaths $tree.directories @('inputs','staged-code','staged-code/core')
    Assert-ExactPaths $tree.files (@(Get-PackagePaths) + 'build-files.sha256')
    Assert-Inventory $Root $ApprovedDigest
}

function Assert-JsonFields($Node, [string[]]$Names) {
    if ($Node -isnot [Xml.XmlElement] -or $Node.GetAttribute('type') -cne 'object') { throw 'Expected seal JSON object' }
    Assert-ExactPaths @($Node.ChildNodes | ForEach-Object { $_.LocalName }) $Names
}

function Assert-JsonValue($Node, [string]$Type, [string]$Value) {
    if ($Node -isnot [Xml.XmlElement] -or $Node.GetAttribute('type') -cne $Type -or
        -not [string]::Equals($Node.InnerText,$Value,[StringComparison]::Ordinal)) {
        throw 'Unexpected build seal identity or value type'
    }
}

function Assert-BuildSeal([string]$Root, [string]$ApprovedDigest, [string]$Attempt, [string]$Job) {
    Add-Type -AssemblyName System.Runtime.Serialization -ErrorAction Stop
    $path = Join-Path $Root 'staged-code/aba/build-seal.json'
    $stream = [IO.File]::Open($path,[IO.FileMode]::Open,[IO.FileAccess]::Read,[IO.FileShare]::None)
    try {
        if ($stream.Length -gt 8192) { throw 'Build seal exceeds its record size bound' }
        $memory = [IO.MemoryStream]::new()
        try { $stream.CopyTo($memory); $bytes = $memory.ToArray() } finally { $memory.Dispose() }
    } finally { $stream.Dispose() }
    # The JSON XML mapping retains repeated properties instead of folding them into a PSObject.
    $reader = [Runtime.Serialization.Json.JsonReaderWriterFactory]::CreateJsonReader(
        $bytes, 0, $bytes.Length, [Text.UTF8Encoding]::new($false,$true), [Xml.XmlDictionaryReaderQuotas]::new(), $null)
    try {
        $doc = New-Object Xml.XmlDocument
        $doc.Load($reader)
    } finally { $reader.Dispose() }
    $seal = $doc.DocumentElement
    Assert-JsonFields $seal @('schema_version','attempt','job','package_sha256','native_write_count','raw','sealed')
    Assert-JsonValue ($seal.SelectSingleNode('schema_version')) 'number' '1'
    Assert-JsonValue ($seal.SelectSingleNode('attempt')) 'string' $Attempt
    Assert-JsonValue ($seal.SelectSingleNode('job')) 'string' $Job
    Assert-JsonValue ($seal.SelectSingleNode('package_sha256')) 'string' $ApprovedDigest
    Assert-JsonValue ($seal.SelectSingleNode('native_write_count')) 'number' '1'
    foreach ($kind in 'raw','sealed') {
        $leaf = if ($kind -eq 'raw') { "$Job.raw.inp" } else { "$Job.inp" }
        $relative = "staged-code/aba/$leaf"
        $identity = $seal.SelectSingleNode($kind)
        Assert-JsonFields $identity @('path','sha256','bytes')
        Assert-JsonValue ($identity.SelectSingleNode('path')) 'string' $relative
        try { $actual = Get-LocalIdentity (Join-Path $Root $relative) -Exclusive }
        catch { throw "Cannot read deck exclusively: $relative ($($_.Exception.GetBaseException().Message))" }
        if ($actual.bytes -le 0) { throw "Empty deck: $relative" }
        Assert-JsonValue ($identity.SelectSingleNode('bytes')) 'number' ($actual.bytes.ToString([Globalization.CultureInfo]::InvariantCulture))
        Assert-JsonValue ($identity.SelectSingleNode('sha256')) 'string' $actual.sha256
    }
}

function Get-AttemptProcesses($Snapshot, [string]$Root, [string]$Job, [int]$SelfPid, $Launch=$null) {
    # A command-line token with the exact root or job proves linkage. The generic
    # noGUI argument alone (or an opaque Abaqus process) proves only ambiguity.
    $rootToken = '(?i)(?:^|[\s"=])' + [regex]::Escape($Root.Replace('/','\')) + '(?:$|[\\\s"])'
    $jobToken = '(?i)(?:^|[\s"=])' + [regex]::Escape($Job) + '(?:$|[\s".])'
    $processes = @($Snapshot | Where-Object { $_.ProcessId -ne $SelfPid })
    $linked = @{}
    if ($Launch) { $linked[[int]$Launch.pid] = $Launch }
    do {
        $changed = $false
        foreach ($p in $processes) {
            $id = [int]$p.ProcessId
            if (-not $p.CreationDate -or $linked.ContainsKey($id)) { continue }
            $created = ([datetime]$p.CreationDate).ToUniversalTime()
            $parent = $linked[[int]$p.ParentProcessId]
            $evidence = if (([string]$p.CommandLine).Replace('/','\') -match $rootToken -or $p.CommandLine -match $jobToken) {
                'exact-attempt-token'
            } elseif ($parent -and $parent.created -and $created -ge $parent.created -and
                      ($null -eq $parent.exited -or $created -le $parent.exited)) { 'launch-tree' }
            if ($evidence) {
                $linked[$id] = @{pid=$id; created=$created; exited=$null; linkage=$evidence}
                $changed = $true
            }
        }
    } while ($changed)
    foreach ($p in $processes) {
        $link = $linked[[int]$p.ProcessId]
        $created = if ($p.CreationDate) { ([datetime]$p.CreationDate).ToUniversalTime() }
        $linkage = if ($link -and $created -eq $link.created) {
            $link.linkage
        } elseif (([string]$p.CommandLine).Replace('/','\') -match $rootToken -or $p.CommandLine -match $jobToken -or
                  $p.CommandLine -match '(?i)\bnoGUI=build_planar_hf\.py\b' -or
                  $p.Name -match '(?i)^(?:abaqus|standard|explicit|pre|package|abq\w*|cae\w*)\.exe$') {
            'ambiguous'
        }
        if ($linkage) {
            [ordered]@{pid=[int]$p.ProcessId; parent_pid=[int]$p.ParentProcessId; created_utc=$created
                name=$p.Name; command_line=$p.CommandLine; linkage=$linkage; cleanup='not attempted'}
        }
    }
}

function Get-BuildDisposition([string]$Root, [string]$ApprovedDigest, [string]$Attempt, [string]$Job,
                              [string]$HostName, [string]$LauncherDigest, [string]$Command,
                              $LauncherExit, [bool]$TimedOut, [object[]]$Processes=@(), [string]$RuntimeError='') {
    $tree = Get-AttemptTree $Root
    [string[]]$files = $tree.files
    [Array]::Sort($files,[StringComparer]::Ordinal)
    $record = [ordered]@{
        schema_version=1; attempt=$Attempt; job=$Job; host=$HostName; cwd=(Join-Path $Root 'staged-code')
        package_sha256=$ApprovedDigest; launcher_sha256=$LauncherDigest; command=$Command
        timeout=$TimedOut; launcher_exit=$LauncherExit; processes=@($Processes); writer_count=$Processes.Count
        files=@($files); status='failed'; reason=''
    }
    try {
        if ($tree.issues) { throw ($tree.issues -join '; ') }
        if ($TimedOut) { throw 'Timed out; retain partial outputs and inspect the writer' }
        if ($RuntimeError) { throw $RuntimeError }
        if ($null -eq $LauncherExit -or $LauncherExit -isnot [int] -or $LauncherExit -ne 0) { throw 'Launcher did not return integer exit 0' }
        if ($Processes.Count) { throw 'Attempt-linked or ambiguous process remains; no cleanup authorized' }
        Assert-ExactPaths $tree.directories @('inputs','staged-code','staged-code/core','staged-code/aba','run-records')
        $required = @(Get-PackagePaths) + @(
            'build-files.sha256', "staged-code/aba/$Job.raw.inp", "staged-code/aba/$Job.inp",
            'staged-code/aba/grid_Membrane.inp', 'staged-code/aba/effective-hf.json',
            'staged-code/aba/applied-nodes.tsv', 'staged-code/aba/native-memberships.json',
            'staged-code/aba/attachments.json', 'staged-code/aba/build-seal.json',
            'run-records/prebuild.json', 'run-records/build-only.launch.json', 'run-records/build-only.completion.json',
            'run-records/build-only.stdout.txt', 'run-records/build-only.stderr.txt'
        )
        # postbuild.json must not yet exist; the entry point creates it exclusively.
        Assert-ExactPaths $tree.files $required @(
            'staged-code/abaqus.rpy', 'staged-code/abaqus1.rec', 'staged-code/abaqus_acis.log',
            'staged-code/hf_run_log.txt', 'staged-code/aba/hf_run_log.txt'
        )
        Assert-Inventory $Root $ApprovedDigest
        Assert-BuildSeal $Root $ApprovedDigest $Attempt $Job
        $record.status = 'passed'
        $record.reason = 'Bounded writer evidence only; independent written verification required'
    } catch { $record.reason = $_.Exception.Message }
    $record
}

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
$Root = 'C:\Users\user\AppData\Local\Temp\sail-hf-edge-native798-v1'
$Attempt = 'hf-edge-native798-v1'
$Job = 'SailState2_HF_Native798_Q2_01'
$Prior = 'C:\Users\user\AppData\Local\Temp\sail-state2-q2-OFEMpg0D-r1\run-approved.ps1'
$Launcher = 'C:\Program Files\SIMULIA\Commands\abaqus.bat'
$Stage = Join-Path $Root 'staged-code'
$Log = Join-Path $Root 'run-records'
$LauncherDigest = '8286c0605b53e257e5f90f2fb14132123775d8ad88c3e1f3a2e8eafce1e02202'
$PriorDigest = '5ee55e19a0395b19def100c8ce782322204218cc603bf0c82d1040d70f3d8bdc'
$ArgsText = 'cae noGUI=build_planar_hf.py'
$Command = '"' + $Launcher + '" ' + $ArgsText
if ($env:COMPUTERNAME -ine 'DESKTOP-6PAE7MV' -or $PSScriptRoot -ine $Root) { throw 'Wrong host or attempt root' }
Assert-PreparedPackage $Root $ApprovedPackageSha256
if ((Get-LocalIdentity $Launcher).sha256 -cne $LauncherDigest) { throw 'Unexpected Abaqus launcher identity' }
$tokens = $null; $parseErrors = $null
# Parse the same bytes whose identity is verified, without executing any other lifecycle code.
$priorStream = [IO.File]::OpenRead($Prior)
try {
    if ((Get-LocalIdentity $Prior).sha256 -cne $PriorDigest) { throw 'Unexpected lifecycle identity' }
    $priorReader = [IO.StreamReader]::new($priorStream)
    $ast = [System.Management.Automation.Language.Parser]::ParseInput($priorReader.ReadToEnd(),[ref]$tokens,[ref]$parseErrors)
} finally { $priorStream.Dispose() }
if ($parseErrors) { throw 'Cannot parse pinned lifecycle' }
foreach ($name in 'HashFile','StartAbaqus','RunAbaqus') {
    $defs = @($ast.EndBlock.Statements | Where-Object {
        $_ -is [System.Management.Automation.Language.FunctionDefinitionAst] -and $_.Name -ceq $name
    })
    if ($defs.Count -ne 1) { throw "Missing or duplicate lifecycle helper: $name" }
    . ([scriptblock]::Create($defs[0].Extent.Text))
}
$before = @(Get-CimInstance Win32_Process -OperationTimeoutSec 10 -ErrorAction Stop)
if ($PID -notin $before.ProcessId) { throw 'Cannot establish a complete prelaunch process inspection' }
$active = @(Get-AttemptProcesses $before $Root $Job $PID)
if ($active) { throw 'Attempt-linked or ambiguous build process already active; no launch' }
New-Item -ItemType Directory -Path $Log | Out-Null
[ordered]@{attempt=$Attempt; job=$Job; host=$env:COMPUTERNAME; cwd=$Stage; command=$Command
    package_sha256=$ApprovedPackageSha256; launcher_sha256=$LauncherDigest; lifecycle_sha256=$PriorDigest
    timeout_minutes=30; processes=$active} |
    ConvertTo-Json -Depth 5 | Set-Content -Encoding UTF8 -LiteralPath (Join-Path $Log 'prebuild.json')
$capture = @{process=$null; created=$null; cleanup_suppressed=$false}
$oldDigest = $env:SAIL_APPROVED_PACKAGE_SHA256
$oldBytecode = $env:PYTHONDONTWRITEBYTECODE
$oldPythonPath = $env:PYTHONPATH
$runtimeError = ''
try {
    $env:SAIL_APPROVED_PACKAGE_SHA256 = $ApprovedPackageSha256
    $env:PYTHONDONTWRITEBYTECODE = '1'
    $env:PYTHONPATH = $null
    & {
        function Start-Process([string]$FilePath,[string]$ArgumentList,[string]$WorkingDirectory,
                               [switch]$NoNewWindow,[switch]$PassThru,
                               [string]$RedirectStandardOutput,[string]$RedirectStandardError) {
            $p = Microsoft.PowerShell.Management\Start-Process @PSBoundParameters
            $capture.process = $p
            $null = $p.Handle
            $capture.created = $p.StartTime.ToUniversalTime()
            $p
        }
        # The extracted timeout path asks for PID/tree termination. Decline rather
        # than racing PID reuse or killing an unverified descendant after timeout.
        function taskkill.exe { $capture.cleanup_suppressed = $true }
        RunAbaqus 'build-only' $Stage $ArgsText 30 ''
    }
} catch { $runtimeError = $_.Exception.Message }
finally {
    $env:SAIL_APPROVED_PACKAGE_SHA256 = $oldDigest
    $env:PYTHONDONTWRITEBYTECODE = $oldBytecode
    $env:PYTHONPATH = $oldPythonPath
}
$exitCode = $null; $launch = $null; $timedOut = $capture.cleanup_suppressed
$processes = @()
try {
    if ($capture.process) {
        $launch = @{pid=$capture.process.Id; created=$capture.created; exited=$null; linkage='launch-handle'}
        if ($capture.process.HasExited) {
            $launch.exited = $capture.process.ExitTime.ToUniversalTime()
            $exitCode = $capture.process.ExitCode
        }
    }
    $timedOut = $timedOut -or $runtimeError.StartsWith('Timeout: build-only.',[StringComparison]::Ordinal)
    $after = @(Get-CimInstance Win32_Process -OperationTimeoutSec 10 -ErrorAction Stop)
    if ($PID -notin $after.ProcessId) { throw 'Cannot establish a complete post-exit process inspection' }
    $processes = @(Get-AttemptProcesses $after $Root $Job $PID $launch)
    if ($launch -and -not $launch.exited -and $launch.pid -notin @($processes | ForEach-Object { $_.pid })) {
        $processes += [ordered]@{pid=$launch.pid; parent_pid=$PID; created_utc=$launch.created
            name='cmd.exe'; command_line=$Command; linkage='launch-handle'; cleanup='not attempted'}
    }
} catch { $runtimeError = (@($runtimeError, "Process/exit inspection failed: $($_.Exception.Message)") -join ' ').Trim() }
$result = Get-BuildDisposition $Root $ApprovedPackageSha256 $Attempt $Job $env:COMPUTERNAME $LauncherDigest `
    $Command $exitCode $timedOut $processes $runtimeError
$logItem = Get-Item -LiteralPath $Log -Force
if ($logItem -isnot [IO.DirectoryInfo] -or ($logItem.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
    throw 'Unsafe run-records directory; cannot write postbuild evidence'
}
$post = 'run-records/postbuild.json'
if ($post -cnotin $result.files) { $result.files = @($result.files) + $post }
[Array]::Sort($result.files,[StringComparer]::Ordinal)
$postStream = [IO.File]::Open((Join-Path $Root $post),[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::None)
try {
    $json = [Text.UTF8Encoding]::new($false).GetBytes(($result | ConvertTo-Json -Depth 8) + "`n")
    $postStream.Write($json,0,$json.Length)
} finally { $postStream.Dispose() }
if ($result.status -cne 'passed') { throw $result.reason }
Write-Output $result.reason
