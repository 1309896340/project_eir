@echo off
setlocal

rem ================================================================
rem  GENCODE human release_50 downloader (resume-capable)
rem
rem  NOTE: keep this file pure ASCII. The previous version used
rem  UTF-8 text + "chcp 65001" and broke under the default Chinese
rem  Windows console codepage (GBK). Do not add non-ASCII text.
rem
rem  Usage:
rem    scripts\download_gencode.cmd          core files only (default)
rem    scripts\download_gencode.cmd all      also optional files
rem
rem  Core files (mRNA design / codon optimization):
rem    gencode.v50.pc_transcripts.fa.gz      protein-coding mRNA incl. UTR
rem    gencode.v50.pc_translations.fa.gz     protein sequences (for QC)
rem    gencode.v50.annotation.gtf.gz         annotation (UTR/CDS boundaries)
rem    MD5SUMS                               checksums
rem
rem  "all" mode adds:
rem    GRCh38.primary_assembly.genome.fa.gz  genome FASTA (for gffread)
rem    gencode.v50.metadata.RefSeq.gz        GENCODE-RefSeq ID mapping
rem    gencode.v50.lncRNA_transcripts.fa.gz  lncRNA negative samples
rem    gencode.v50.polyAs.gtf.gz             polyA sites
rem
rem  Downloader: wget -c if found on PATH, otherwise curl -C -
rem  (curl.exe ships with Windows 10+). Both resume partial files and
rem  skip already-complete files, so re-running is always safe.
rem
rem  Do NOT gunzip before all downloads finish: unzipping deletes the
rem  .gz file, so a later re-run would re-download it from scratch.
rem
rem  Output dir: data\gencode\human
rem ================================================================

set "BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50"
set "OUTDIR=%~dp0..\data\gencode\human"
set "MODE=%~1"
if not defined MODE set "MODE=core"
set "FAIL=0"

set "DL="
where wget >nul 2>nul && set "DL=wget"
if not defined DL (
    where curl >nul 2>nul && set "DL=curl"
)
if not defined DL (
    echo [ERROR] Neither wget nor curl found on PATH.
    exit /b 1
)

if not exist "%OUTDIR%" mkdir "%OUTDIR%"
echo Downloader: %DL%
echo Output dir: %OUTDIR%
echo Source:     %BASE%
echo Mode:       %MODE%
echo.

rem ---------------- core files (always) ----------------
call :download gencode.v50.pc_transcripts.fa.gz
call :download gencode.v50.pc_translations.fa.gz
call :download gencode.v50.annotation.gtf.gz
call :download MD5SUMS

rem ---------------- optional files ("all" mode) ----------------
if /i "%MODE%"=="all" (
    call :download GRCh38.primary_assembly.genome.fa.gz
    call :download gencode.v50.metadata.RefSeq.gz
    call :download gencode.v50.lncRNA_transcripts.fa.gz
    call :download gencode.v50.polyAs.gtf.gz
)

echo.
if "%FAIL%"=="0" (
    echo All downloads finished. Verify checksums in Git Bash with:
    echo   cd data/gencode/human ^&^& md5sum -c MD5SUMS --ignore-missing
) else (
    echo %FAIL% file(s) not finished. Re-run this script to resume.
)
endlocal
exit /b %FAIL%

:download
rem subroutine: fetch one file with resume; keep going on failure
echo === Downloading %~1 ===
if "%DL%"=="wget" (
    wget -c --tries=5 --timeout=60 -P "%OUTDIR%" "%BASE%/%~1"
) else (
    curl --fail --location --retry 5 --connect-timeout 60 -C - -o "%OUTDIR%\%~1" "%BASE%/%~1"
)
if errorlevel 1 (
    echo [WARN] %~1 not finished. Partial file kept, re-run to resume.
    set /a FAIL+=1
) else (
    echo [DONE] %~1
)
echo.
goto :eof
