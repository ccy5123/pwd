#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
apply_fixes.py — bsa_docking_desolv_pipeline_v2.py 에 FIX①④ 안전 적용.
- anchor 문자열 유일성 assert(count==1) → 오적용 차단
- GSOLV_SCALE 존재 여부로 멱등 가드(이미 적용 시 SKIP)
- 치환 후 ast.parse 로 문법 검증 후에만 파일 기록
- FIX②③은 box·coords 유실(재실행 필요)이라 여기서 자동적용 안 함(pipeline_fixes.md 참조)
실행: (작업폴더에서) python apply_fixes.py
"""
import ast, sys
from pathlib import Path

F = Path("bsa_docking_desolv_pipeline_v2.py")
if not F.exists():
    print(f"[STOP] {F} 없음. 작업폴더 확인."); sys.exit(1)
src = F.read_text()

if "GSOLV_SCALE" in src:
    print("[SKIP] 'GSOLV_SCALE' 이미 존재 → FIX①④ 이미 적용됨. 아무것도 안 함."); sys.exit(0)


def repl(s, old, new, tag):
    n = s.count(old)
    assert n == 1, f"[STOP] [{tag}] anchor count={n} (기대 1). 파일이 이전과 다름 → 중단."
    return s.replace(old, new)


# --- FIX④: GSOLV_SCALE 상수 추가 (L67 뒤) ---
src = repl(
    src,
    "LOG_M_BSA = math.log10(M_BSA)  # = 1.823",
    "LOG_M_BSA = math.log10(M_BSA)  # = 1.823\n"
    "# [FIX4] xTB Gsolv 가산 스케일. 1.0=기존. Vina 스코어는 수용액 결합에 경험적 학습됨 ->\n"
    "#        xTB Gsolv(|mean|~4.4 kcal/mol) 1:1 가산은 수상 desolvation 이중계상 의심\n"
    "#        (reference-state 불일치). knob은 desolvation 감도분석용(재해석은 reinterpret.py).\n"
    "GSOLV_SCALE = 1.0",
    "FIX4 const",
)

# --- FIX①④: 변환 블록 (L1009-1010) ---
src = repl(
    src,
    "                dG_corr = dG_contact - gsolv\n"
    "                logK_pred = -dG_corr / (2.302585 * RT) - LOG_M_BSA",
    "                # [FIX1] LOG_M_BSA(-1.823)=몰->질량기준 변환. albumin_exp가 질량기준\n"
    "                #        log K_BSA/w (L/kg)일 때만 정합(Endo&Goss2011/UFZ 정의 확인 필요).\n"
    "                #        offset은 순수 가산상수 -> pearson/calR2 불변, absR2/bias만 이동.\n"
    "                # [FIX4] GSOLV_SCALE=1.0 기본. 이중계상 caveat는 상단 상수 주석 참조.\n"
    "                dG_corr = dG_contact - GSOLV_SCALE * gsolv\n"
    "                logK_pred = -dG_corr / (2.302585 * RT) - LOG_M_BSA",
    "FIX1/4 block",
)

ast.parse(src)  # SyntaxError 나면 여기서 중단, 파일 미기록
F.write_text(src)
print("[OK] FIX1/4 적용 완료 + ast.parse 통과.")
print("     검증: grep -c 'GSOLV_SCALE' bsa_docking_desolv_pipeline_v2.py  → 3 (정의 L71 + 주석 + 사용 L1017)")
