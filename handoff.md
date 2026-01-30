# PyLCM LEM Branch Handoff

## Branch: `feature/LEM-collision`

## 작업 내용

### 1. LEM (Linear Eddy Model) SGS Mixing 구현

**새로 추가된 파일:**
- `PyLCM/sgs_mixing.py` - LEM SGS mixing 핵심 로직
- `PyLCM/condensation_lem.py` - LEM 적용 응결 루틴
- `LEM_IMPLEMENTATION_NOTES.md` - 구현 문서

**구현 내용:**
- Krueger (1993, JAS) 기반 Linear Eddy Model
- 1D 배열에서 분자/난류 확산 + triplet map rearrangement
- 입자별 T_lem, q_lem → 개별 과포화도(S') 계산
- SAM 연구코드(micro_sgs_mixing.f90) 참조

### 2. 검증 결과

| 지표 | Original | LEM |
|------|----------|-----|
| Peak S | 1.42% | 1.42% |
| S_std | 0% | 0.53% |
| r_std | 1.55 um | 1.86 um |

- S_std 0.53%는 문헌값 범위 내 (0.1-2%)
- 이론적 예측(0.68%)과 유사 (비율 0.78)

### 3. 연구코드와 차이점

| 항목 | 연구코드 | 박스모델 |
|------|----------|----------|
| 난류 소스 | LES tk_LCM | 사용자 eps |
| eta 추적 | 절대값 (kg/kg) | 섭동 S' |
| 격자 | 3D LES | 단일 파셀 |

## 사용법

```python
from PyLCM.condensation_lem import drop_condensation_lem

# 시간 루프 내에서:
# 1. 단열 변화를 T_lem에 적용
for p in particles_list:
    p.T_lem += (T_parcel - T_parcel_old)

# 2. LEM 응결 호출
particles_list, T, q, S_lst, ... = drop_condensation_lem(
    particles_list, T_parcel, q_parcel, P_parcel,
    nt, dt, air_mass, S_lst, rho_aero, True,
    con_ts, act_ts, evp_ts, dea_ts, True,
    diss_rate=1e-4,  # 난류 소산율
    L_domain=100.0   # LEM 도메인 크기
)
```

## 초기화 필요사항

```python
# 입자 초기화시 LEM 변수 추가
p.T_lem = T_parcel + np.random.normal(0, 0.1)  # T fluctuation
p.q_lem = q_parcel + np.random.normal(0, 1e-5)  # q fluctuation
p.lem_id = i  # LEM 배열 인덱스
```

## 주요 파라미터

- `diss_rate`: 난류 소산율 (m²/s³), 기본값 1e-4
- `L_domain`: LEM 도메인 크기 (m), 기본값 100.0
- `tau_nudge`: 넛징 시간스케일 (s), 코드 내 900초

## TODO / 향후 개선

1. [ ] 명시적 eta 추적 (연구코드처럼)
2. [ ] 시간 변화하는 fluctuation source
3. [ ] eps의 물리적 파라미터화
4. [ ] Entrainment 이벤트 시 fluctuation 주입

## 테스트 방법

```bash
# 브랜치에서 테스트 실행
git checkout feature/LEM-collision
python -c "
from PyLCM.sgs_mixing import sgs_mixing_lem
from PyLCM.condensation_lem import drop_condensation_lem
print('LEM modules imported successfully')
"
```

## 관련 참고자료

- Krueger (1993, JAS): LEM 원논문
- SAM LCM: `~/SAM6.10.10.LCM_JS/SRC/MICRO_LAGRANGE/micro_sgs_mixing.f90`
- Grabowski & Abade (2017): 과포화도 변동과 스펙트럼 확장
