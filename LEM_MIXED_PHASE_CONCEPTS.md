# LEM for Mixed-Phase Clouds: Conceptual Approaches

## Core Challenges

### Challenge 1: Sedimentation
```
현재 LEM: 입자들이 1D 배열에서 T', q' 섭동을 가짐
문제: 침강하면서 입자들이 다른 속도로 이동 → 배열 구조 깨짐
     큰 입자는 빨리 낙하 → 어떤 S'를 경험하는가?
     빠져나간 입자의 S' slot은 어떻게 채우는가?
```

### Challenge 2: Ice Particle Statistics
```
물방울: ~1000개 super-droplet, 각각 multiplicity ~10^6
빙정:   ~10개 super-particle, 각각 multiplicity ~10^3

LEM triplet map을 10개 입자로 할 수 없음
통계적으로 의미 있는 mixing이 안 됨
```

### Challenge 3: Phase Coupling
```
Wegener-Bergeron-Findeisen:
  - S_liq ≈ 0 (물 포화)
  - S_ice > 0 (빙정은 과포화 상태)
  - 같은 T, q에서 다른 성장률

LEM에서 T', q' fluctuation → S'_liq, S'_ice 둘 다 영향
빙정 성장이 주변 q 소비 → 물방울 증발
이 coupling을 어떻게 표현?
```

---

## Possible Approaches

### Approach A: Phase-Separated LEM

```
┌─────────────────────────────────────────┐
│           Liquid LEM Domain              │
│   [T'₁, q'₁] [T'₂, q'₂] ... [T'ₙ, q'ₙ]  │
│      ↓ diffusion + triplet map           │
│   S'_liq for each liquid droplet         │
└─────────────────────────────────────────┘
                    ↓
          mean(S'_liq), std(S'_liq)
                    ↓
┌─────────────────────────────────────────┐
│           Ice Particles                  │
│   Sample S'_ice from N(0, σ) based on   │
│   liquid phase statistics                │
└─────────────────────────────────────────┘
```

**장점:**
- 물방울 LEM 통계 유지
- 구현 간단

**단점:**
- 빙정-수증기 coupling 무시
- 빙정 S' 독립성 가정

---

### Approach B: Spatial LEM (Fortran 연구코드 방식)

```
물리적 1D 도메인 (z 방향):
    z=0  ────────────────────────  z=L
     │ T(z), q(z) 스칼라 장       │
     │ ← triplet map + diffusion → │
     └─────────────────────────────┘
           ↑ 입자들이 z 위치 가짐
           │ 침강: z → z - v*dt
           │ cyclic BC (modulo)
```

**장점:**
- 물리적으로 일관됨
- 침강이 자연스럽게 포함

**단점:**
- Cyclic BC는 비물리적 (낙하한 입자가 다시 나타남)
- 빙정 수가 적으면 여전히 sparse sampling
- 구현 복잡

---

### Approach C: Lagrangian S' History

```
각 입자가 자신의 S' history를 추적:

입자 i:
  S'(t) = S'(t-dt) * exp(-dt/τ_corr) + noise

τ_corr: 난류 상관 시간 (~Kolmogorov timescale)
noise: 새로운 fluctuation 주입

침강 효과:
  빠른 입자: 빠르게 다른 eddy로 이동 → τ_corr 감소
  느린 입자: 같은 eddy에 오래 머묾 → τ_corr 증가
```

**장점:**
- 침강 자연스럽게 처리
- 입자 수에 무관
- Lagrangian 관점과 일관

**단점:**
- τ_corr 파라미터화 필요
- triplet map의 공간 구조 없음

---

### Approach D: Conditional Statistics

```
아이디어: 입자 크기/종류에 따라 다른 S' 분포

1. LEM에서 S'(z) 계산
2. 침강 속도에 따라 S' "경험" 시간 조정:
   - 빠른 입자: S' 평균으로 relaxation 빠름
   - 느린 입자: 현재 S' 오래 유지

effective S':
  S'_eff = S'_local * f(v_sedi/w_turb)

where f → 0 as v_sedi >> w_turb (fast fall = mean field)
      f → 1 as v_sedi << w_turb (slow fall = local S')
```

**장점:**
- 물리적 직관에 기반
- 연속적 전환

**단점:**
- f 함수 형태 불확실
- 검증 어려움

---

### Approach E: Hybrid (Recommended for Research)

```
상황에 따라 다른 방법 적용:

Phase 1: 활성화/초기 성장 (t < t_collection)
  └── 침강 미미 → Full LEM for liquid
  └── 빙정 없거나 적음 → mean field for ice

Phase 2: 혼합상 발달 (t_collection < t < t_glaciation)
  └── 침강 시작 → Lagrangian S' history
  └── 빙정 증가 → sample from liquid S' distribution

Phase 3: 완전 빙결 (t > t_glaciation)
  └── 물방울 없음 → mean field only
  └── 빙정끼리 aggregation 중심
```

---

## Research Questions

1. **Sedimentation-LEM coupling:**
   - 침강하는 입자가 통과하는 S' field를 어떻게 정의?
   - Open BC vs Cyclic BC의 물리적 의미?

2. **Ice particle statistics:**
   - 최소 몇 개의 빙정이 있어야 LEM이 의미있는가?
   - Variable multiplicity가 triplet map과 충돌하는가?

3. **Phase coupling:**
   - WBF 과정에서 S' fluctuation의 역할?
   - 빙정 성장이 주변 물방울에 미치는 영향의 공간적 범위?

4. **Validation:**
   - DNS나 LES 결과와 비교 가능한가?
   - 관측에서 S' fluctuation 측정 가능한가?

---

## Suggested Starting Point

**Step 1:** Lagrangian S' history 구현 (Approach C)
- 간단하고 침강 문제 자연스럽게 해결
- τ_corr = L_turb / sqrt(TKE) 로 시작

**Step 2:** 물방울에 대해 기존 LEM과 비교
- S' 분포 유사한지 검증

**Step 3:** 빙정에 적용
- 물방울 S' 통계를 기반으로 sampling

**Step 4:** WBF 효과 추가
- 빙정 주변 q 감소 → 물방울 증발 가속
