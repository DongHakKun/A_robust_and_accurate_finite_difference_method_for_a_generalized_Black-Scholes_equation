# Generalized Black–Scholes equation · Unofficial C++ implementation

**비공식 구현체입니다.** Cen–Le (2011)의 논문을 재현하기 위해 작성한 C++17 수치실험 노트북이며, 논문 저자가 배포한 공식 코드가 아닙니다.

Unofficial C++17 notebook implementation reproducing the numerical experiments in Cen–Le (2011), *A robust and accurate finite difference method for a generalized Black–Scholes equation*.

Cen–Le (2011), *A robust and accurate finite difference method for a generalized Black–Scholes equation*의 §5에 제시된 Test 1과 Test 2를 각각 재현한다. 각 노트북에서 문제의 설정, Figure·Table 계산과 결과 해석을 연결한다. 수치방법과 구현에 대한 설명은 별도 참고 노트북에 모았고, 정리는 실험에 필요한 결과만 사용한다.

## 노트북

| 파일 | 내용 |
|---|---|
| [5_1_Test_1_Figure_1_Table_1.ipynb](ipynb/5_1_Test_1_Figure_1_Table_1.ipynb) | Test 1: 시공간 의존 변동성, Figure 1·Table 1 재현 |
| [5_2_Test_2_Figure_2_Table_2.ipynb](ipynb/5_2_Test_2_Figure_2_Table_2.ipynb) | Test 2: 공간 의존 변동성, Figure 2·Table 2 재현 |
| [Numerical_method_and_experiment_setup.ipynb](ipynb/Numerical_method_and_experiment_setup.ipynb) | 수치방법 참고: 표기, payoff smoothing, 격자, 차분식, Thomas 풀이, 오차 정의 |

세 노트북 모두 실제 C++ 실행 출력과 그림을 저장했다. 각 노트북은 독립적으로 위에서부터 실행할 수 있다. 수식에는 `x_i`, `t_j`, `h_i`, `U_i^j`, `A_i^j`, `B_i^j`, `C_i^j`, `L^{N,M}` 표기를 사용한다. 시간은 만기 `j=M`에서 `j=0`으로 진행한다.

## 실행

WSL/Linux에서 Conda를 준비한 뒤 프로젝트 루트에서 환경을 만들고 JupyterLab을 연다.

```bash
conda env create -f environment.yml
conda activate black-scholes-cpp
jupyter lab ipynb/
```

노트북에서 설치된 xeus-cpp의 **C++17** 커널을 선택하고 **Restart Kernel → Run All**을 실행한다. 저장된 커널 이름은 `black-scholes-cpp`이므로, 새 환경에서는 커널을 다시 선택해야 할 수 있다. VS Code를 쓰는 경우 WSL 창에서 연다. 커널의 작업 디렉터리는 `ipynb/`이다.

설정은 [environment.yml](environment.yml)에 있다. 기존 Python 환경과 분리하여 xeus-cpp, Jupyter, gnuplot을 설치한다. [xeus-cpp 공식 설치 안내](https://xeus-cpp.readthedocs.io/en/latest/InstallationAndUsage.html)를 참고할 수 있다.

수치해, 보간, 오차, 표는 C++로 계산하고 그림은 C++에서 gnuplot을 호출해 생성한다. 셀을 수정한 뒤에는 커널을 재시작하고 전체 실행하면 된다.

## 논문 재현 설정과 결과

공통 설정은 `T=1`, `r=0.06`, `K=25`, `Smax=100`, `epsilon=1e-4`이다.

- Test 1: `sigma(x,t)=0.2*(1+t*exp(-x))`, `alpha=0.04`.
- Test 2: `sigma(x,t)=0.4*(2+sin(x))`, `alpha=0.16`.
- 두 실험 모두 `betaStar=0.06`. `alpha`, `betaStar`는 주어진 계수의 최선의 하한·상한으로 선택했고, 이 선택으로 논문 표가 재현된다.
- Figure 1·2: `N=M=64`.
- Table 1·2: `M=1024`, `N=128,256,512,1024`; 참조해 `N=2048, M=1024`.

참조해를 같은 시간층에서 공간 선형 보간하고 **모든 공간·시간 격자점의 최대오차**를 계산한다. 현재 시점 `t=0`의 오차만 구하면 논문 표와 달라진다.

| N | Test 1 계산 오차 | 관측 수렴률 | Test 2 계산 오차 | 관측 수렴률 |
|---:|---:|---:|---:|---:|
| 128 | 6.80267592e-2 | 1.122 | 5.71239536e-2 | 1.551 |
| 256 | 3.12581899e-2 | 1.312 | 1.94981844e-2 | 1.919 |
| 512 | 1.25878300e-2 | 1.880 | 5.15498502e-3 | 2.284 |
| 1024 | 3.41885173e-3 | — | 1.05863181e-3 | — |

여덟 오차값과 여섯 수렴률 모두 원문 Table 1·2의 인쇄 정밀도와 일치한다. 수렴률은 현재 행과 다음 행의 오차 비율로 계산한다. 참조해도 수치해이므로 이 값은 정확해에 대한 총오차나 독립적인 시간 수렴률을 의미하지 않는다. 비균일 격자 접합점의 국소 절단오차에 관한 주의사항도 노트북에 설명했다.

## 코드와 결과 파일

핵심 계산은 `thomas_class`의 간단한 클래스, 기본 배열, `new[]/delete[]`, 명시적 반복문과 `printf` 스타일을 따른다.

| 파일 | 역할 |
|---|---|
| [include/TridiagonalSolver.h](include/TridiagonalSolver.h) | Thomas algorithm, 작업 배열 재사용 |
| [include/BlackScholesSolver.h](include/BlackScholesSolver.h) | 격자·payoff·계수·경계조건·시간 적분 |
| [include/ErrorAnalysis.h](include/ErrorAnalysis.h) | 선형 보간, 시공간 최대오차, 상수 계수 정확해 검증 |
| [include/ExperimentOutput.h](include/ExperimentOutput.h) | CSV와 gnuplot 그림 출력 |
| [include/NotebookDisplay.h](include/NotebookDisplay.h) | xeus-cpp 셀 안의 SVG 표시 |

노트북을 실행하면 `results/`를 자동으로 만들고 `test1_errors.csv`, `test2_errors.csv`, 시간 단면 CSV, 표면도·단면도·수렴 그림의 SVG 및 PNG를 저장한다. `common_mesh.*`, `common_payoff.*`는 공통 설명용 그림이다. `.dat`, `.gp`는 그림 생성을 위한 재생성 가능한 중간 파일이다. `results/`는 Git 추적에서 제외하며, 저장된 실행 출력과 그림은 노트북 안에서 확인할 수 있다.

## 검증

Thomas 풀이, 격자와 payoff, 경계조건·시간 방향·잔차, 내부 M-matrix, 비중첩 격자 보간, 전체 시간 최대오차, 논문 표 재현, 상수 변동성 Black–Scholes 정확해에 대한 오차 감소를 확인했다. 전체 테스트와 AddressSanitizer·UBSan 검사를 통과했다. xeus-cpp에서도 세 노트북의 모든 코드 셀과 그림 출력을 검증했다.

## 참고 논문

Cen–Le (2011), *A robust and accurate finite difference method for a generalized Black–Scholes equation*. [DOI: 10.1016/j.cam.2011.01.018](https://doi.org/10.1016/j.cam.2011.01.018).

논문 PDF, 개인 참고 자료(`docs/`), 작업 이력(`.archive/`)은 이 저장소에 포함하지 않는다.
