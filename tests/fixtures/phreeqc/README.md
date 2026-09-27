# PHREEQC 공식 예제 (일부)

USGS PHREEQC 3.8.6 배포본(`phreeqc-3.8.6-17100.tar.gz`, https://water.usgs.gov/water-resources/software/PHREEQC/)의
`examples/ex1, ex3, ex4, ex5, ex9` 를 **수정 없이** 복사했다. 사용 조건은 같은 폴더의 `NOTICE`(USGS User Rights Notice).

`tests/test_phreeqc_examples.py` 가 이 입력을 msl A4 경로(phreeqpython IPhreeqc)로 돌려,
공식 phreeqc 실행 파일의 결과(`kb/validation/phreeqc_examples_3.8.6.json`)와 비교한다.
기준값 다시 만들기: `uv run msl bench phreeqc --dist <배포본 폴더> --binary <빌드한 phreeqc> --save kb/validation/phreeqc_examples_3.8.6.json`
