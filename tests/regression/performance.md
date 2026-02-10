# Checkpoint benchmark on 2026-01-31
After this checkpoint I will start the refactoring of the state and then see if there is improvement in the time.

```txt
✨ Pixi task (regression in test): python tests/regression/test_output_regression.py                                                                                                                                                           
Running 6 test case(s) with 6 worker(s)...

✓ hupselbrook: regression ok (annual stats match fixture) [2.25s]
✓ surfacewater: regression ok (annual stats match fixture) [3.48s]
✓ grassgrowth: regression ok (annual stats match fixture) [3.79s]
✓ salinitystress: regression ok (annual stats match fixture) [7.65s]
✗ oxygenstress: 
  Mismatches found (tolerance=1e-02):
    Year     Variable       Expected         Actual         Diff
  ------------------------------------------------------------
    1993        MOWDM     12066.3000     12066.9000       0.6000
    1995        MOWDM     13374.9000     13459.9000      85.0000
    1996        MOWDM     13938.7000     13940.7000       2.0000
    1997        MOWDM     12330.0000     12331.0000       1.0000
    1999        MOWDM     14251.0000     14252.1000       1.1000
    2002        MOWDM     12668.3000     12668.4000       0.1000
✓ macropore: regression ok (annual stats match fixture) [170.88s]

============================================================
Results: 5 passed, 1 failed
Total execution time: 170.88s

Individual test timings:
  macropore            170.88s
  oxygenstress         111.47s
  salinitystress         7.65s
  grassgrowth            3.79s
  surfacewater           3.48s
  hupselbrook            2.25s
```