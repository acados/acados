## Differential Drive Robot
Run tests:
```python ./diff_drive/main.py```

Timings of solving one OCP with fast zoRO (BLASFEO) when $\mathrm{max\_iter} = 2$:
![plot](./figures/timings_diff_drive_blasfeo.png)

Timings of solving one OCP with zoRO (numpy) when $\mathrm{max\_iter} = 2$:
![plot](./figures/timings_diff_drive_numpy.png)

## Hanging Chain

Run tests:
```python ./zoro-NMPC-2021/main.py```

CPU time per OCP when $n_\mathrm{mass} = 5$ and $\mathrm{max\_iter} = 50$:
![plot](./figures/timings_nm5.png)

Closest distance to the wall:
![plot](./figures/constraint_violation_nmass_5_seeds_20.png)


## Continuous Stirred-tank Reactor (CSTR)