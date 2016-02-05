# tmp_SEIFR
(temporary repo)

Simple Ebola compartmental stochastic model.

S --> E --> I  (infection through a latent stage)

I --> R  (recovery after infection, not dead)

or

I --> F --> B (dead after infection; funerals (F) where transmissions can happen, then removed once buried (B))

Stochastic (Gillespie) simulations use `adaptivetau` R package. ABC fit uses `EasyABC` R package.

Script `seifr_gillespie.R` codes the SEIFR model. The main function to call to simulate is `SEIFR.sim(...)`. See (and run) `example.R` for an example of how to use.

For a back-testing on multiple synthetic data sets, run `Rscript backtest_syntdata.R`

