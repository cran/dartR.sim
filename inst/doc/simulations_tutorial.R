## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 6)
# Files written by the examples go to a temporary directory
knitr::opts_knit$set(root.dir = tempdir())
dartR.base::gl.set.verbosity(0)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("dartRverse")
# library(dartRverse)
# dartRverse_install("dartR.sim", rep = "CRAN")

## -----------------------------------------------------------------------------
library(dartR.sim)

ref_table <- gl.sim.WF.table(
  file_var = system.file("extdata", "ref_variables.csv", package = "dartR.sim"),
  interactive_vars = FALSE)

res_sim <- gl.sim.WF.run(
  file_var = system.file("extdata", "sim_variables.csv", package = "dartR.sim"),
  ref_table = ref_table,
  interactive_vars = FALSE)

res_sim$iteration_1$generation_10   # a genlight object

## -----------------------------------------------------------------------------
file.copy(system.file("extdata", c("ref_variables.csv", "sim_variables.csv"),
                      package = "dartR.sim"), getwd())

## -----------------------------------------------------------------------------
res_sim <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                         interactive_vars = FALSE, phase1 = TRUE,
                         number_pops_phase1 = 2, population_size_phase1 = "200 200",
                         number_pops_phase2 = 2, population_size_phase2 = "50 50")

## -----------------------------------------------------------------------------
res_sim <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                         interactive_vars = FALSE,
                         number_pops_phase2 = 2, population_size_phase2 = "200 200",
                         ne_phase2 = "50 100", dispersal_phase2 = FALSE,
                         gen_number_phase2 = 30, every_gen = 5)
res_sim$iteration_1$generation_30@other$sim.vars$ne_expected   # "50 100"

## -----------------------------------------------------------------------------
res_sim <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                         interactive_vars = FALSE,
                         number_pops_phase2 = 3, population_size_phase2 = "50 50 50",
                         dispersal_phase2 = TRUE, number_transfers_phase2 = 2,
                         transfer_each_gen_phase2 = 2)

## -----------------------------------------------------------------------------
gl.sim.create_dispersal(number_pops = 4, dispersal_type = "all_connected",
                        number_transfers = 1, transfer_each_gen = 1,
                        outpath = getwd())
# edit dispersal_table.csv, then:
res_sim <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                         interactive_vars = FALSE, number_pops_phase2 = 4,
                         population_size_phase2 = "50 50 50 50",
                         dispersal_phase2 = TRUE,
                         file_dispersal = file.path(getwd(), "dispersal_table.csv"))

## -----------------------------------------------------------------------------
ref_table_sel <- gl.sim.WF.table(file_var = "ref_variables.csv", interactive_vars = FALSE,
                                 loci_deleterious = 50)
res_sim <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table_sel,
                         interactive_vars = FALSE, selection_phase2 = TRUE,
                         number_offspring_phase2 = 10)

## -----------------------------------------------------------------------------
res_sim <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                         interactive_vars = FALSE, sib_mating_phase2 = 0.3,
                         population_size_phase2 = "100", number_offspring_phase2 = 20,
                         gen_number_phase2 = 20)
# expected F = 0.3 / (4 - 0.9) = 0.097

## -----------------------------------------------------------------------------
library(dartR.sim)
x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
x <- gl.filter.monomorphs(x, verbose = 0)
x@chromosome <- as.factor(x@other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
x@position <- as.integer(x@other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1)

## -----------------------------------------------------------------------------
ref_x <- gl.sim.WF.table(file_var = "ref_variables.csv", x = x, interactive_vars = FALSE,
                         real_freq = TRUE, chunk_neutral_loci = 0)
res_x <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_x, x = x,
                       interactive_vars = FALSE,
                       real_pops = TRUE, real_freq = TRUE,
                       real_freq_shrink = "auto",          # founders' FST = x's FST
                       real_migration = TRUE,              # keeps it over time
                       dispersal_phase2 = TRUE,
                       population_size_phase2 = "200 200 200",
                       ne_phase2 = "60 60 60",             # e.g. from gl.LDNe()
                       real_sample_size = TRUE,            # stores 23, 17 and 41
                       real_inbreeding = TRUE, sib_mating_phase2 = 0,
                       store_founders = TRUE,
                       gen_number_phase2 = 10, every_gen = 5)

## -----------------------------------------------------------------------------
res_fam <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                         interactive_vars = FALSE, population_size_phase2 = "100",
                         number_offspring_phase2 = 20, sample_percent = 40,
                         sample_families = "9 15 9", sample_parents = TRUE,
                         store_founders = TRUE, store_pedigree = TRUE,
                         gen_number_phase2 = 2, every_gen = 1)
table(res_fam$iteration_1$generation_1@other$ind.metrics$sample_role)

## -----------------------------------------------------------------------------
res <- gl.sim.WF.run(file_var = "sim_variables.csv", ref_table = ref_table,
                     interactive_vars = FALSE, number_pops_phase2 = 2,
                     population_size_phase2 = "100 100", dispersal_phase2 = TRUE,
                     gen_number_phase2 = 50, every_gen = 5)
diag <- gl.diagnostics.sim(res, Ne = 100)
tail(diag$fst)

