Simulation choices for lottery model
================

## Overview

The lottery model is based on Ellner et al. (2022), which itself
formalizes the classic lottery-competition idea of Sale (1977) and
Chesson and Warner (1981). At each time step some adults die, vacating
cells in a spatial lattice, and the vacancies are won at random by
propagules in proportion to how many each species has produced. Ellner
et al. (2022) give each species unique responses to a fluctuating
environment and allow for overlapping generations such that many species
can coexist via the storage effect (Warner and Chesson 1985; Chesson
1994).

The simulation code lives in `scripts/simulations/` as four files:

- **`terrestrial_lottery_sim_functions.R`** – the mechanics of one time
  step of the lottery model on an `M x J` lattice: counting neighbors,
  computing fecundity under local competition and the environment,
  building the seed rain each species casts across the lattice, and
  killing/replacing adults.
- **`lottery_param_functions.R`** – generates a biologically-motivated
  set of species parameters and a competition matrix, plus a
  regional-immigration term and post-hoc diagnostics.
- **`lottery_model_sims.R`** – wires the two files above together to run
  a full simulation and plot cover trajectories through time.
- **`subsample_lottery_plots.R`** – re-runs a simulation and subsamples
  the resulting landscape into small, non-overlapping plots, producing a
  long-format data frame (species, plot, and time-step columns) of
  lagged cover (`cover_t` vs. `cover_tm1`) suitable for fitting a
  dynamic occupancy/cover model.

### One time step of the lottery model

Each time step of `lottery_model_sims.R` (coded in
`terrestrial_lottery_sim_functions.R`) applies four operations to the
lattice `X`, an `M x J` matrix of species IDs:

1.  **Count neighbors.** `kernel_count_fast()` (a vectorized version of
    `kernel_count()`) tabulates, for every cell and every species, how
    many conspecific/heterospecific individuals fall within a square
    neighborhood of radius `nbrhood_radius` (Pacala and Silander 1985).
    `count_neighbors()` computes the (edge-corrected) size of that
    neighborhood for each cell.
2.  **Fecundity.** `fecundity_ll()` gives each individual a Poisson
    number of seeds with mean \(\lambda_i(t) / (1 + C_{i,t})\), where
    \(\lambda_i(t) = \lambda_{\text{max},i} \exp(-\sigma_i (v_t - \mu_i)^2)\)
    is maximized when the shared environmental variable `v_t` matches
    the species’ optimum \(\mu_i\), and the competition term \(C_{i,t}\)
    is built from the neighbor counts and the competition matrix `A_mat`
    (Watkinson 1980; Levine and HilleRisLambers 2009).
3.  **Dispersal / seed rain.** `seed_rain_array_fast()` (vectorized
    version of `seed_rain_array()`) spreads each cell’s seeds outward
    over concentric square rings up to `d_max`, using a truncated,
    discretized exponential kernel `P_d_exp()` so that seed density
    declines monotonically with distance (Nathan and Muller-Landau 2000;
    Clark et al. 1999). Seeds landing in a ring are split evenly among
    the in-lattice cells of that ring. The result is an `M x J x S`
    array of seed rain into every cell from every species.
    `add_immigration()` then adds a small, uniform trickle of seeds from
    a regional species pool to every cell (rate `m` relative to the mean
    local seed rain), which prevents global extinction and keeps
    closed-population artifacts from dominating the dynamics (Hubbell
    2001; Mouquet and Loreau 2003).
4.  **Death and replacement.** `die_replace()` kills each individual
    independently with its species’ death probability, then replaces
    each dead individual with a new one drawn with probability
    proportional to the local seed rain from each species in that cell –
    the actual “lottery.”
5.  **Environment update.** Outside these functions,
    `lottery_model_sims.R` advances a single shared AR(1) environmental
    covariate, `env[t+1] = phi * env[t] + rnorm(sd = sigma_env)`,
    representing a positively autocorrelated climate signal.

## Initializing parameters

Parameterizing the simulations takes care in order to generate
‘realistic-looking’ landscapes with many coexisting species. In an
effort to initialize parameters in such a way as to promote
biologically-realistic patterns, such as long-lived, large species that
are clustered in space, short-lived annuals with high dispersal that can
colonize new openings, etc, we have taken the following approaches.

### Perennial - annual axis

Rather than drawing every species parameter independently,
`lottery_param_functions.R` places each of the `nsp` species at a
position \(z \in [0,1]\) on a perennial (\(z=0\)) to annual (\(z=1\))
life-history axis (`draw_life_history()`, evenly spaced or
Beta-distributed), motivated by the fast-slow continuum in plant life
histories (Salguero-Gómez et al. 2016; Franco and Silvertown 1996; Grime
1977; Stearns 1992). `life_history_params()` then derives every other
species-level parameter from this draw \(z_i\):

- **Death probability** is defined as \(P(\delta_{i,t}) = z_i^2\), where
  \(z_i = 1\) gives a strict annual and \(z_i \to 0\) results in a
  long-lived perennial since \(P(\delta_{i,t}) \to 0\).

- We tie **Fecundity** (\(\lambda_\text{max}\)) to the death rate,
  \(\lambda_i \propto \delta_i\) because the lottery invasion
  criterion is \(\lambda_i/\delta_i > \lambda_r/\delta_r\) for an
  invading species \(i\) and a resident species \(r\), holding
  \(\lambda/\delta\) roughly constant across species keeps fitness
  differences small, so that the *stabilizing* mechanisms in the model
  (intra- \> inter-specific competition, local crowding, the storage
  effect, immigration) can maintain diversity (Chesson 2000; Adler et
  al. 2007; Moles et al. 2004). Fecundity is further divided by the
  fraction of seeds that actually leave the parent’s cell
  (\(1 - P(d=0)\)), where \(d\) is dispersal distance in discrete grid
  steps, since seeds that fall back into an occupied cell are wasted,
  and `fitness_sd` adds log-normal noise so species aren’t exactly
  equal.

- **Dispersal rate** decreases log-linearly from a short-dispersing,
  clumped perennial to a far-dispersing annual, representing
  clonal/lateral spread in perennial herbs versus wind/seed dispersal in
  annuals (Thomson et al. 2011; Klimešová and Bello 2009; Bolker and
  Pacala 1999).

- **Environmental sensitivity** (`sigma`) increases with `z`: annuals
  track the environment more closely, while perennials are
  demographically buffered against bad years
  (<span class="nocase">Morris et al.</span> 2008; Pfister 1998).

- **Environmental optima** are spread evenly over the stationary
  distribution of the shared AR(1) environment, independent of `z`, so
  that every life-history type has some species favored in any given
  year (Adler et al. 2006).

- **Intraspecific competition strength** (`alpha`) is drawn around a
  common mean with log-normal noise.

`comp_matrix_gradient()` turns `z` and `alpha` into a full competition
matrix, where the effect of species `j` on species `i`‘s fecundity
combines two multiplicative terms on top of intraspecific competition :
a limiting-similarity term that peaks when `z_i` and `z_j` are close
(MacArthur and Levins 1967), and an asymmetry term that lets perennials
(low `z`) suppress annuals’ fecundity more than the reverse
(<span class="nocase">Kunstler et al.</span> 2016; Seabloom et al. 2003;
Tilman 1994; Levine and Rees 2002). (`comp_matrix()` and
`comp_matrix2()` are older, simpler alternatives that assign a fixed
number of randomly- or structurally-chosen “non-generic” strong
competitors per species rather than basing competition on the
life-history gradient; they are not used in `lottery_model_sims.R`.)

## Running a simulation and checking it

`lottery_model_sims.R` sets lattice size, neighborhood radius, dispersal
truncation, number of time steps/burn-in, immigration rate, and
environmental AR(1) parameters; draws `nsp` species along the
life-history axis; builds their parameters and competition matrix; and
then loops the four-step time-step logic above over a
randomly-initialized lattice, saving every intermediate lattice and
environment value.

After the run, `landscape_diagnostics()` summarizes, per species (after
discarding `burn_in` steps): mean and final percent cover, whether the
species persists to the last time step, and a clumping index (the
fraction of a species’ own 8-cell neighborhood that is conspecific,
divided by its landscape-wide cover – values above 1 indicate spatial
clumping beyond what random placement would produce). The script then
plots each species’ cover trajectory through time, colored by its
life-history position `z`. Because the model is spatially explicit, its
equilibrium behavior isn’t guaranteed to match the non-spatial
lottery/storage-effect prediction (Durrett and Levin 1994), and because
perennials have long expected lifespans, some species can still be
drifting slowly after the chosen burn-in (Hurtt and Pacala 1995) – both
reasons to inspect `landscape_diagnostics()` and the cover-trajectory
plot rather than assume convergence.

## Subsampling into plots

`subsample_lottery_plots.R` treats one simulated landscape as a stand-in
for a real study area that is only ever observed through small, discrete
plots. It runs `lottery_model_sims.R` to get a full time series of
lattices, then:

1.  `sample_plot_corners()` randomly picks `n_plots` non-overlapping
    `plot_size x plot_size` sub-grids (plots) from the lattice.
2.  `plot_cover()` computes each species’ proportional cover within each
    plot at every time step (including explicit zeros).
3.  `lag_cover()` reshapes this into a long data frame keyed by
    `plot_id`, `time`, and `species`, pairing each time point’s cover
    (`cover_t`) with the previous time point’s cover (`cover_tm1`) and
    the environmental value that drove that transition (`env`, matching
    the `t-1 -> t` convention used in
    `scripts/stan/ricker_ss_hurdle_lognormal.stan`).

The resulting `plot_df` is the kind of plot-level, lagged cover data a
dynamic cover/occupancy model would be fit to, and the script’s final
plot (`cover_t` vs. `cover_tm1`, faceted by species, colored by `env`)
is a quick visual check that the simulated dynamics look sensible.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-adler2006climate" class="csl-entry">

Adler, Peter B., Janneke HilleRisLambers, Phaedon C. Kyriakidis,
Qingfeng Guan, and Jonathan M. Levine. 2006. “Climate Variability Has a
Stabilizing Effect on the Coexistence of Prairie Grasses.” *Proceedings
of the National Academy of Sciences* 103 (34): 12793–98.
<https://doi.org/10.1073/pnas.0600599103>.

</div>

<div id="ref-adler2007niche" class="csl-entry">

Adler, Peter B., Janneke HilleRisLambers, and Jonathan M. Levine. 2007.
“A Niche for Neutrality.” *Ecology Letters* 10 (2): 95–104.
<https://doi.org/10.1111/j.1461-0248.2006.00996.x>.

</div>

<div id="ref-bolker1999spatial" class="csl-entry">

Bolker, Benjamin M., and Stephen W. Pacala. 1999. “Spatial Moment
Equations for Plant Competition: Understanding Spatial Strategies and
the Advantages of Short Dispersal.” *The American Naturalist* 153 (6):
575–602. <https://doi.org/10.1086/303199>.

</div>

<div id="ref-chesson1994multispecies" class="csl-entry">

Chesson, P. 1994. “Multispecies Competition in Variable Environments.”
*Theoretical Population Biology* 45 (3): 227–76.
<https://doi.org/10.1006/tpbi.1994.1013>.

</div>

<div id="ref-chesson2000mechanisms" class="csl-entry">

Chesson, Peter. 2000. “Mechanisms of Maintenance of Species Diversity.”
*Annual Review of Ecology and Systematics* 31 (1): 343–66.
<https://doi.org/10.1146/annurev.ecolsys.31.1.343>.

</div>

<div id="ref-chesson1981lottery" class="csl-entry">

Chesson, Peter L., and Robert R. Warner. 1981. “Environmental
Variability Promotes Coexistence in Lottery Competitive Systems.” *The
American Naturalist* 117 (6): 923–43. <https://doi.org/10.1086/283778>.

</div>

<div id="ref-clark1999dispersal" class="csl-entry">

Clark, James S., Miles Silman, Ruth Kern, Eric Macklin, and Janneke
HilleRisLambers. 1999. “Seed Dispersal Near and Far: Patterns Across
Temperate and Tropical Forests.” *Ecology* 80 (5): 1475–94.
[https://doi.org/10.1890/0012-9658(1999)080\[1475:SDNAFP\]2.0.CO;2](https://doi.org/10.1890/0012-9658(1999)080[1475:SDNAFP]2.0.CO;2).

</div>

<div id="ref-durrett1994discrete" class="csl-entry">

Durrett, R., and S. Levin. 1994. “The Importance of Being Discrete (and
Spatial).” *Theoretical Population Biology* 46 (3): 363–94.
<https://doi.org/10.1006/tpbi.1994.1032>.

</div>

<div id="ref-ellner_2022_spmct" class="csl-entry">

Ellner, Stephen P., Robin E. Snyder, Peter B. Adler, and Giles Hooker.
2022. “Toward a ‘Modern Coexistence Theory’ for the Discrete and
Spatial.” *Ecological Monographs* 92 (4): e1548.
<https://doi.org/10.1002/ecm.1548>.

</div>

<div id="ref-franco1996fastslow" class="csl-entry">

Franco, Miguel, and Jonathan Silvertown. 1996. “Life History Variation
in Plants: An Exploration of the Fast-Slow Continuum Hypothesis.”
*Philosophical Transactions of the Royal Society of London. Series B:
Biological Sciences* 351 (1345): 1341–48.
<https://doi.org/10.1098/rstb.1996.0117>.

</div>

<div id="ref-grime1977strategies" class="csl-entry">

Grime, J. P. 1977. “Evidence for the Existence of Three Primary
Strategies in Plants and Its Relevance to Ecological and Evolutionary
Theory.” *The American Naturalist* 111 (982): 1169–94.
<https://doi.org/10.1086/283244>.

</div>

<div id="ref-hubbell2001neutral" class="csl-entry">

Hubbell, Stephen P. 2001. *The Unified Neutral Theory of Biodiversity
and Biogeography*. Vol. 32. Monographs in Population Biology. Princeton
University Press.

</div>

<div id="ref-hurtt1995recruitment" class="csl-entry">

Hurtt, George C., and Stephen W. Pacala. 1995. “The Consequences of
Recruitment Limitation: Reconciling Chance, History and Competitive
Differences Between Plants.” *Journal of Theoretical Biology* 176 (1):
1–12. <https://doi.org/10.1006/jtbi.1995.0170>.

</div>

<div id="ref-klimesova2009clopla" class="csl-entry">

Klimešová, Jitka, and Francesco de Bello. 2009. “CLO-PLA: The Database
of Clonal and Bud Bank Traits of Central European Flora.” *Journal of
Vegetation Science* 20 (3): 511–16.
<https://doi.org/10.1111/j.1654-1103.2009.01050.x>.

</div>

<div id="ref-kunstler2016traits" class="csl-entry">

<span class="nocase">Kunstler, Georges, Daniel Falster, David A. Coomes,
et al.</span> 2016. “Plant Functional Traits Have Globally Consistent
Effects on Competition.” *Nature* 529 (7585): 204–7.
<https://doi.org/10.1038/nature16476>.

</div>

<div id="ref-levine2009niches" class="csl-entry">

Levine, Jonathan M., and Janneke HilleRisLambers. 2009. “The Importance
of Niches for the Maintenance of Species Diversity.” *Nature* 461
(7261): 254–57. <https://doi.org/10.1038/nature08251>.

</div>

<div id="ref-levine2002annual" class="csl-entry">

Levine, Jonathan M., and Mark Rees. 2002. “Coexistence and Relative
Abundance in Annual Plant Assemblages: The Roles of Competition and
Colonization.” *The American Naturalist* 160 (4): 452–67.
<https://doi.org/10.2307/3079234>.

</div>

<div id="ref-macarthur1967limiting" class="csl-entry">

MacArthur, Robert, and Richard Levins. 1967. “The Limiting Similarity,
Convergence, and Divergence of Coexisting Species.” *The American
Naturalist* 101 (921): 377–85. <https://doi.org/10.1086/282505>.

</div>

<div id="ref-moles2004seeds" class="csl-entry">

Moles, Angela T., Daniel S. Falster, Michelle R. Leishman, and Mark
Westoby. 2004. “Small-Seeded Species Produce More Seeds Per Square Metre
of Canopy Per Year, but Not Per Individual Per Lifetime.” *Journal of
Ecology* 92 (3): 384–96.
<https://doi.org/10.1111/j.0022-0477.2004.00880.x>.

</div>

<div id="ref-morris2008longevity" class="csl-entry">

<span class="nocase">Morris, William F., Catherine A. Pfister, Shripad
Tuljapurkar, et al.</span> 2008. “Longevity Can Buffer Plant and Animal
Populations Against Changing Climatic Variability.” *Ecology* 89 (1):
19–25. <https://doi.org/10.1890/07-0774.1>.

</div>

<div id="ref-mouquet2003sourcesink" class="csl-entry">

Mouquet, Nicolas, and Michel Loreau. 2003. “Community Patterns in
Source-Sink Metacommunities.” *The American Naturalist* 162 (5): 544–57.
<https://doi.org/10.1086/378857>.

</div>

<div id="ref-nathan2000dispersal" class="csl-entry">

Nathan, Ran, and Helene C. Muller-Landau. 2000. “Spatial Patterns of
Seed Dispersal, Their Determinants and Consequences for Recruitment.”
*Trends in Ecology & Evolution* 15 (7): 278–85.
<https://doi.org/10.1016/S0169-5347(00)01874-7>.

</div>

<div id="ref-pacala1985neighborhood" class="csl-entry">

Pacala, Stephen W., and J. A. Silander. 1985. “Neighborhood Models of
Plant Population Dynamics. I. Single-Species Models of Annuals.” *The
American Naturalist* 125 (3): 385–411. <https://doi.org/10.1086/284349>.

</div>

<div id="ref-pfister1998variance" class="csl-entry">

Pfister, Catherine A. 1998. “Patterns of Variance in Stage-Structured
Populations: Evolutionary Predictions and Ecological Implications.”
*Proceedings of the National Academy of Sciences* 95 (1): 213–18.
<https://doi.org/10.1073/pnas.95.1.213>.

</div>

<div id="ref-sale1977" class="csl-entry">

Sale, Peter F. 1977. “Maintenance of High Diversity in Coral Reef Fish
Communities.” *The American Naturalist* 111 (978): 337–59.
<https://doi.org/10.1086/283164>.

</div>

<div id="ref-salguerogomez2016fastslow" class="csl-entry">

Salguero-Gómez, Roberto, Owen R. Jones, Eelke Jongejans, et al. 2016.
“Fast–Slow Continuum and Reproductive Strategies Structure Plant
Life-History Variation Worldwide.” *Proceedings of the National Academy
of Sciences* 113 (1): 230–35. <https://doi.org/10.1073/pnas.1506215112>.

</div>

<div id="ref-seabloom2003invasion" class="csl-entry">

Seabloom, Eric W., W. Stanley Harpole, O. J. Reichman, and David Tilman.
2003. “Invasion, Competitive Dominance, and Resource Use by Exotic and
Native California Grassland Species.” *Proceedings of the National
Academy of Sciences* 100 (23): 13384–89.
<https://doi.org/10.1073/pnas.1835728100>.

</div>

<div id="ref-stearns1992evolution" class="csl-entry">

Stearns, Stephen C. 1992. *The Evolution of Life Histories*. Oxford
University Press.

</div>

<div id="ref-thomson2011dispersal" class="csl-entry">

Thomson, Fiona J., Angela T. Moles, Tony D. Auld, and Richard T.
Kingsford. 2011. “Seed Dispersal Distance Is More Strongly Correlated
with Plant Height Than with Seed Mass.” *Journal of Ecology* 99 (6):
1299–307. <https://doi.org/10.1111/j.1365-2745.2011.01867.x>.

</div>

<div id="ref-tilman1994competition" class="csl-entry">

Tilman, David. 1994. “Competition and Biodiversity in Spatially
Structured Habitats.” *Ecology* 75 (1): 2–16.
<https://doi.org/10.2307/1939377>.

</div>

<div id="ref-warner1985storage" class="csl-entry">

Warner, Robert R., and Peter L. Chesson. 1985. “Coexistence Mediated by
Recruitment Fluctuations: A Field Guide to the Storage Effect.” *The
American Naturalist* 125 (6): 769–87. <https://doi.org/10.1086/284379>.

</div>

<div id="ref-watkinson1980density" class="csl-entry">

Watkinson, A. R. 1980. “Density-Dependence in Single-Species Populations
of Plants.” *Journal of Theoretical Biology* 83 (2): 345–57.
<https://doi.org/10.1016/0022-5193(80)90297-0>.

</div>

</div>
