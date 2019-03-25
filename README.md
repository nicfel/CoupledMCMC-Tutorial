---
author: Nicola F. Müller
level: Intermediate
title: AIM in StarBeast2 v0.15.2 Tutorial
subtitle: Setting up an Bayesian Skyline with coupled MCMC
beastversion: 2.5.2
---


# Background

Coupled MCMC (also called parallel tempering or MC^3) is a Bayesian approach that uses heated chains in order to traverse unfavourable intermediate states more easily and in order to parallelise analyses.
Coupled MCMC works by having 1 cold chain which works exactly the same as a standard MCMC chain and one or more heated chains.
The chains are heated in temperature increments defined in this implementations as *deltaTemperature*.
These heated chains have increase acceptance probabilities, making it more easily for them to explore parts of the state space that are less likely.
In turn however, these heated chains to not explore the posterior probability space such that the frequency at which they visit a state is proportional to their probability.
The can however be used to propose new states.
In order to do so, coupled MCMC proposes to exchange the states of two random chains after each chain is run for some amount of iterations in what amounts to a MCMC move.
After how many iterations the states of two random chains are proposed to be swapped is in this implementation called *resampleEvery*.

The webpage [https://darrenjw.wordpress.com/2013/09/29/parallel-tempering-and-metropolis-coupled-mcmc/](https://darrenjw.wordpress.com/2013/09/29/parallel-tempering-and-metropolis-coupled-mcmc/), gives a good introduction into how coupled MCMC works.


----

# Programs used in this Exercise

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 ([http://www.beast2.org](http://www.beast2.org)) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees. This tutorial is written for BEAST v{{ page.beastversion }} {% cite BEAST2book2014 --file AIM-Tutorial/master-refs %}.


### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### TreeAnnotator

TreeAnnotator is used to summarise the posterior sample of trees to produce a maximum clade credibility tree. It can also be used to summarise and visualise the posterior estimates of other tree parameters (e.g. node height).

TreeAnnotator is provided as a part of the BEAST2 package so you do not need to install it separately.


----

# Practical: Setting up an analysis with coupled MCMC

{% cite Mueller348391 --file AIM-Tutorial/master-refs.bib %}



In this tutorial, we will describe the two different ways to setup a BEAST2 analysis to run with coupled MCMC.
To do so, we will setup a Bayesian Skyline plot analysis by following analogue to the tutorial on [skyline plots](https://taming-the-beast.org/tutorials/Skyline-plots/).

Setting up an analysis in BEAUTi is currently only possible for analyses that only use the standard template, i.e. such for which setting up an analysis does not require to load a template.
All other analyses have to be setup by editing one line in the `*xml` file.


## The Data
The dataset consists of an alignment of 63 Hepatitis C sequences sampled in 1993 in Egypt {% cite Ray2000 --file Skyline-plots/master-refs %}. This dataset has been used previously to test the performance of skyline methods {% cite Pybus2003, Drummond2005, Stadler2013 --file Skyline-plots/master-refs %}.

With an estimated 15-25%, Egypt has the highest Hepatits C prevalence in the world. In the mid 20^(th) century, the prevalence of Hepatitis C increased drastically (see [Figure 1](#fig:prevalence) for estimates). We will try to infer this increase from sequence data.

<figure>
	<a id="fig:prevalence"></a>
	<img style="width:50%;" src="figures/Estimated_number_hcv.png" alt="">
	<figcaption>Figure 1: The estimated number of Hepatitis C cases in Egypt {% cite Pybus2003 --file Skyline-plots/master-refs.bib %}.</figcaption>
</figure>
<br>


### Download CoupledMCMC
First, we have to download the CoupledMCMC package by using the BEAUTi package manager. Go to `File >> Manage Packages` and download the package and CoupledMCMC StarBeast2.

<figure>
	<a id="fig:example1"></a>
	<img style="width:50%;" src="figures/Package.png" alt="">
	<figcaption>Figure 1: Download the CoupledMCMC package</figcaption>
</figure>

## For analyses that don't require to load a template

If for your analysis, no loading of a template is required, we can setup the analysis in BEAUTi:

### Loading the template

Next, we have to load the BEAUTi template from `File`, select `Template >> CoupledMCMC`.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/Template.png" alt="">
	<figcaption>Figure 2: Load the Coupled MCMC template to setup an analysis.</figcaption>
</figure>

This template is exactly the same as the standard BEAUTi template, but uses coupled MCMC instead of regular MCMC.


### Setting up the analysis with Bayesian Coalescent Skyline (similar to [skyline plots](https://taming-the-beast.org/tutorials/Skyline-plots/))

To import the aligned sequences into BEAUti, use `File > Import Alignment` to select the `*.nexus` file.

BEAUti will recognize the sequences from the `*.nexus` file as nucleotide data. It will do so for sequence files with the character set of **A | C | G | T | N**, where **N** indicates an unknown nucleotide. As soon as other non-gap characters are included (e.g. using **R** or **Y** to indicate purines and pyramidines) BEAUti will not recognize the data as nucleotides anymore, unless the type of data is specified in the `*.nexus` file.

After we have loaded the sequences into BEAUti, we have to specify the evolutionary model. We will be using the very general GTR model ([Figure 3](#fig:model)), which estimates transition probabilities between individual nucleotides separately, meaning that transition probabilities between e.g. **A** and **T** will be inferred separately to the ones between **A** and **C**. Additionally, we should allow for rate heterogeneity among sites. We can do this by changing the Gamma Category Count to 4 (normally between 4 and 6).

<figure>
	<a id="fig:model"></a>
	<img src="figures/choose_gtr.png" alt="">
	<figcaption>Figure 3: Set GTR as a site model. Also use a Gamma Category Count of 4.</figcaption>
</figure>
<br>

As we use sequences that were sampled at the same point in time, we need to fix the clock rate (for more information on this please refer to the tutorial on molecular clocks). We will use an estimate inferred in {% cite Pybus2001 --file Skyline-plots/master-refs %} to fix the clock rate. In this case all the samples were contemporaneous (at the same time) and the clock rate works as a mapping of the estimated tree branch lengths into calendar time.

We will keep the strict clock model and will set `Clock.rate` to 0.00079.

Next, we need to go the the `Priors` tab and set the Bayesian Coalescent Skyline as a tree prior ([Figure 4](#fig:coalescent)).

<figure>
	<a id="fig:coalescent"></a>
	<img src="figures/choose_coalescentSkyline.png" alt="">
	<figcaption>Figure 4: Choose the Coalescent Bayesian Skyline as a population prior.</figcaption>
</figure>
<br>

The Bayesian Coalescent Skyline works by dividing the time between the present and the root of the tree into intervals, thus the number of these intervals has to be defined. Each interval will have a different effective population size.
The Bayesian Coalescent Skyline will estimate the number of coalescent events within each interval (which is captured in the Group Size parameter) as well as the effective population size for that interval. The number of intervals is equal to the dimension specified. If we have {% eqinline d %}  intervals, the effective population size is allowed to change {% eqinline d-1 %} times. To specify the number of dimensions, we need to first go to the initialization panel. This is by default not visible `View > Show Initialization Panel`.

For this analysis we will set the number of dimensions to 4 (the default value is 5). Keep in mind that one has to change the dimension of `bPopSizes` as well as `bGroupSizes`. The dimension of both parameters has to be the same ([Figure 5](#fig:dimensions)).

<figure>
	<a id="fig:dimensions"></a>
	<img src="figures/set_dimension.png" alt="">
	<figcaption>Figure 5: Set the dimension of the two parameters, bPopSizes and bGroupSizes, to 4.</figcaption>
</figure>
<br>

Choosing the dimension for the Bayesian Coalescent Skyline can be rather arbitrary. If the dimension is chosen too low, not all population changes are captured, if it is chosen too large, there might be too little information in an interval to support an estimate of a population size. There are implementations in BEAST of the coalescent skyline that either sample dimensions (Extended Bayesian Skyline {% cite Heled2008 --file Skyline-plots/master-refs %}) or do not require dimensions to be specified (Skyride {% cite Minin2008 --file Skyline-plots/master-refs %}).

We can leave the rest of the priors as they are and go to the `Coupled MCMC panel`.
In contrast to regular MCMC, we have to define a few more things.
First, we have to define the number of chains.
This number should be equal to the number of threads you can run, resp. the number of CPU cores.
The `resampleEvery`, defines after how many iterations, two chains should be proposed to exchange states.
If this is too low, the chains will communicate a lot and the exchanging of states itself can consume more time then the actual MCMC.
If it's too high, the states of chains are exchanged and the ESS per time isn't as high as it could be.
We here propose to exchange states between the 2 chains every 1000 iterations, which will allow for at most 10 000 swaps of states.

The next parameter we have to set is the `deltaTemperature`, the higher this value, the hotter chains become and the more time they spend in unlikely parts of the posterior probability space.
Hotter chains on the other hand are more easily able to cross unlikely intermediate states and can therefore help chains to move out of local optimas.
Here, we use a value of 0.05.
This value should be different depending on the dataset, the analysis and the number of chains.
Overall, it should be chosen such that the acceptance probability of an exchange of states between chains is between 0.25 and 0.6 {% cite altekar2004parallel --file Skyline-plots/master-refs %}.

<figure>
	<a id="fig:dimensions"></a>
	<img src="figures/CoupledMCMCpanel.png" alt="">
	<figcaption>Figure 6: Setting up the parameters fo the coupled MCMC.</figcaption>
</figure>
<br>

### Running the xml


Next, we can run the xml.
The output to the screen of a Coupled MCMC run looks slightly different then the one of a standard MCMC run.
The column called *sample* describes at which iteration of the coupled MCMC we are. The column *swapsColdChain* denotes how many times the one cold chain (the chain that runs just like a regular MCMC chain) has been swapped with another chain. The *swapProbability* denotes how likely it is that a swapping between two chains is accepted. This vaue should be somewhere between *0.2* and *0.6*. A low values indicates that the heated chains are running too hot and are not efficiently exploring the posterior. A too high values indicates that the heated chains are not running hot enough and are thus exploring parameter space that are too similar to the one of the cold chain.

```
sample    swapsColdCain    swapProbability
10000    0    0.0 --
20000    1    0.5 3m15s/Msamples
30000    1    0.3333333333333333 2m56s/Msamples
40000    1    0.25 2m34s/Msamples
50000    1    0.2 2m29s/Msamples
60000    1    0.16666666666666666 2m24s/Msamples
70000    1    0.14285714285714285 2m22s/Msamples
80000    1    0.125 2m20s/Msamples
90000    1    0.1111111111111111 2m15s/Msamples
100000    1    0.1 2m12s/Msamples
110000    1    0.09090909090909091 2m9s/Msamples
120000    1    0.08333333333333333 2m8s/Msamples
```


### Exploring the results of Bayesian Coalescent Skyline analysis

For the reconstruction of the population dynamics, we need two files: the `hcv.log` file and the `hcv.trees` file.
The log files contain the information about the group size and the population size.
The group size specifies how many intervals are combined to have the same effective population size.
The two files are logging the states of the one cold chain.
There are however two more files `chain1hcv.log` and `chain1hcv.trees`.
These two files, log the states of the 1st heated chain and are not needed for the post-processing.
If we look at the inferred distributions between the two files `hcv.log` and `chain1hcv.log`, we see that these don't match up and that `chain1hcv.log` explores less optimal states with higher frequency.

<figure>
	<a id="fig:dimensions"></a>
	<img src="figures/HeatedChains.png" alt="">
	<figcaption>Figure 7: The heated chain does no correctly explore the posterior probability space.</figcaption>
</figure>
<br>



After the runs have finished, load the finished `hcv.log` file into Tracer. Alternatively you can use the `hcv.log` files and the `hcv.trees` files you downloaded with the data. To run the analysis, open the finished `hcv.log` file into Tracer, then go to `Analysis > Bayesian Skyline Reconstruction`. From there open the finished `hcv.trees` file. To get the correct years in the analysis we should specify the `Age of the youngest tip`. In our case it is 1993, the year where all the samples were taken. If the samples were taken through time, the age of the youngest tip is the time when the most recent sample was taken. If you now press the `Ok` button, the reconstruction of the past population dynamics will be performed ([Figure 8](#fig:trees)).

<figure>
	<a id="fig:trees"></a>
	<img style="width:50%;" src="figures/open_trees.png" alt="">
	<figcaption>Figure 8: Reconstructing the Bayesian Skyline plot in Tracer.</figcaption>
</figure>
<br>

The output will have the years on the x-axis and the effective population size on the y-axis. By default, the y-axis is on a log-scale. If everything worked as it is supposed to work you will see a sharp increase in the effective population size in the mid 20^(th) century, similar to what is seen on [Figure 9](#fig:skyline).

<figure>
	<a id="fig:skyline"></a>
	<img style="width:50%;" src="figures/skyline_analysis.png" alt="">
	<figcaption>Figure 9: Bayesian Coalescent Skyline analysis output. The black line is the median estimate of the estimated effective population size (can be changed to the mean estimate). The two blue lines are the upper an the lower estimates of 95% interval. The x-axis is the time in years.</figcaption>
</figure>
<br>

There are two ways to save the analysis, it can either be saved as a `*.pdf` or as a tab delimited file. To save it as a tab delimited file, you can go to `File > Export Data`. The exported file will have five rows, the time, the mean, median lower 95% interval and the upper 95% interval of the estimates, which you can use to plot the data with other software (R, Matlab, etc).




## Setting up the same analysis without using BEAUTi

In order to setup the analysis to run with coupled MCMC without using BEAUTi, we can create and xml that is supposed to run with regular MCMC.
After this is done and the `*.xml` file (here `hcv_mcmc.xml`), has been saved, open the `*.xml` file in a Text Editor (e.g. TextEdit in MAC or notepad, but not word! in Windows).

Next, go to the line that contains the code `id="mcmc"`. In `hcv_mcmc.xml`, this is going to be the following line

```
<run id="mcmc" spec="MCMC" chainLength="10000000">
```

Next, replace that line with
```
<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" chainLength="10000000" deltaTemperature="0.1" chains="2" resampleEvery="10000">
```

*  `chainLength="100000000"` defines for how many iterations the chains is run

*  `deltaTemperature="0.1"` defines the temperature difference between the chain *n* and chain *n-1*.

*  `chains="2"` defines the number of parallel chains that are run. The first chain is the one that explores the posterior just like a normal MCMC chain. All other chains are what's called *heated*. This means that MCMC moves of those chains have a higher probability of being accepted. While these heated chains don't explore the posterior properly, they can be used to propose new states to the one cold chain.   


Next, save the xml file again and run it with BEAST. In order to post-process the output, just use the log files that do *not* start with `chain...` since these are the log files for the heated chains.


### A note on resuming coupled MCMC runs

Coupled MCMC runs can be resumed just like any other BEAST.
It is however possible that the different log files of the different chains are at different iterations, which will return an error if you try to resume this run.
If this error happens, load the log files of the cold chain and all hot chains (the chains that start with `chain...`) and look at the lowest value of the iteration.
Next, open the log files of all chains which had a higher iteration and remove the last lines until the last line has the same sample number as the lowest iteration of any chain.

----

# Useful Links

- Blog post explaining coupled MCMC [https://darrenjw.wordpress.com/2013/09/29/parallel-tempering-and-metropolis-coupled-mcmc/](https://darrenjw.wordpress.com/2013/09/29/parallel-tempering-and-metropolis-coupled-mcmc/)
- Coupled MCMC source code: [https://github.com/genomescale/starbeast2](https://github.com/genomescale/starbeast2)
- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file AIM-Tutorial/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users)

----

# Relevant References

{% bibliography --cited --file AIM-Tutorial/master-refs %}
