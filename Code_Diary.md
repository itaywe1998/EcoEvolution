---
editor_options: 
  markdown: 
    wrap: sentence
---

## Historic Progress (Code Diary)

> This is a segment of free text which address the work difficulties with ecoevo from the start at early June 2023.

#### Introduction to original data structure

```         
In main file (ecoevo_main.R) line 26 indicated there are 4 optional models for this code."baseline", "trophic", "Tdep", or "Tdep_trophic", I will look them up in here

Somehow got negative densities on 30 degrees / 1e5 yrs , baseline model, 5 species, in the after_cc there are 15 out of 250 initia entrees for which n\<0 ??

In Tdep, tstart = -4000, during = 2E+6, after= 2.5E+4 , S=5 & L=20, time was 2 minutes There were no differences between results for during = 2E+4 and during 2E+6 -\> only temprature depndent, meaning species arrived equilibrium for during = 2E+3 they havnt reached equilibrium yet (change too fast for default vbar and dbar, 4 out of 5 speceis went extinct by tE).

So vbar is quite sensitive.
when reduced from 1e-5 to 1e-6 the gaussian on the mid-latitude went much narrower, and when reduced to 1e-7 the only before-thriving species now had to find shelter in the poles, slowly coming back towards the equator when the temp rise stopped.

Now the code starts to make sense!
Tdep is therefore a favorite mode to work with.
I guess this is Speciation time, since I invented my 5 initial species and went beyond abiogenesis.
For that matter if I will generate just 1 species he will start at the poles and most certainly wont make it.
But let's check.
Yep, he died.
Poor lad.
But hey first overall extinction !
Congrats!
I might need to work on the generation on a random patch rather than the poles for a single inital species, or even get the chances higher with the equator.
It will be nice to think soon on updating some real tempratures other than Earth poles vs equator gap.
```

#### Speciation and Speciation Event implementation

```         
It is possible to use the m values on dat (mean trait values for each species in each patch) We can make a Speciation Event Call (SEC, now i came up with it) when the \|(m(i+1) - T(i+1)\| \< \|m(i) - T(i)\| when i is the origin patch out of the two neigboring patches.
When this happens n(i+1) needs to be re-associated with the new species.
We assume the ones that did not adapt will eventually die, and so we leave n(i+1) for the ancestor to be zero for the moment, and let the spatial dispersion do it thing with the traits from outer patches.


Alright, Temp gives the temp over all locations at a given time inside rhs_eval.
I need to find the right place to : A.
Check serially if the SE occured B.
Update all variables with S as a dimension and update S C. If the ancestor to die in this patch is the only one left, then must S-- (even though not very realistic since n\>0 everywhere, need to think if it matters) D.
Associate n with the new species

Updating species number S from previous timestep From the moment will be performed here, might be in the end instead Might not matter on large scale 

OK, so the implementation has shifted to fully utilize state vector = dvdt = ic The thing is (this is a vector which stores in the first S*L places the densities(n(i,k)) and the second S*L for the traits (m(i,k)) now that we modify S, each state vector will be in a different lenght, which is fine till the part of getting it all in one matrix (before orginize data) or datatable (after orgdata) will figure this out!!!

debugg log (with much prints in cpp): state size is 121 S now is 3 dvdt size is 121 state size is 121 S now is 4 dvdt size is 161 state size is 121 S now is 5 dvdt size is 201 state size is 121 S now is 6 Error in evalq({ : Index out of bounds: [index=3; extent=3].

So, as it turns out ode (from desolve package) accepts only constant dimensions to its state vector, which i guess is determined by the initial condition vector.
For this reason, I must go on plan B - which is predetermine maximum species number, and padd with zeroes, that is , always have the longest state vector possible, with dynamical S.
Let us have for debugging manners an initial MAX_S = 100, so we will deal with relatively short runtimes.

Good, now we have origins column, and a MAX_S oriented state vector.
I am going to start speciation now.
I had an idea about a log file, relating to ancestry and death of species to follow the simulation history.
Extinction will be noted by dat\$o = 0 for the relevant species, ancestry will be reported dynamically on the creation of a new species.
Origins must be kept as data and not log since it is used to determine a SE.
The log is nice to have and not top-priority right now, hopefully will be implemented soon after speciation is complete.

Well, I guess that for simplicity, I will treat MAX_S as the maximal amount of species !!!T
HROUGOUT THE ENTIRE HISTORY!!!
the alternative is bad.
reindicing must be perfect in order to not lose continuity in the ODE.
Another dificulty is to keep track , for instance: There are 7 species, and species #3 died.
When the next species is to be born, and we want to save space, we designate it as the new species #3, so, from the simulation perspective it seems the old one came back to life ?
We do document it on log file, and the origin field will probably be different (what happens when its not?) but it will be easy to confuse when dealing with a large quantity of species in the future.
It means our data will carry the history throughout the process, and the display will be a bit less costly.
It does limit the code efficiency, but for the moment we will go on the tidy option, and optimize later :) In this implementation, S only grows, never decreasing.

Another thing to do now, before Speciation kicks in, is to review all the params which are S dependent, and probably make them fatter :( But then again, unlike m,n,o - These variables are not te be altered, and maybe there is a way to pass them in a different way. For a no-consumers scenerio (god forgive me for the day I'll need to get consumers done) these are: (v, d, rho, a)

These are not calculated but rather randomized from a uniform distribution !
So it cannot be done inside eval (cpp) since it is non-consistent.

Ok, let us have a settlement.
It was crazy enough to make the vectors this big, but since they can theoretically can really fill up, it has an excuse.
On the other hand, some random variables from a uniform distribution is not important enough so we generate 1000 of them and carry with us on all the vectors (like origins), although the vectors are MAX_S*L+MAX_S now, and the suggestion will enlarge this to additional MAX_S* MAX_S from a alone, and additional MAX_S from each of the others.
THIS IS UNACCEPTABLE!
So the compremise is as follows.
I will generate a generous size for a,d,v,rho which will be 20.
If the species number is greater than this SIZE, then the effective index will be i%20.
I assume 20\*20 matrice or 20 long vector is diverse enough, and species will face different challenges (m and n) so the behaviour will not feel too repetative.
If repetition is felt, then a larger SIZE will be set by me in the future.

Apperenatly, the seed is determining the values of each randomization in the run, that is why for each run with SIZE = 20 (now noted as P in the code) I will get the exact same values, down to the 8th digit after the dot, in each value of each vector and matrix in (a,d,v,rho) When P is changed , the randomization changes and so the vlaues, and henceforth the simulation yields different results which are visually different in the plots.
So - Do not panic!
it is technically ok since nothing is wrong with a certain configuration or the other, just remember not to change the seed every time because then it will be hard to tail bugs and mistakes, since the results are very seed-sensitive.

Well, I am deep into speciation implementation.
All that is (Theoretically) left is to smooth the density(n) in the place which was previously occupied with the old species, and create the gaussian denisty profile for the new one as well.
But - it is in cpp.
I am going to write a function without ability to debug?
maybe i should take a time-out to a different script or something, just to see if it behaves as intended.

Debugging the cpp with prints and constant exceptions is annoying, but it seems that I covered much misfunctionings in the code so far.
Needs more runs to get to working code.
Last thing I added was dvdt[place] = origins[i], which was missing before.

Hi there!
Im back.
For some reason the timestamp for NewS -\> above MAX_S is doubles around -4000 (initial time!) so the population explosion allegadly happens instantly.
I will assume this is not a real data driven case but wrong code writing by me.
In this point print debugging on c++ file doesn't help no more, and we shift to 'gdb' debugging from the linux console.
Good Luck!

Got a little bit tires of preparing a debugger for cpp (emacs or other option) so i figured out i will get back to some new print debug.
In the process of remembering why I blocked the code from enlarging S and putting an exception for growth S \> MAX_S under the first time step !
(time = -3999.5 ) I printed the origins vector, and as a part of the state vector it seems to change values, from [1 20] to \~ [1.03 20.61].
Id this is the problem, it is bad since it means I cant use the origins vector as a state variable and might have to calculate through a functions like Temp, and not sure how it can be done.

Well, several observations so far.
1.
the dvdt vector is the deriatives (rate of change) for the states vector.
Thats why origin=20 changed faster, because you also put its deriative to be 20.
2.
The eqs is printable in several points arount t=time.
The if that controls the n.m change around a speciation (line 141) was not done right because now the prints in there are not activated.
3.
For the former reasons (no update at neither o and n.m) the speciation occured repeadatly around a true m change for the first species in the second patch.
The condition might be too easy, or the genetic variance too large.
This is a secondary effect and partially a scientific and not technical question.
For those reasons and I am sure that more that yet to be detected, the S runaway (S -\> \> MAX_S) was possible in this random seed.
I believe the first step is to move origins outside of the statets vector since it cannot be altered in the desired/intuitive way at all.
There is an option in ode function for substitution rather than rate change, but then I lose the entire logic of the code for m and n, which is a practical disaster.
I am considering a file to be called from R and managed by cpp, or maintaining some kind of a global string that encodes the information in some sort.
Both sounds not so good, but I need to see whats faster in runtime.
For the record - I did try to assume that 1 change rate = +3% change, and throught this assumption micro-manage the origin values to desired, but it failed, at least at the micro step that the simulation did not crush in (t circa tstart=-4000).

VERY VERY IMPORTANT WARNING I realized my intention to assign new values to n & m after Speciation event are not happening in the present code, since I only change the n & m locally , not in the states vector.
The thing that does get changed is the deriatives, such that instead of spliting , lets say, half-and-half a certain parent population, I just reduce the deriatiative such that it follows a n/2 value (not lineariliy i suppose) and invent a new deriative to a population zero of a new specie.
It is not bad, but its not what I meant, and its certainly not a discrete event, but a process.
It might be not very awful, at second though, since the bdf_d ode soving method for stiff equations takes a tini-tiny step of several digits after the zero of a time period (0.001 years).
The time steps I asked for are just the savepoints, the documented entries in a table, a graph resolution.
The true information is much more refined which makes much more sense.
At least I want to believe this is true, and I will be happy to understand more about it in the future.
This means that if I dont change anything and make the simulation work, I might accidently do better than intended since the solver will stay continuios, and in such equations it might be the only option for solution.
On the other hand , there is an option that I will go for an experiment on a copy of the directoty before modif (not the original from stockholm, the one I finetuned to my needs).
There I will try to recreate the model success through the "iteration" option of ode, and so directly assign values to m & n (and yes, in that case o is allowed back in).
The Very Great Risk in this is that I will fail to perform a reliable alternative under a sane time period, and it feels to me very wastefull.Altogether it might not work as I intend to break continuity on purpose.

Well , I guess I can forget about direct substitution , since in a simple print for the time in each eqs call - the timesteps are not equal!
for instance : time is 34589.6 time is 34589.6 time is 34591.5 time is 34595.3 time is 34595.3 time is 34599 time is 34626.3 time is 34626.3 time is 34626.3 time is 34653.6 time is 34653.6 time is 34653.6 time is 34653.6 .
. . . time is 43881.8 time is 43832.8 time is 43832.8 time is 43849.1 time is 43881.8 time is 43881.8 time is 43914.5 time is 43914.5 time is 43914.5 time is 43914.5 time is 43914.5 time is 43990.6 time is 43990.6 time is 43990.6 time is 44066.8

Whats is this ?
Nothing I can predict, so i don't know to directly translate the dvdt to dv since dt is changing?
Something is off.

So , assuming I am bound to use the more natural way of change rates, maybe it is a good thing.
I can decide that Speciation Event will be a turn-on event of a process with positive flow to new species and negative flow to the father species.
I would like to believe that new species starting with population zero will take time to evoke another SE in a different (probably neighbor) cell.
After this will work I believe there is a place to compensate on the simplicity of the model with harsher SE condition (that is , a certain trait distance abs(m(i,k)-T(k)) \< Threshold in addition to the existing condition, this will further signify that the candidate species truly adapted to its local enviroment.
In later stages it will ofcourse get harder to provoke a SE with futher complexity included in the model (additional trait types, bringing back predators layer, maybe in a completly follow up simulation, we'll see).
First thing, I still think origin is a vital information in order to keep track on the species indentity.
First thing to be done is get origin out of state vec, and implement an alternative.

There is a big question mark here if I want to keep the exact species identity at the end , or just see if there is anyone who is still alive.
If the latter option is the case, then I can use my space more efficiently and re-orginize a species to the place of a dead one.
In tables there are no worries, but in cpp level it is a pain in the arse to move everything, and prone to error.
I will need to see in the future if my datatable is mostly sparse and if I run to complexity issues, then it will be the first hting to do.
At the moment no ne dies so it is not relevant yet.

Pay attention that in current implementation, split only controls the favor in growth for the descendant speices int the very first step of its creation (that is, to start up some gentic migration) In the next call of eqs the groth terms will be calculated by the n itself, which will probably will not shift much from full population for ancestor and zero for descandant.
For the descandant that aint much of an issue, maybe a harsher condition to really come to surface.
The problem is that the original species will not vanish neccesarily, and so the S will boom with tiny species all around.
We will see.

After short consultation with the gen3is rep in github (<https://github.com/project-gen3sis/R-package/blob/master/R>) under files gen3sis_species.R & speciation.R, in create_species_from_existing in the species script, lines 67-68 indicates that my original thought was alright, that in the cells of birth of a new species n(son)=n(parent),and the same for m.
Sadly, I could not directly substitute in my case, and so i must think of a way to stay in an event or re-labeling format.

All for the worse, I have a difficulty to copy m, since it start with 0, I will have an entirly wrong trait for the newborn species!!
VERY CRUICAL!
I have no choice but to dig in deSolve src, and find the way dt is determined.
With that knowledge and the hope that dt is somehow constant per run (now from the deSolve docs I am sure that the "times" input is only the output displayed, not the entire timepoints calculated. That is why the results are unchanged by "times". When I find dt I will activate the "iteration" option in ode, and sustitute new values instead of groth rates. I think its for the better to recreate results with directories without speciation, so I will know the code itself wasn't buggy before I started to change it. We will have a new directory for that.

First of All - Congratulations on your very first working Speciation-Ecoevolutionary code!
It is preferable now to work through the correct parameters and time settings.
Not the greatest fear of mine.
I remembered that the original paper by Nornberg mentioned 'generations' and the code gave an option to plot_richness which occured a suspect that maybe I missed a speciation code in the original one and was doing nonesense.
In the plot_richness function under plot_functions.R file line 247 the "richness" is just the different densities segmented by their local zone (temperate/polar/tropic) with the condition that densities \< 1e-6 does not count, and so richness may shift (started with 50 species with max density 1, all of them became sparse and some "vanished" according to this count).
```

#### Back to original - looking for extinction vs survival cases

```         
Hagai Plan: The ORIGINAL CODE is good enough!
Maybe one day we will use the Speciation option, for the moment the addition of unknowns is not compulsary since we hav enough of them.
The runs to do are: With and without time to adapt: Generate 1 species at one source (or 10 identical, see a couple of configurations all at once) , and : 1.
Give 100M yrs to disperse and adapt, then perform climate change (say 300M yrs) , of 30 degrees.
2.
Give minimalistic time (<1M yrs), so it doesn't disperse much, then perform climate change.
Arrange all the constants so the behaviour deifference will be on this scale.
Obviously needs to change methods to "rk4" or something so I can control the dt (I don't want the simulation to run forever).

As for the Thesis Suggestion , I should write in 0.5-1 pages in word (english) what is the aim of the research, what tools (the EcoEvolution code) and what we intend to do.

So far I had a nice chance at balancing the model parameters, but its yet to be good.
I will have to suffice with the adaptive timestep, since rk4 goes negative n(i) for step > 0.01.
For tE~1e5 it runs ~2 minutes.
Runtime dependency is not exactlt linear but not very far from it, so for tE~3e8 it will take ~6000 min = 100h = 4 days!
not worth it

The rtol and atol was nice.
Got to speedy computation on the right timescales, with baseline model.
But, when ic heck for the random seed, aperantly the ode success (getting to the end) is sensitive to seed choice.
I will continue with the baseline to make the model immune to this sensitivity.

Back again (19.09.23) after a while.
For comparison, I got access to Nornberg's code in 5th of June 2023.
So I figured out that the original code by Nornberg is not generic in reference to the initial mean temprature over the planet, that is, for a long simulation with Tmin = Tmax = 35: (add image from screenshots) the population at t=0 is 4 times larger than initial time, compared to long simulation with Tmin=Tmax=15: (add image from screenshots) where the t=0 actually shows a decrease in population compared to initial, up to several tens of %.
The reason must be connected to hard-coded numbers in the code.
It might come from the constants in smoothStep, or it probably originates in the constants main sends to eqs.
It needs to be cleanly tested.
Also - it might be the right time to move to cleaner platforms.
It should be a good idea to invest a little time in how to push to git some updaes since i change the code all the time, and some versions reconstruction might come at handy.
As for this notes document, it can move to a markdown format, where it can be divided to sub-parts for readability.

In order not to break the format of keeping track of all deprecated , I decided for now to even keep the note to create this new markdown!
Git is up and running, and probably the Temprature dependency is in aw parameter. When set this slope to 0, the dependency is seemingly handled, but further runs on full use case are yet to be done.

I believe this code is still small so a proper doc is not needed yet, but it can be a good practice to update the README file. Plus, thesis suggestion is confirmed (2.10.23).
Update : aw = 0 consumes breaks down the R Session (memory consumption) in the case of  Tmax != Tmin. Should check why when uniform temprature there is no effect and what is the meaning of aw parameter.  

It appears to be (04.10.23) that the first stage is non-varying, and so it is fast. The more gradients are present in the ODE, the longer the calculation takes, and smaller steps are required. 
Now, an "interesting" development may take (even though its not the best way rate such a thing) few minutes for the common configuration we currently have. Changing run id to "MakeBeforeVarying" , because my current interest shift to solve this issue as a key to achieve cases separation.

Both images under 'ecoevo_original/noted' demonstrate the same thing (maybe id2 was wrong, will check) that adaptation time does not change the final result. Could see what is the difference in evolve towards this final common state.

Carry on from ecoevo_original/ecoevo_main.R

Congratulations!
The simulation id "Promising" bear fruit. For short adaptation time exctinction occurs, and for long adaptation time life prevails.
Several comments:
1. What I did not realize before qas that I did not touch every single variable of the code, and the specific relatively constant scale was the Climate change time compared to the short adaptation time. 
Once scaling them both with the vbar and dbar through scale variable `y` , this degree of freedom was utilized.
2. This config is not the most general one. Competition was reduced as far as possible, dbar and vbar identical across all species, and most of all, the CC (Climate Cahnge) was not equally drstic through the entire planet, and a temprature range which existed at t=0 remained (somewhere else on the planet, but still) up until the end.
Namely for the first success case:
Tmax <- 15.0 # initial mean temperature at equator
Tmin <- Tmax-30 # initial mean temperature at poles
Cmax <- 30 # projected temperature increase at poles
Cmin <- 12 # projected temperature increase at equator
Specifically , for CC duration of 2e6 years, all adaptation times smaller than 1e8 years were extinct


Little Remainder : There is no point in setting different random seeds if all "runif" in the code bound to equal min and max. Runif is a random number sampler, with syntax - (sample_num, min, max).

Multiple runs ares still wanted, now that the format is clearer (e.g):
Rscript ecoevo_main.R Tdep False 2 bashTrial

I guess check for more interesting temperature range.Than add complexity step by step.
For some reason when CC is uniform, the population goes to inf, crashes R. I wonder why is that so and if this is a sensitive case or a more general one.

writing a dual run bash (small and large automatically one after another ) is good, but a dual demo is also needed. I want to show the small time failure , and so a certain data needs to come from there. Currently when throwing an error the data up to this point is lost. Should think of a way to save it.
Near future (08.10.23) task - find a family/neighborhood of inputs so that case separation occurs.

(09.10.23) Done with data display for FAILED runs. Now just a dual demo.
Demo ready.Going for little grid search on bash.

(10.10.23) Bash loop on y variable is implemented.[multi_run.bash]

Remember! patch #1 - poles, patch maximal number  - equator.
For reading paramters saved from past runs, use the "load" R command from any R console with the proper file path.


Hagai meeting summary:
Very good progress. For article we also need 
1. periodic climate change. Maybe do few (4-5) periods over entire CC time
2. start technical description in overleaf. Introduction text and mainly the equations behind the code.
Nice!

Let is start with periodic CC.
Periodic CC first attempt, needs to configure a successful survival scenario for long adaptation since now for fully CC shifts (30 degrees) the CC is effectively cycles*2 faster, so needs to answer this time dependencies.

First successful (16.10.23) periodic CC! Survival at long adaptation, mediate CC shift (Cmin=7 and Cmax=15) and cycles=3.
Need to see:
1. Generate cases where short time fail when long succeed
2. Consider sin() and not just abs(sin) for temp profile [in rhs_eval.cpp, periodic_smoothstep function], that is, low tempratures beyond initials

*by the way, good to know that with the updated tolerances (1e-14) the sin profile is a non issue. maybe it is a run duration thing and not deriatives problem.
3.Allow more periods safely, enlarging tE is a good direction but not trivial,
didn't go smoothly at first attempt

17.10.23
Tried to compare (using cpp std::chrono) the speed of periodic (sin) compared to 
regular CC (polynomial). In the single call level, sin was even faster (20ns compared to 30ns in average)
for some reason the runs never ended even after 25 minutes, and the output files reached 2GB (not mega, giga!) before interrupted by user.
Lesson learned is not trying to print each time such a small function is called. Run time skyrocketed for this prints alone.
Plus, there is no way to compare the amount of times of function was called, but by the logic of the rhs_eval function structure, it is supposed to be indentical call number.
So for now the check is enough, and the run time extension is probably not 
due to sin() itself, but the complexity the frequent CC is loading on the ODE solver.

made progress with periodic CC, very prospering and good demo. Just make sure changing the tolerances doesn't affect the final population in order of magnitudes


22.10.23
Learned from bitter huge batch runs, so examined manually the surroundings of eligible configuration.
It appears that changing S is the same as having a complete different "dice throwing", as for each 'runif' function alters the number of samples, and so the entire values.

As for updown, in the configuration relevant to today, the population skyrocketing doesn't change much, and not periodically at all.

The tolerances modifications are important in order for ODE to yield a successful integration, as when the run code = -1 in diagnostics the integration failed, and so the results are not valid (for most of the times failure will occur on the very first step of ODE, i.e for before it will be t = -9.997e7)

In case of adding Predators (I was looking for ways to deal with the ever-growing population in this particular configuration) - it is still (TODO) to be found how to hold a proper integration in this more complex scenario, since no tolerance I tried helps, all integrations failed (code -1).

Now for the most intriguing part. It seems the sweet spot for now is [Cmax,Cmin]=[30,20], For 20,10 the short time survived, and for [50,40] both failed and, more importantly , long time lost its advantage, and both demised at exact same time after CC start (same t).

Fairly interesting. The safe bet now for case manufacturing is just go with the said configuration and alter only the seed.
I guess it is quite boring but it will potentially give several periodically successful cases. Problem is the periodicity does not show with this population explosion, so I think retrieving the past dual configuration is preferred, and go changing seeds from there.

Would like to add an aid - if lsode returns -1, the program entirely is not valid, and exec must stop (something like exit() in C)


Comment :in the end it is not that "FAILED" label is not applied on relevant cases (when small time fails) but the coordination between demo and main script in regard of the vbar and dbar values is just the worst!
Should coordinate between them automatically , or give up the mentioning in the name.

For now , since this is more than just these two variables, the manual coordination will continue, but I intend to dedicate next workday to this alone, since this is risky and highly prone to errors.

Factory is at least functioning, just did not find the grain case yet. It will come soon enough.

25.10.23
Coordinated between demo and main R script through shared inputs from multi_run.bash.
Due to this implementation heavily modified input/output format (up to 12 argument to a single script)
But now all is controllable from multi_run.bash, with option to replace scalars with arrays , if nneded in the future.

Plus discovered a serious BUG in R main - periodic was set to false during mid process, no reason I can think of it will be justified.Removed.


27.10.23
Had a bug when adding 'plyr' library to the include section. Did not see the warning regarding its use with dplyr (needed to arange them with plyr first).
This is an example to why supress warning is not a good idea, unless you already trsut the package.
For next time - do include the library outside the SupressWarning block, and only after check you can move it inside.

For the good news - we have a first case for factory ! Will run now on multiple seeds, at the very least. Then try to find additional cases.


In manual run (Rstudio) all is well, but fro some reason when passing the lat two time arguments through bash scripts:
// Bash echo $@ output
3690 PeriodicCC_Factory_cycles4_seed3690_AbsSin_25_10 .00003000000000000000 .00000010000000000000 4 FALSE 25 10 -100000 -100000000 2000000
// ecoevo_main.R print(clargs) output
 [1] "Tdep"                                            
 [2] "FALSE"                                           
 [3] "3690"                                            
 [4] "PeriodicCC_Factory_cycles4_seed3690_AbsSin_25_10"
 [5] ".00003000000000000000"                           
 [6] ".00000010000000000000"                           
 [7] "4"                                               
 [8] "FALSE"                                           
 [9] "25"                                              
[10] "10"                                              
[11] "36900"                                           
[12] "36901"    

which are following the seed for some reason. Problem was avoided with demo since I found a way to not have demo necessary informed from bash, but from the data loaded.
If no need to coordinate between demo and ecoevo this things, maybe it can go back being an internal ?
But I still want to keep the capability to easily manipulate them in different runs.

Three possible ways for solving:
1. Check for plyr blame - re implement the round_any function, shouldn't be long. See
if without the package things go back to normal
2. Try to write the relevant times to sys.env variable, and write/read minimally.
3. Keep investigating 

30.10.23
Done 1 - not the solution.
Investigating - for some reason when run.bash gets to the 10th argument things start to go wrong.
May use a single string to describe all three times, in a format "[small,large,final]"
Maybe in the future will make a more compact argument list so won't cross 9 arguments.

Ok, now it is working with manual parsing, both bash-bash & bash-R.

31.10.23
I am quite satisfied with the current examples of periodic CC.
It is however important to keep note for future case production for how to manage the current configuration in order to keep getting successful results -efficiently and without unnecessary searching time.

For the configuration 
vbar=$(bc -l <<<"3/(10^5)")
dbar=$(bc -l <<<"1/(10^7)")
Cmax=25
Cmin=10
small_time=$((-1*(10**5)))
large_time=$((-1*(10**8)))
final_time=$((2*(10**6))) 

(all is abs(sin) and not full sin)
One can have 3 or 4 cycles with good results.
At 5 cycles the CC is too fast.
So, for instance a solution like this can be done:
  vbar <- 3e-4  
  dbar <- 1e-6 
  cycles <- 5
  updown <- FALSE
  Cmax <- 25 # projected temperature increase at poles
  Cmin <- 10 # projected temperature increase at equator
  tstart <- if (small) -1e5 else -1e8 
  tE <- 2.5e6
  
Where the experienced user will notice a slight tE extension, but mainly a uniform rise at the vbar,dbar.
One could also try lowering the temperatures (Cmac,Cmin) but not directly successful like the aforementioned suggestions.

Another thing is to pay attention to difference between went extinct but run , vs failed to compute (code -1 from lsode).
When code -1, try to make the tolerances smaller. When too small (slowing down the process considerably) or even too much steps were taken (up to MAXSTEP) consider moderating the CC with ny of the suggestions above.

Even so, sometimes it will not work out, but if you get stuck under many attempts, lower the species number, it will hasten run times and make the process more accessible. 
```
#### Towards First Paper & Kozai-Lidov Oscillations

```    
05.11.23
Hagai meeting summary:
1. Define in the introduction what is habitability.
2. Define in each of the equations in methods what is the meaning of each of the terms (heratability for instance)
3. In the introduction needs to explain the idea of our research
4. In methods needs to present the program (tool) i am working on
5. Attach some graphs as is to results section, then we will see which ones to chose from
6. Summarize in a table what series of runs we did precisely, what are the inputs of the model and what range of each of them was tested
7. *Code Related*: Find successful cases also for Kozai Oscillations - average temperature from real systems (outer simulation) , will be a nice addition to the paper
8. Format - for now aim for APJ (Astrophysical Journal). I guess it is for my own good but he said it in the last minute so it is not top priority.

11.11.23
I have worked but did not document.
In the middle of typing the equations from naoz 2016 appendix for the secular equations.
Need the other 7 out of 8 state variables rate equations.
Good luck and take it easy, it will take some time to make sure its all correct.
* I believe all ^2 in trigonometric terms mean (cos(itot))^2 and not cos(itot^2), it just does not make sense to have the square of the inclination.

24.11.23
got periodic eccentricity, but changes are very small and with periodicity of ~60K yrs, and so the temperature shifts between -107.35 and -107.45 Celsius, not very interesting.
The eccentricity changes between 0.22 to 0.205, which is more than the error on the initial value (+-0.03) which is significant.
I might have the chance to see a proper process when tweaking the numbers.
Should consult with Hagai if taking a real system is good or take other parameters.
For now let us try to find more significant configuration.

25.11.23
When running data from naoz2016 figure 4
large eccentricity differences where found, but on e1 !!!
that's in accordance with the article, which means i have to find a configuration for highly changing e2 from a different place.

We kind of have a design problem.
The significant lidov-kozai oscillations are on the inner test particle, a.k.a e1. For that to be a planet, one must place a nearby host star, and a far perturbator. 
If the m3 (perturbator) is the most massive, then it will emit most of the radiation and will dominate the heating of the planet. On main sequence stars for instance , L is proportionate to M^3.5.
Stellar BH are even more luminous from accretion disks.

Maybe we can consider a sun-near planet-far planet relation, like the sun-earth-jupiter setting.

In this current config - m3 is the parameter who controls the cycles - that is how many periods of temprature profile there are in a certain time range. Make it smaller for fewer cycles, and vice versa."
[main 7cb92fd] in this current config - m3 is the parameter who controls the cycles - that is how many periods of temperature profile there are in a certain time range. Make it smaller for fewer cycles, and vice versa.

26.11.23
proud to say kozai is working.
I opened a new directory in order to massively adapt the main R and cpp scripts to an outside temperature profile.
The attitude is now completely different. There is no before climate change and during one.
There is no small or large preparation time.
The entire temperature history is given by the astrophysics from the start to the end.
There are several modifications to be done:
1. Cpp needs to have access to all temperatures in any given time, even though it chose its time jumps on its own. This is going to be a nightmare, since besides the accessibility to a non R script (passed as parameter is the best solution) one must match the times inside eqs to a time given in the Temprature table. I guess a round to nearest will be fine as long as the table is with enough resolution from one side, but not to heavy on the other side.
2. Many script arguments are now deprecated, such as small (boolean), cycles, periodic, Cmax, Cmin, Tmax, Tmin, tstart, tE ...
Keep only model, dbar, vbar, id, seed.
Less complex on this aspect.
3. Need to adapt Cpp not to generate its own temperature profile, but to read from a given one, in the manner described in clause 1. The only profile to be made is the spatial interpolation, which in turn can be also better described than the linear one if a scientific source will be found. For now its good enough, but patches are still an inner choice and not astrophysical.

P.S dont forget to do source to kozai in the new ecoevo_main.R

28.11.23
Wrote a sketch. Need to be checked.Only ecoevo_main.R and cpp, not all the bash around and not demo. Plus need to see about plot, not decided yet.
01-02.12.23
successful kozai - even too good because our model is well trained in temperature deltas bigger than a few degrees, and the temperature change rate is surprisingly not challenging enough but easy for the populations. Will try harsher situations
longer kozai was performed, fully adapted. Will run lower vbar dbar to see failure or change the profile to be denser. Need to consult Hagai on the matching cases in the kozai context

02.12.23
From the inline documentation of kozai.R to the kozai() function:
A little "Controls & sensitivities" Guide
   1.if desired to decrease the gap between the pole and equator - increase obliquity (eps)
   even possible to revert the relation above eps=45 degrees, making the pole hotter than the equator
   2. the Jupiter-like perturbuter's mass (m3) affects the order of behavior for the eccentricity (e1)
    certain values will generate a periodic (ordered) pattern, others will create a rising frequency
   3. i1 (inclination for planet) heavily impacts - increasing makes the eccentricity much more extreme and hence the temperature difference over time
   4. lowering i2 as well depresses to a well ordered and low difference eccentricity, so it is the combination of both inclinations who rules  
   5. Additionally, i1-i2 matters, the larger it is the less moderate the system (e1) is
   6. As planets' mass tend to be lower compared to the other masses, it is more vulnerable to extreme changes
   the bigger it is the more ordered behavior occurs. In small masses - rising frequency
   7. Omegas don't change much (as expected), just the initial point of the same behaviour
   8. initial e1 changes the stability of the e1 solution, leading to various resulting 
   trends, not very linearly responding.
   9. Star mass does not change much of its own, but in order to keep sense,
   the mass-luminosity relations obliges to increase much of Ls with minor changes to Ms
   10. In a similar manner - a1 also controls the average heat, similarly to Ls
   and so both parameters are recommended to stay the same.
   11. Albedo is a rather straight forward tuning to average heat 
   12. the more distant (a2) the perturbater is the more ordered behavior 
   extremely close Jupiter-like breaks the solution (as expected)


So now - we have two options to simulate the parallel to "long vs short" adaptation time:
micro-scale : The frequency of temperature change is high vs low, over the entire profile
macro-scale : a chirp-like profile, with rising frequency, but a small vs large rate of frequency rise - or in simpler words, starts ok and then gets crazy in both, but for one the crazy part gets much much sooner.
I believe keeping the 5e8 years is a good amount of time since it is biologically significant and not too slow to compute on its own (just T profile , not the ecoevo model, which is also to be considered when enlarging tE)
Plus it is comparable to the temperature history we know on earth.

Until I discuss this with Hagai, let us first try to deal with:
1. A certain T profile with large T difference over time and smaller Tgap between pole and equator (already achieved, now need to find a configuration in which ecoevo survives)
2. The micro-scale comparison (pick a T frequency in which it dies vs one in which it lives)
3. The macro-scale comparison (slow advancing chirp vs fast advancing chirp)

A certain difficulty I am having is this - the indicating difference (maximal single step (2 neighboring points on profile) temprature rise/fall) is mostly close to the max - min of the entire profile. that is - the maximal difference in temprature occurs in a single step ! Maybe even in less than that.
It is not possible to expect the ecosystem to survive a total on-off situation in which the maximal CC occurs in a very short timescale. On the other hand, the ecosystem should be challenged by a bigger max-min difference than the Tgap between the equator and pole in order to show genetic adaptability and not only dispersive. Maybe in a different obliquity it will be simpler to get this cases.

Achieved a rather ok profile for kozai, but still extinct. should modify vbar to find the value which allows survival. 
Another thing which is annoying - the ggplot save is quenching the plot to be very small and not details can be seen above 2 subplots together, should figure how to improve this.

09.12.23
Hi I'm back again. It was a display problem, see dpi, height and width options for ggsave().

Ok, so for v=300 and d=1e-5 they managed to be not bothered by the kozai CC at all ! (talking about the
"First_lowGap_highDelta_Succeed" configuration). Proves that this is still only a numeric consideration and not something about the kozai process itself. Now let us lower the bar to find where does the minimal requirement lie.
The bar is circa 50, but results are boring - constant on the very first step (5 million years in this case)

Well (on v=30) turns out it was a matter of tolerances! Decreased from 1e-5 to 1e-8 each and doubled the fail time from 2.3e8 to 4.6e8! Further tuning will help me yield proper solutions.
TODO : Find a reasonable compromise on performance and tolerances, reached 1e-12 and still got to fail time 9.85e8.
Make sure there is nothing else, and it is wierd that the major derivatives are behind and still the calculation encounters difficulties in the late process.

30.12.23

Hi, this is Itay from late - December here. Sorry for being out for quite long, Prague was nice and required my attention.
Looking at the last efforts I think the direction of stretching the model into a high vbar model is quite problematic, and it might be better to invest a little bit more time into relevant temperature profile, change it smartly (with regular vbar dbar) until a flip happens (easy success to fast failure) and then it will be feasible to search for the right vbar dbar for a case seperation.
I am really looking forward to the time this procedure will be automated and not a matter of art, but for a first step in the field I guess it is ok and indeed more complex than it seems with it being an already highly dimensional problem. 

31.12.23

Turns out increasing tolerances actually helps the solution to come, and the failures so far where not for population, but for numeric reasons.
So far did not manage to kill them all, will keep trying.

4.1.24
When I come to think of it, CCmoderateSearch5 is quite interesting. It indicates the living creatures on the planet do not care about the total dT over a billion years, but just about the derivative (rate) of change.
There , the maximal rate is 2.85 deg/Myr. Even when the dT=~60deg, it is manageable by a dbar=1e-7, vbar=1e-5.
The only way to kill them is to increase dT/dt relatively to vbar,dbar.
Let us examine the limit of change.

Newsflash:
Atol and Rtol effect the result, namely the cutoff (fail) time, which is never of the n=0 type, but rather a sharp (non differntiable function) density profile, on high values (order of magnitude 1e1 or 1e2). This is an old issue solved usally by small atol and rtol, but here these small values (1e-10) yeilds a failtime=0.
Whats interesting is for atol=0,rtol>0 there are 2 species at final frame, but for atol>0 ,rtol=0 there is only one surviving.
Error log is a bit more specific now, seperating the mean(n)<0 and the max(n)<1e-5 cases.

5.1.24
So now with the temperature col in dat it is clearer that the species that have significant in a certain patch follow very well the Tenv, it is just the numerical calculation that stopped at a certain point since the population did not undergo an extinction.
If not fixed today, TODO : Figure out why the global assignments from the error handling section are not working (fail time not updated to >0, for instance, outfile as well).

6.1.24
I will get to the globals fix shortly, now I am more interested in tolerances search since it became a non-scientific but rather technical obstacle which need to be addressed NOW.
Plus, for some reason not every run is saved, which is good now, but later needs to be checked.

A little summary:
For the rather difficult (to calculate!) profile of "Target2Dead?" (Oscillating and ascending eccentricity up to 1, with T ranges from -25 to +70 at Equator over 1Gyr)
All tolerances over 1e-7 lead to an exploding solution at 7.75e8 yr,
 and all tolerances under 1e-7 lead to pre-cutting the solution at earlier times, the smaller the tol the sooner our model stops, down to failtime=0 somewhere below tol=1e-12.
Tried varying the maxsteps around tol=1e-7, did not effect the failtime.
Maybe the profile itself is so extreme the population really dies? But there is no logic in sharp edges, this is a numerical runaway.
(when looking a few step before the last time stamp it is clear that there was no unusual change in climate, relatively, but the max population density is unprecedented throughout the entire simulation in an order of magnitude, and the peak is sharp)

There is a way to examine this issue. Increase the dbar, vbar to a bizzare value , and see if the system manages to get along the entire CC, or is it something in the derivatives of the T_kozai itself that won't allow convergence.
At least with the tol >1e-7 (i.e 1e-6) we have results/dat to look at.With the lower tol we don't even have those.
When increased vbar:1e-5->1e-2 , failed at exact same point.
Maybe the dbar was the surviving mechanism all along? Thats odd, but plausible, giving the vbar is sufficient.Will check as well.

Alright, nuliffied the vbar dbar differences between species.
The reason Species #2 is always winning , up to a final bizarre explosion, is the competition matrix a, which is quite hard to understand without further investigation of rhs_eval, but rho, the resource growth-tolerance tradeoff parameter, indicates a clear discrimination - rho is maximal for #2, then #4 comes second both in density and rho, and so on.

Alright, new conclusions- the explosion occurs on the most competitive species, only at rises to maximal tempratures, and becomes more and more peaky(sharp) and higher when the peak temprature oncreses from cycle to cycle.
It occurs also between 6.98e8 and 7.12e8, the peak rises and diminishes, and also between 6.27e8 and 6.35e8. In the earlier period peak reached n=40, in the latter crossed n=60, and the last one is around 7.77e8, which is the current fail time. I assume the peak was even higher at the last time (using v=1e-2, for v=1e-5 it was accordingly lower for all, but same principle)
This has to do with the Tdep model in my opinion, but the rather strange thing is that the peak goes strong where other species are pretty much extinct, temproraily, while in regions where others still live the densities of the dominant species are much more modest, forming a exponential rise towards the more empty regions.
It might not seem strange at all, but the question is who causes who. Was the absence of other species the one that allowed unproportioal thrive? Or the other way around? And if that is the case, why stop at a sub-equatorial region and not go all the way to the equator where no one else lives? It might be too hot there.
To answer this I must check the trait of species #2 at those periods, even though it is a fully responsive variable to the density in the patch, and vice-versa. But- maybe it will shed some light anyway.

From the lookup at T_kozai, at the time of the last peak cycle which was completed (circa 7.02e8) the maximal T at equator was 46 deg. At failtime T=50deg, and it was not the maximum of the relevant cycle. So the density runaway is Temprature dependent, when 50 degrees were unprecedented, maybe it can be killed synthetically by a very low vbar, but it will just set the initial conditions a little bit back and purspone the problem, maybe up to full 1Gyr but this is an unwanted effect since it is unreliable to thrive in harsher conditions.
It should be mentioned that when Tequator=46 or 50 no one was there, only for the earlier cycle, t= 6.27e8 to 6.35e8, there was a species #2 presence in the equator, around 39 degrees at most.
Throughout the simulation it seems adaptability is not something these creatures are hard dealing with, quite the opposite. Over 1 step , mathematically speaking, the living creatures change their entire nature quite in sync with the environment, in certain patches ofcourse, but still.
The thing is a weird mechanism which sharpens the density profile as Temperature goes higher up to 50.
Another way to check if this is an absolute T issue is just raise the entire profile by a couple of Celsius and see if the problem occurs earlier. The complement of lowering and later fail time is possible as well but since fail time is already quite late and intial T is quite cold, I prefer the first direction.

Just run it, the fail time did come earlier , but it was not at a peak or even an extermum at all , quite ambient T - 8 to 13 degrees, with mild change rate, less than 3deg/Myr. I am not sure what is going on, and so this is a great time to stop.

11.1.24
Hello there 25.11.23 note! I see a certain recreation of Naoz2016 figure 4 was made, but it is still interesting for me to reconstruct another example. I need also to make sure that I do not perhaps have a Sanity check which corresponds to Naoz.
Got return code -5 for all cases. For starters, I though figure 5 data is causing Inf values and hence the non-convergence , but when shifting to solar scale examples still the convergence problem sticks.
I am to track the version from November in github and see what is different.Tolerances are always first suspects but I changed them from 1e-2 to 1e-15 with no use.

12.1.24
After I thought I had good reconstructions of paper results, e.g Fig4 Blue from Naoz 2016, in the Fig 18 recreation the results where similiar but not identical. Then the dilema of comapring different equations systems arouse, and a simpler example was examined, namlely Fig3 Red, which also appears in Naoz 2013 as Fig4.
cos(i1) surpasses 1 (like in other cases in the past) and so code -5 occurs. This is normally occuring for inclinations out of range for the Kozai to work, but this time no angle at all managed to yield the correct solution, even when the profile was not cut and ode was succsesful.
I started to review my equations one by one. My suspicion forst layed on 2 places where I assumed Sin^2(itot) and not Sin(itot^2) where it was written as sinitot^2.
The alternative did not change at all Nan / fail time.
I will continue reviewing my equations, and if all correct, then there is not much left for me to do and move on. All verified ones are marked with "#V" on line-end in kozai_SanityRecreate.R.


15.01.23
An Error Was found and fixed during verfication! It was cos(w1) instead of the correct cos(2w1) at the (3-5cos(2w1)) term in dw1/dt.
In code, abrevaited as co1 and c2o1 accordingly.
Corrected also at sanityCheck and main kozai. Only at Recreate are the #V noted.
And it did not solve the cos>1 problem...what is going on?
It also did not effect the regular kozai much, which is positive from the side of not needing to re-run everything, but also odd from the aspect of not knowing the problem's origin.
It might be that comparing full octupole order with test particle wuatropole and the liking is just not going to work in some cases.
Then again, if I write all the TPQ equations and test them, that does not prove anything on the FO (full octupole) , just letting me use the TPQ which I am not sure is better than the most general set of equations.
I will try to look for FO results in other papers tommorow, and hopefully will find and avial to recreate one example with cutting evidence.


19.01.24
Sanity Check is Complete!
Turns out My Full Octopule Order yields slightly different results than the paper's Quadropole order, as demonstrated at the small window in Fig18 at Naoz2016.
Additionally I understood that to get Naoz results , i1 should be very close to itot,(i.e , itot=80, i1=79.9). Only then did things start to look familiar.
As for understanding the approximations, those are fairly simple - assymetric aprox is just e2 = 0. Test particle approx is m2--> 0.
Some configurations of the full set I have can get close to these limits' but not all of them. That is why Fig18 was my main proof for eqaution validity, in addition for a second review line-by-line.
The -5 return code error was avioded with big masses as I rewrote C3 calculation such that no elemnt will exceed e300. (Namely, changing m^8 to m^7 * m3 in the C3 calculation line).

Next step is take any kozai profile, and lower the genetics.
Even though I have a bountiful of kozai profiles, even the most simple of them all (2 cycles over 1e8 yrs) is failing and succeeding with very exagerated vbars (1e-3 vs 1e-5).
Tolerances helped a little to lower calc times (down to 3 min with 1e-8) but it still much longer than in the past. What is the reason I am not sure.
Plus, the deriative is the one failing the population, not the shape of the profile, and so it is not very interesting because there is no periodicity in the failure. Then again, maybe it was never so, but just a barrier deriative.

21-22.01.24
For comparison, the ecoevo_main.R from original directory is still running on 20 sec scale of time.
At first I thougth the difference is th timespan, so kozai also produced new prifles in 1e6 instead of 1e9 years.
Then deriatives were lowered since code -1 kept happening and no ode was solved yet in the old-new t profiles.
Notably, I dicovered that in nearly circular orbits, the oscilations are well sinusiadal (good point to interface with past runs) and not rising in eccentricity with offset. Jupiter mass\distance  controls the frequency and i1 controls the max e1 and T accordingly.
This way the "large" vs "small" adaptaiton time is quite straight forward constructed,see examples in kozai_parameters.

v2e-03_d1e-05idLargeAdaptTimeMildAlive13844_FAILED is a very intersting output.
First of all it represents runs that did not actually failed but only after tE, which is a phenomenon I still do not know the reason for and how to avoid it.
In peak tempratures the polar population is reduced to almost below threshold and survives the extreme. In later periods the distribution is already fairly sparse and so the middle and eqautorian populations can benefit from the climate change once established. In this output specifically the dbar was small enough so dispersion wont happen too fast, but will eventually succeed.
Still need to run the SmallAdaptTimeMild and see a big failure (real one), but since the timings arent solved for now, I will start to treat it as an actual complex computation, and so insert parrallelism to my cpp (using fork()), and see if runtimes get better.

22.01.24
As you might expect, fork theory is very intersting.
You have installed the dependencies for pasl library but yet to understand how to include it properly.
Using the standart unistd.h fork() will be nasty and unsafe, i would prefer thinking about how to integrate pasl or a similiar package succesfully.
In addition, after a "uname -a" command in terminal, my kernel is indeed supporting SMP (symmetric multiproccesing) which means each new process will be schedule by the kernel to (optionally) different CPU cores, but not mandatory.

23-24.01.24
So the efforts for threading are progressing.
Turns out the flag added as comments in the export above Rcpp  function does the work , but for some reason when all RcppThread is included in the cpp file, weird errors about ending paranthesis missing pop up, even though when not including it did not appear (and even without using the parralel code itself) so it is something about the c++11 inclusion? not sure.
When it is solved, it can be tested on the first loop, and then implemented on all loops in the script, starting from the inner ones creating L=20 threads at a time. Since our OS is SMP supporting, it is supposed to be shared over available proccesors even though I write threads and not procceses. 

```
