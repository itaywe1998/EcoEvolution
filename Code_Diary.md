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
```Tmax <- 15.0 # initial mean temperature at equator
Tmin <- Tmax-30 # initial mean temperature at poles
Cmax <- 30 # projected temperature increase at poles
Cmin <- 12 # projected temperature increase at equator```
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
Lesson learned is not trying to print each time such a small function is called. Runtime skyrocketed for this prints alone.
Plus, there is no way to compare the amount of times of function was called, but by the logic of the rhs_eval function structure, it is soppused to be indentical call number.
So for now the check is enough, and the run time extension is probably not 
due to sin() itself, but the complexity the frequent CC is loading on the ODE solver.


```
