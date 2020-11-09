---
title: "Bayesian mentalizing (k-ToM)"
---

{:toc}

When it comes to interpreting others' behaviour, we almost irrepressibly engage in the attribution of mental states (beliefs, emotions...). Such "[mentalizing](https://en.wikipedia.org/wiki/Mentalization)" or [Theory of Mind](https://en.wikipedia.org/wiki/Theory_of_mind) (ToM) can become very sophisticated, eventually endowing us with highly adaptive skills such as convincing, teaching or deceiving. Here, sophistication can be captured in terms of the depth of our [recursive beliefs](https://en.wikipedia.org/wiki/Hierarchy_of_beliefs), as in "I think that you think that I think...". In [Devaine et al. (2014b)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992), we show that such sophisticated recursive beliefs subtend learning in the context of social interactions (e.g., repeated [dyadic](https://en.wikipedia.org/wiki/Dyad_(sociology)) games). Critical here is the notion that update rules that account for recursive beliefs of this sort (hereafter: *k-ToM* models) induce behavioural responses that are very different from typical [trial-and-error learning](https://en.wikipedia.org/wiki/Trial_and_error) schemes.

The next sections explain how to use VBA's k-ToM Bayesian mentalizing models.


## Summary principles of k-ToM mentalizing models

First, recall that, in its simplest form, a [game](https://en.wikipedia.org/wiki/Game_theory) is formally defined in terms of a utility table $$U(a,b)$$, which yields the payoff player $$a$$ gets when making his decision, depending on what player $$b$$ decides. Incentives can be arbitrarily chosen to capture different forms of social exchanges or transactions, such as competitive or cooperative interactions. In what follows, we will only deal with (repeated) 2x2 games, i.e. games with two players, both of which have only two alternative options.

Second, according to Bayesian [decision theory](https://en.wikipedia.org/wiki/Decision_theory), agents aim at maximising expected payoff $$E\left[U(a)\right]$$, where the expectation is defined in relation to the agent's uncertain prediction $$P(b)$$ about his opponent's next move. The repeated observation of his opponent's behaviour gives the agent the opportunity to learn this prediction. Theory of Mind comes into play when agents consider that the opponent's behavioural tendency is motivated by his hidden beliefs and desires. This is what k-ToM learners do (when $$k>0$$):

- **0-ToM** tracks her opponent's probability of picking any alternative option, as it changes over time or trials. 0-ToM does not mentalize: her learning rule essentially operates a (nonlinear) [Kalman filter](http://en.wikipedia.org/wiki/Kalman_filter) on her opponent's moves. This is a slightly more sophisticated variant of so-called "[fictitious play](https://en.wikipedia.org/wiki/Fictitious_play)" strategy, which best responds to the opponent's frequency of play.
- **1-ToM** assumes she is facing a 0-ToM agent. She mentalizes, i.e. (in Bayesian terms) she holds an uncertain belief about her opponent's prediction about her own actions. To update this belief, she has to learn how her opponent learns, by estimating (online) hidden parameters of 0-ToM. This makes 1-ToM a meta-Bayesian agent ([Daunizeau et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015554)), i.e. a Bayesian observer of a Bayesian agent. 
- **2-ToM** mentalizes also, and assumes she is facing either a 0-ToM or a 1-ToM, and learns her opponent's hidden parameters, as well as her opponent's ToM sophistication level. The difference between 2-ToM and 1-ToM is that 2-ToM can deal with mentalizing agents. 
- More generally, **k-ToM** mentalizes and assumes she is facing a k'-ToM agent with $$k>k'$$, and learns her opponent's hidden parameters and sophistication. 

Here, the recursion depth ($$k$$) induces distinct ToM sophistication levels, which differ in how they update their (recursive) beliefs and the ensuing trial-by-trial predictions about their opponent's next move. Note that one can derive the learning rule of any k-ToM agent recursively, from the learning rule of 0-ToM.

> The mathematical details of k-ToM learning models are described in [Devaine et al. (2014a)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0087619) and [Devaine et al. (2014b)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992). Empirical applications of this model to ethology and computational psychiatry can be found in [Devaine et al. (2017)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005833) and [Devaine et al. (2020)](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1007700).


## k-ToM mentalizing models in VBA

In VBA, there is a pair of generic evolution and observation functions for k-ToM agents, namely: `f_kToM.m` and `g_kToM.m`. These can be used in conjunction with the appropriate inputs (previous players' choices) to fit the series of observed choices in a 2x2 game. Constraints and comments for applying the existing VBA's k-ToM models are as follows:

- the game has to be a 2X2 game (2 agents, 2 actions)
- the structure of the utillity table of k-ToM's opponent is fixed, but does not need to be identical to k-ToM's.


### Running a k-ToM model inversion
Setting k-ToM models may at first seem tedious, because it requires specific optional arguments. But in fact, VBA provides a simple function that generates simulation and/or inversion input arguments, namely: `prepare_kToM.m`. Below is an example piece of code that would run a k-ToM model inversion in VBA (with $$k=2$$):
```matlab
% prepare specific dim and options structures for k-ToM model inversion
K = 2; % depth of k-ToM's recursive beliefs
payoffTable = cat(3,[1,0;0,1],[0,1;1,0]); % game payoff matrix (here: hide-and-seek)
role = 1;  % subject's 'role' (here: 1=seeker, 2=hider)
[options,dim] = prepare_kToM(K,payoffTable,role,0);
% run inversion
f_fname = @f_kToM; % k-ToM model evolution function
g_fname = @g_kToM; % k-ToM model observation function
y = SubjectChoices; % sequence of choice data 
u = [zeros(2,1),[OpponentChoices(1:end-1);SubjectChoices(1:end-1)]]; % sequence of players' actions (at the previous trial)
options.skipf = [1,zeros(1,size(y,2)-1)];  % skip 1st trial (no learning until trial #2)
options.binomial = 1; % inform VBA about binomial data
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
```
Here, `SubjectChoices` and `OpponentChoices` are the binary actions of the participant and of his/her "opponent", respectively.

Note that:

- `SubjectChoices` enters as data to be fitted (`y`), but also as an input to the model (second row of `u`). This is because k-ToM's learning rule requires the knowledge of her own action at the previous trial...
- the incentive rule of the game is completely specified by the payoff table (here, the game is "[hide-and-seek](https://en.wikipedia.org/wiki/Hide-and-seek)"). In our example, `payoffTable` is a 3D-array, which is such that `payoffTable(:,:,role)` (resp., `payoffTable(:,:,3-role)`) provides the payoff for the participant (resp., for his/her opponent) for all possible action pairs (2x2 possibilities) in the game, where `role` is set appropriately (1="seeker", 2="hider").
- The states' indexing of the k-ToM model is automatically set by `prepare_kToM.m` (see below).



### Eyeballing inversion results
Once the model has been inverted, one may be willing to recover and interpret the estimated hidden states. Note that `prepare_kToM.m` stores the states' indexing in the structure `options.inG.indlev`. It's a bit tricky because both the types and the number of states depend upon $$k$$ (`f_kToM.m` is recursive). Nevertheless, the hidden states of a k-ToM agent can be simply eyeballed using the function `unwrapKTOM.m`:
```matlab
hf = unwrapKTOM(posterior.muX,options.inG);
```
where `posterior` is the output of the above VBA inversion of the k-ToM model

An example graphical output is appended below:

![]({{ site.baseurl }}/images/wiki/unwrapkTom.bmp)

> Example of evolving belief of a 2-ToM learner (unknowingly playing against random noise). **Top**: 2-ToM's belief about her opponent's sophistication (here, either 0-ToM or 1-ToM). **Middle**: 2-ToM's trial-by-trial prediction about her opponent's next move, for each possible opponent's sophistication. NB: 2-ToM's overall prediction is a weighted average of these two conditional predictions, where the weights are 2-ToM's belief about her opponent's sophistication. **Bottom**: 2-ToM's trial-by-trial estimate of her opponent's unknown parameters (learning rate, behavioural temperature and bias), for each possible opponent's sophistication (`unwrapKTOM.m` proposes to toggle between each possible virtual opponent type, here: for a virtual 1-ToM opponent).

Editing `unwrapKTOM.m` enables users to extract any internal state (or a transformation thereof) of k-ToM learners. Nevertheless, for the sake of transparency, we provide below a summary of k-ToM's states' indexing:

- k-1 first hidden states : sufficient statistics (truncated log-odds) of the agent's posterior belief about her opponent's sophistication level. More precisely, passing these states through `sigm.m` yields the posterior probability that the opponent's sophistication is k'=0,1,...,k-2 (by construction, the probability for k'=k-1 is `1-sum(sigm(posterior.muX(1:k-2,:)))`). NB: this is only for 2-ToM learners and above (because 0-ToM and 1-ToM do not estimate their opponent's sophistication).
- the next hidden states are the represented hidden states of the agent's virtual opponents, for all admissible sophistication levels (here: k'=0 or 1). The corresponding indices are stored in the structure `options.inG.indlevel`. For example, `options.inG.indlev(2).X` contains the indices of the represented hidden states of a virtual 2-ToM opponent (only posible if the agent is himself at least 3-ToM). Note that the number of these indices is increasing with the sophistication level (because of the recursive nature of k-ToM beliefs).
-	then `f` (where `sigm(f)` is the probability that the agent's virtual opponent will pick the first option) and its partial derivatives w.r.t. learning parameters (typically volatility and behavioural temperature). Note: we keep track of the latter only to limit the computational burden of k-ToM learning. The corresponding indices are stored in `options.inG.indlev(2).f` and `options.inG.indlev(2).df` (for a 2-ToM virtual opponent).
- Finally, the [sufficient statistics](https://en.wikipedia.org/wiki/Sufficient_statistic) of the agent's posterior belief about his virtual opponents' learning parameters (mean and variance for each learning parameter: [M1,V1,M2,V2,...]). The corresponding indices are stored in `options.inG.indlev(2).Par` (for a 2-ToM virtual opponent).

We acknowledge that the states' indexing of k-ToM models might seem excessively complex at first sight. But, in our experience, this standard has the virtue of enabling one to extract very rich information from model inversions...


### Adapting k-ToM models to one's experimental needs
Many optional features of k-ToM learners can be changed to adapt the model to one's specific experiment. Some of these features are in fact inputs to the function `prepare_kToM.m`:

- the recursion depth $$k$$ of k-ToM learners
- the game's payoff table $$U(a,b)$$
- the possibility to include "belief dilution" effects (e.g. forgetting effects)

Other features have to be reset by editing the function `prepare_kToM.m`, or by directly modifying the appropriate optional inputs stored in `options.inF` and `options.inG`. For example, one can reset the field `options.inF.dummyPar`, which is a 3x1 binary vector that specifies which hidden parameter(s) of her opponent k-ToM assumes to vary across trials (volatility, temperature and bias). By default, it is set to `[1;0;0]`, which means that k-ToM assumes that her opponent's volatility is the only thing that changes across trials. The impact of setting any of these entries to 1 is that the induced learning rule includes some form of "forgetting effect" which results from the dilution of belief precision from one trial to the next.


> Some additional behavioural forces (e.g., [perseveration](https://en.wikipedia.org/wiki/Perseveration) and/or directed [exploration](https://en.wikipedia.org/wiki/Exploration)) may be easily inserted in k-ToM's observation function.


The script `demo_recur.m` exemplifies the simulation and inversion of a k-ToM learner (using a competitive 2X2 game, namely: "hide-and-seek"), where some of these features are edited...

> Practical experience with k-ToM models show that the statistical power, in terms of both estimating unknown model parameters and/or comparing different learning models, critically depends upon the experimental design. In [Devaine et al. (2014b)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992), we had participants play "hide-and-seek" against on-line k-ToM artificial learners. This was shown to yield much better model discriminability than, e.g., if participants had played cooperative and/or coordination games (e.g., "[battle-of-the-sexes](https://en.wikipedia.org/wiki/Battle_of_the_sexes_(game_theory))"). We encourage users of k-ToM models to perform some form of a [confusion analysis](https://en.wikipedia.org/wiki/Confusion_matrix) to chack that their intended experimental design eventully provide reasonable statistical power.







