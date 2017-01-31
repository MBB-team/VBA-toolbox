---
title: "Bayesian mentalizing (k-ToM)"
---

{:toc}

When it comes to interpreting others' behaviour, we almost irrepressibly engage in the attribution of mental states (beliefs, emotions...). Such "[mentalizing](https://en.wikipedia.org/wiki/Mentalization)" or [Theory of Mind](https://en.wikipedia.org/wiki/Theory_of_mind) (ToM) can become very sophisticated, eventually endowing us with highly adaptive skills such as convincing, teaching or deceiving. Here, sophistication can be captured in terms of the depth of our [recursive beliefs](https://en.wikipedia.org/wiki/Hierarchy_of_beliefs), as in "I think that you think that I think...". In [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992), we show that such sophisticated recursive beliefs subtend learning in the context of social interaction (e.g., repeated [dyadic](https://en.wikipedia.org/wiki/Dyad_(sociology)) games). Critical here is the notion that update rules that account for recursive beliefs of this sort (hereafter: *k-ToM* models) induce behavioural responses that are very different from typical [trial-and-error learning)(https://en.wikipedia.org/wiki/Trial_and_error) schemes.

The next sections explain how to use VBA's k-ToM Bayesian mentalizing models.


## Summary principles of k-ToM mentalizing models

First, recall that, in its simplest form, a [game](https://en.wikipedia.org/wiki/Game_theory) is formally defined in terms of a utility table $$U(a,b)$$, which yields the payoff player $$a$$ gets when making his decision, depending on what player $$b$$ decides. Incentives can be arbitrarily chosen to capture different forms of social exchanges or transactions, such as competitive or cooperative interactions. In what follows, we will only deal with 2x2 games, i.e. games with two players, both of which have only two alternative options.

Second, according to Bayesian [decision theory](https://en.wikipedia.org/wiki/Decision_theory), agents aim at maximising expected payoff $$E\left[U(a)\right]$$, where the expectation is defined in relation to the agent's uncertain prediction $$P(b)$$ about his opponent's next move. The repeated observation of his opponent's behaviour gives the agent the opportunity to learn this prediction. Theory of Mind comes into play when agents consider that the opponent's behavioural tendency is motivated by his hidden beliefs and desires. This is what k-ToM learners do (when $$k>0$$):

- **0-ToM** tracks her opponent's probability of picking any alternative option, as it changes over time or trials. 0-ToM does not mentalize: she essentially operates a nonlinear [Kalman filter](http://en.wikipedia.org/wiki/Kalman_filter) on her opponent's moves.
- **1-ToM** assumes she is facing a 0-ToM agent. She mentalizes, i.e. (in Bayesian terms) she holds an uncertain belief about her opponent's prediction about her own actions. To update this belief, she has to learn how her opponent learns, by estimating (online) hidden parameters of 0-ToM. This makes 1-ToM a meta-Bayesian agent ([Daunizeau et al. 2010](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0015554)), i.e. a Bayesian observer of a Bayesian agent. 
- **2-ToM** assumes she is facing either a 0-ToM or a 1-ToM, and learns her opponent's hidden parameters, as well as her opponent's ToM sophistication level. The difference between 2-ToM and 1-ToM is that 2-ToM can deal with mentalizing agents. 
- More generally, **k-ToM** assumes she is facing a k'-ToM agent with $$k>k'$$, and learns her opponent's hidden parameters and sophistication. 

Here, the recursion depth ($$k$$) induces distinct ToM sophistication levels, which differ in how they update their (recursive) beliefs and the ensuing trial-by-trial predictions about their opponent's next move. Note that one can derive the learning rule of any k-ToM agent recursively, from the learning rule of 0-ToM.


## k-ToM mentalizing models in VBA

In VBA, there is a pair of generic evolution and observation functions for k-ToM agents, namely: `f_kToM.m` and `g_kToM.m`. These can be used in conjunction with the appropriate inputs (previous players' choices) to fit the series of observed choices in a 2x2 game, whose utility table is given. The script `demo_recur.m` exemplifies the simulation and inversion of such k-ToM learner (using a competitive 2X2 game, namely: "[hide-and-seek](https://en.wikipedia.org/wiki/Hide-and-seek)").

Constraints and comments for applying the existing VBA's k-ToM models are as follows:

- the game has to be a 2X2 game.
- the structure of the utillity table of k-ToM's opponent is fixed, but does not need to be identical to k-ToM's.
- states indexing has to be set according to a standard. In the demo, this is done using the function `defIndlev.m`.
- some additional behavioural forces (e.g., [perseveration](https://en.wikipedia.org/wiki/Perseveration) and/or directed exploration) may be easily inserted in k-ToM's observation function.

Once the model has been inverted, one may be willing to recover and interpret the estimated hidden states. The indexing is a bit tricky here, essentially because it depends upon k (`f_kToM.m` is recursive):

- k-1 first hidden states : sufficient statistics of the agent's posterior belief about her opponent's sophistication level. More precisely, passing these states through `sigm.m` yields the posterior probability that the opponent's sophistication is k'=0,1,...,k-2 (by construction, the probability for k'=k-1 is `1-sum(sigm(posterior.muX(1:k-2)))`).
- the next hidden states are the represented hidden states of the agent's vitual opponents, for all admissible sophistication level. The corresponding indices are stored in the structure `options.inG.indlevel`. For example, `options.inG.indlev(2).X` contains the indices of the represented hidden states of a virtual 2-ToM opponent (only posible if the agent is himself at least 3-ToM). Note that the number of these indices is increasing with the sophistication level (because of the recursive nature of k-ToM beliefs).
-	then `f` (where `sigm(f)` is the probability that the agent's virtual opponent will pick the first option) and its partial derivatives w.r.t. learning parameters (typically volatility and behavioural temperature). Note: we keep track of the latter only to limit the computational burden of k-ToM learning. The corresponding indices are stored in `options.inG.indlev(2).f` and `options.inG.indlev(2).df` (for a 2-ToM virtual opponent).
- Finally, the [sufficient statistics](https://en.wikipedia.org/wiki/Sufficient_statistic) of the agent's posterior belief about his virtual opponents' learning parameters (mean and variance for each learning parameter: [M1,V1,M2,V2,...]). The corresponding indices are stored in `options.inG.indlev(2).Par` (for a 2-ToM virtual opponent).


Below are optional inputs that can be changed to adapt the model to one's specific experiment:

- `options.InF.game`: 2x2x2 payoff table (`game(:,:,1)` is the first player's payoff, whereas `game(:,:,2)` is the second player's payoff).
- `options.inF.player`: this specifies which 2x2 entry of `options.InF.game` refers to the k-ToM agent (1 or 2).
- `options.inF.lev`: this is the player's ToM sophistication level ($$k$$).
- `options.inF.dummyPar`: this is 3x1 binary vector that specifies which hidden parameter of her opponent k-ToM assumes to vary across trials (volatility, temperature and bias). By default, it is set to `[1;0;0]`, which means that only k-ToM's opponent's prior volatility is assumed to change across trials.


All relevant mathematical details are described in [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992).







