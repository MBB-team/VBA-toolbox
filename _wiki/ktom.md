---
title: "Bayesian mentalizing (k-ToM)"
---

{:toc}

When it comes to interpreting others' behaviour, we almost irrepressibly engage in the attribution of mental states (beliefs, emotions...). Such "mentalizing" or Theory of Mind (ToM) can become very sophisticated, eventually endowing us with highly adaptive skills such as convincing, teaching or deceiving. Here, sophistication can be captured in terms of the depth of our recursive beliefs, as in "I think that you think that I think...". In [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992), we show that such sophisticated recursive beliefs subtend learning in the context of social interaction (e.g., repeated dyadic games). Critical here is the notion that belief update rules that account for recursive beliefs of this sort (hereafter: k-ToM models) induce behavioural responses that are very different from typical trial-and-error based learning schemes. 
This section summarizes the prequisites for using k-ToM models in dyadic games.

First, recall that, in its simplest form, a game is defined in terms of a utility table U(a,b), which yields the payoff player a gets when making his decision, depending on what player b decides. Incentives can be arbitrarily chosen to capture different forms of social exchanges or transactions, such as competitive or cooperative interactions. In what follows, we will only deal with 2x2 games, i.e. games with two players, both of which have only two alternative options.

Second, according to Bayesian decision theory, agents aim at maximising expected payoff E[U(a)], where the expectation is defined in relation to the agent's uncertain prediction p(b) about his opponent's next move. The repeated observation of his opponent's behaviour gives the agent the opportunity to learn this prediction. Theory of Mind comes into play when agents consider that the opponent's behavioural tendency is motivated by his hidden beliefs and desires. This is what k-ToM learners do (when k>0):

- 0-ToM tracks her opponent's probability of picking any alternative option, as it changes over time or trials. 0-ToM does not mentalize: she essentially operates a nonlinear Kalman filter on her opponent's moves.
- 1-ToM assumes she is facing a 0-ToM agent. She mentalizes, i.e. (in Bayesian terms) she holds an uncertain belief about her opponent's prediction about her own actions. To update this belief, she has to learn how her opponent learns, by estimating (online) hidden parameters of 0-ToM. This makes 1-ToM a meta-Bayesian agent (Daunizeau et al. 2010), i.e. a Bayesian observer of a Bayesian agent. 
- 2-ToM assumes she is facing either a 0-ToM or a 1-ToM, and learns her opponent's hidden parameters, as well as her opponent's ToM sophistication level. The difference between 2-ToM and 1-ToM is that 2-ToM can deal with agents that can mentilize. 
- More generally, k-ToM assumes she is facing a l-ToM agent with k>l, and learns her opponent's hidden parameters and sophistication. 

Here, the recursion depth (k) induces distinct ToM sophistication levels, which differ in how they update their beliefs and the ensuing prediction about their opponent's next move. Note that one can derive the learning rule of any k-ToM agent recursively, from the learning rule of 0-ToM.


In VBA, there is a pair of generic evolution and observation functions for k-ToM agents, namely: `f_kToM.m` and `g_kToM.m`. These can be used in conjunction with the appropriate inputs (previous players' choices) to fit the series of observed choices in a 2x2 game, whose utility table is given. The script `demo_recur.m` exemplifies the simulation and inversion of such k-ToM learner.

Below are optional inputs that can be changed to adapt the model to one's specific experiment:

- options.InF.game: 2x2x2 payoff table (game(:,:,1) is the first player's payoff, whereas game(:,:,2) is the second player's payoff). For exmaple (as in the demo script), setting `game=cat(3,[1,0;0,1],[0,1;1,0])` induces a competitive ("hide-and-seek") game, where the first player is the seeker and the second is the hider.
- options.inF.player: 1 or 2. In the above example, this flag effectively switches the role (hider or seeker) of the player, whose behaviour we are trying to fit.
- options.inF.lev: this is the player's ToM sophistication level (k).
- options.inF.dummyPar: this is 3x1 binary vector, which states which hidden parameter k-ToM assumes to vary across trials (volatility, temperature and bias). By default, it is set to [1;0;0], which means that only k-ToM's opponent's prior volatility is assumed to change across trials.

Note that some states indexing has to be set according to: `options.inG.indlevel=defIndlev(level, NtotPar)`, where `level` is the player's ToM sophistication level, and `NtotPar` is the total number of hidden evolution and observation parameters for k-ToM players (here, `NtotPar=3`: volatility, temperature and bias).

All relevant mathematical details are described in [Devaine et al. (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003992).







