---
title: "Bayesian sequence learning"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


Obviously, Bayesian learning is not limited to simple cue-action-outcome associations.
In particular, one may rely on Bayes' rule to learn, on-line, the structure of sequences, i.e. patterns of outcomes that are specified in terms of transition probabilities. This is what "Bayesian Sequence Learners" (BSL) do.

In VBA, BSL is simply tracking the log-odds of $$P\left(u_t=1\mid u_{t-K}\right)$$, where $$u$$ is a binary outcome. This variable is updated according to a Laplace-Kalman filter, yielding 2 sufficient statistics (mean and avariance) per combination of past outcome. BSL can learn sequences of arbitrary depth ($$K$$). For example:

- if $$K=1$$, then BSL tracks 2 probabilities, namely: $$P\left(u_t=1\mid u_{t-1}=1\right)$$ and $$P\left(u_t=1\mid u_{t-1}=0\right)$$. In this case, BSL needs to know about the previous outcome $$u_{t-1}$$.
- if $$K=2$$, then BSL tracks 4 probabilities, namely: $$P\left(u_t=1\mid u_{t-1}=1,u_{t-2}=1\right)$$, $$P\left(u_t=1\mid u_{t-1}=0,u_{t-2}=1\right)$$, $$P\left(u_t=1\mid u_{t-1}=1,u_{t-2}=0\right)$$ and $$P\left(u_t=1\mid u_{t-1}=0,u_{t-2}=0\right)$$. In this case, BSL needs to know about the two previous outcomes $$u_{t-1}$$ and $$u_{t-2}$$.
- etc.

More generally, BSL tracks $$2^K$$ probabilities. In this scheme, the only evolution parameter is BSL's prior volatity about the log-odds.

Note: unsampled sequences will eventually be "forgotten", since the prediction step in the Laplace-Kalman update will dilute any previously sampled evidence.

The script `demo_BSL.m` can be used to simulate and invert k-BSL algorithms. Note that, in addition evolution (`f_BSL.m`) and observation (`g_BSL.m`) functions, VBA provides a simple tool for displaying BSL's evolving belief about transition probabilities (`unwrapKBSL.m`). This is exemplified below:

![]({{ site.baseurl }}/images/wiki/unwrapBSL.bmp)

> Example of evolving belief of a 1-BSL learner (unknowingly tracking random noise). **Top**: 1-BSL's trial-by-trial bet about the upcoming outcome. **Bottom**: 1-BSL's trial-by-trial estimate of the outcome's transition probabilities $$P\left(u_t=1\mid u_{t-1}\right)$$ (given each possible previous outcome).
