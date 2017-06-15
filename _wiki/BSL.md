---
title: "Bayesian sequence learning"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}


Obviously, Bayesian learning is not limited to simple cue-action-outcome associations.
In particular, one may rely on Bayes' rule to learn, on-line, the structure of sequences, i.e. patterns of outcomes that are specified in terms of transition probabilities. This is what "Bayesian Sequence Learners" (BSL) do.

In VBA, BSL is simply tracking the log-odds of $$P(u_t=1|u_{t-1})$$, where $$u$$ is a binary outcome. This variable is updated according to a Laplace-Kalman filter, yielding 2 sufficient statistics (m and V) per combination of past outcome. BSL can learn sequences of arbitrary depth ($$K$$). For example, if $$K=1$$, then BSL tracks 2 probabilities, namely: $$P(u_t=1|u_{t-1}=1)$$ and $$P(u_t=1|u_{t-1}=0)$$. More generally, BSL tracks $$2^K$$ probabilities. In this scheme, the only evolution parameter is BSL's prior volatity about the log-odds.

Note: unsampled sequences will eventually be "forgotten", since the prediction step in the Laplace-Kalman update will dilute any previously sampled evidence.
