![hidden logo](logo.png)
 
*0.1.0*
 
`hidden` is an `R` package for simulation of hidden Markov models (HMMs).
I wrote this in a hurry because I was confused about HMMs and needed a way to understand them properly.

It implements the filtering, prediction and smoothing algorithms via recursion and the Viterbi algorithm via dynamic programming.
I think it might be useful mostly for learning purposes.
All code is written in R and it uses some simple matrix algebra to get things done.

## Install

Use `devtools` to install this package.
In `R`, run

``devtools::install_github('davnovak/hidden')``

Then, load the package by entering

``library(hidden)``

## Run

Let's take the classic example: you're a hard-working PhD student confined to a basement office, on a diet of coffee and instant ramen noodles.
You have no idea what the weather outside is.
Except that your advisor arrives to work everyday, either with an umbrella or without (an observable state variable).
That gives you a hint as to whether it is raining outside or not (a hidden state variable).

We suppose that rain can be modelled as a 1st-order Markov process (hidden state on day *t* is dependent only on the hidden state on day *t-1*) and that the sensory assumption hold true (observable state is entirely determined by the hidden state at the same time and conditionally independent of the hidden state and the observable state at any other time).

We use the transition model (a matrix) to tell us how probable rain is on day *t* given the probability distribution of the hidden state variable (of *'Rain'* and *'No Rain'*) at *t-1*.
We use the emission model (also a matrix) to tell us how probable it is that the advisor bring an umbrella on day *t* given the probability distribution of the hidden state variable on the same day (*t*).
In addition to this we need an initial probability distribution of the hidden state variable and we are good to go.

Let's start by filling in the transition model and emission model matrices, as well as a column vector for the initial probability distribution.
In `hidden`, we create an S3 object (similar to an instantiation of a class in C++, Java, Python, etc.) of the class `HiddenMarkovModel`.

```
Tr <- Matrix(c(.7, .3, .3, .7), nrow = 2) # transition model: if it's raining, there's a 70-percent chance of the rain persisting;
                                          # if not, a 30-percent chance of rain the next day
Em <- Matrix(c(.9, .2, .1, .8), nrow = 2) # emission model: if it's raining; there's a 9-to-2 chance we'll see and umbrella;
                                          # if it isn't we only have a 1-to-8 chance
X0 <- Matrix(c(.5, .5), nrow = 2)         # initial probability distribution: fifty-fifty whether it rains or not

H <- HiddenMarkovModel(Tr, Em, X0,
                       hidden_states = c('Rain', 'No rain'),
                       observables   = c('Umbrella', 'No umbrella'),
                       name          = "My Umbrella Conundrum")
```

The next thing we do is supply the model with some observations.

```
SetObservations(H, c('Umbrella', 'Umbrella', 'No umbrella', 'No umbrella', 'No umbrella', 'No umbrella', 'Umbrella', 'No umbrella'))
```

We can then compute the posterior probability distributions of it raining on day *t* given

* the evidence up to day *t* or

* the evidence up to day *t* and futher (looking back and incorporating what we have observed since.

In the first case, we want to run an algorithm called *filtering*.
In the second case, we want to apply *smoothing*.
In addition, some observations can be set as `NA`, meaning they're in the future and we don't know yet.
In that case, we can use filtering except for having fewer data points to do *prediction*.

We can also find the most likely sequence of hidden state variable values using the *Viterbi algorithm*.

The package `hidden` can do all of these things.
To find more, just pull up the documentation for those algorithms that I just listed.

```
?Filter.HiddenMarkovModel
?Smooth.HiddenMarkovModel
?Viterbi.HiddenMarkovModel
```

Additionally use `print(H)` to view the internals of your model and `PrintStats(H)` to see the results of what you have already computed.

Have fun.



