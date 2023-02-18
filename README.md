# special-potato
BioInformatics II homework sets

The project now is about creating a hidden markov model....

A sequence of length 1200 was generated with a hidden Markov model that allowed some regions to be GC-rich and some regions to be  AT-rich.  The model used to generate the data was very similar to the one that we studied in class... Position 1 in the sequence was forced to be in a AT-rich region. Thereafter, the probability that position i+1 was in a GC-rich region given that position i was in a GC-rich region was set to 0.95. 

The probability that position i+1 was in an AT-rich region given that position i was in an AT-rich region was set to 0.98.

In GC-rich regions, the probabilities of G and C occupying sites were each 0.3 whereas the probabilities of A and T in these GC-rich regions were each 0.2.

In AT-rich regions, the probabilities of A and T occupying sites were each 0.3 whereas the probabilities of G and C in these AT-rich regions were each 0.2.

The sequence of length 1200 that was generated according to this model can be found as 'hwk2_sequence.txt'

1.  Implement the forward algorithm to calculate the logarithm of the probability of observing this sequence given the model.  When calculating this  log-likelihood, do it for the true values of the parameters that are listed above. What log-likelihood value results?

2.  Implement the forward algorithm again to analyze the same data set.  But, this time use the true values of the parameters except 0.5 should be the probability that position i+1 was in a GC-rich region given that position i was in a GC-rich region and 0.8 should be the probability that position i+1 was in an AT-rich region given that position i was in an AT-rich region.  What value of the log-likelihood results?

3. Implement the forward algorithm again to analyze the same data set.  But, this time use the true values of the parameters except 0.51 should be the probability that position i+1 was in a GC-rich region given that position i was in a GC-rich region and 0.51 should be the probability that position i+1 was in an AT-rich region given that position i was in an AT-rich region.  What value of the log-likelihood results?
 

4.  For the parameters values listed in Question 1, implement the Viterbi algorithm to analyze the simulated sequence data.  The implementation should report the most probable path of AT-rich and GC-rich states through the hidden Markov model.  It should also report the logarithm of the probability of generating the sequence data with this specific path.

5.  Answer Question 4 again except now use the parameter values that were used in Question 2.  In a few sentences, try to explain why the reconstructed path from this question differs from the path reconstructed for Question 4.

6.  Answer Question 4 again except now use the parameter values that were used in Question 3.  In a few sentences, try to explain why the reconstructed path from this question differs from the path reconstructed for Question 4.

This homework assignment should be returned to me via email before class on Friday February 24.  As with the last assignment, please email to me the computer code that you write,  the output that answers the questions, and the commands that you used to compile and run your program.  When doing this assignment, you can assume  that the length of the sequence was known a priori to be 1200 sites and that the first position in the sequence was known a priori to be in an AT-rich region.

(hwk instuctions were found at: http://statgen.ncsu.edu/thorne/bioinf2hwk2.html)
