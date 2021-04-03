#!/usr/bin/env python

"""prop.py: 
    Builds dictionary that can be used to 
    1. Create exact sampling distribution of the proportion.
    2. Calculate the expected proportion, 
    explanation below
"""



#Note: this is a more efficient numerical routine to calculate the expected proportion
#relative to what appeared in the previous working paper (based on the joint distribution of runs)
#the code leverages a recursion idea suggested by Michael Wiener, and extends it to the general case

#Explanation:
#Consider the proportion of hits after 1 hit
#Each sequence has a tuple (#miss after hit, #hit after hit) 
#Because tuples are immutable, they can serve as keys in a dictionary.

#With 3 shots we have the following possibilities.
#sequence, dictionary = {key=(#miss after hit, #hit after hit): value=# sequences}
# 000       {(0,0):1}
# 001       {(0,0):1}
# 010       {(1,0):1}
# 100       {(1,0):1}
# 110       {(1,1):1}
# 101       {(1,0):1}
# 011       {(0,1):1}
# 111       {(0,2):1}

# we can add these dictionaries together to get:
#     B[0][3]= {(0,0):2,(1,0):3,(1,1):1, (0,1):2}
#were for B[r][c], 
# r=# of consecutive hits to the left, and c
# c=# of remaining shots to the right

# With this dictionary, it is trivial, and quick to:
# 1. Compute expected proportion (function defined below) 
# 2. Plot a histogram of its values (i.e. the sampling distribution)

# The trick then to define a recursion to build this dictionary
# The function outcome_and_frequency_dictionary below does this by
# initializing a few dictionaries, and then building on them,
# for each additional shot that is added, hit or miss.




#Use Counter dictionary subclass
#Reason: when combining 2 dicts, items with matching keys add their values.
from collections import Counter 

#With the dictionary dict_on_a_streak consisting of sequences
#in which adding the next outcome will create a hit after a streak of hits
#or a miss after a streaks of hits, then we need to increase the count of each

# Helper code for updating dictionary keys and values
def update_count_and_probability(hit_streak,hit,prob_hit,dict_on_a_streak):
    p=prob_hit
    q=1-p
    updated_dict=Counter({})
    for key, value in dict_on_a_streak.items():
        if hit_streak:
            return Counter({(key[0] + 1-hit, key[1] + hit):value*p**hit*q**(1-hit) \
                for key, value in dict_on_a_streak.items()})
        else:
            return Counter({(key[0], key[1]):value*p**hit*q**(1-hit) \
                for key, value in dict_on_a_streak.items()})


    
# Build the Matrix of dictionaries (Counter subclass). Defined recursively.
# For each dictionary the items correspond to a unique sequence "type" defined by:
# dict key=(#misses after k hits, #hits after k hits)
# dict value = probability. 
def outcome_and_frequency_dictionary(number_of_shots,streak_length,probability_of_hit):
   
    p=probability_of_hit
    q=1-p
    k = streak_length
    n = number_of_shots
    
    #Set Collection of empty dictionaries: columns (r = right) x rows (l = left)
    B = [[Counter({}) for r in range(n+1)] for l in range(k+1)]

    for n in range(0,number_of_shots+1):
        L = min(k,n)
        for l in range(L,-1,-1):
            r = n - l
            if r == 0:
                # If there is as streak of r hits on the left (r for row), and no more shots to the right
                # Then there is just one way for this to happen
                B[l][0]=Counter({(0,0):1})              
            else: #r>0
                if l == k:
                    B[k][r] = update_count_and_probability(1,0,p,B[0][r-1]) \
                    + update_count_and_probability(1,1,p,B[k][r-1])
                else: #l<k
                    B[l][r] =update_count_and_probability(0,0,p,B[0][r-1]) \
                    + update_count_and_probability(0,1,p,B[l+1][r-1])
    return B

# Calulate the expected proportionw with the matrix of dictionaries B as input
def expected_proportion(number_of_shots,streak_length,B):    
    num=0
    den=0
    for number_of_events_of_each_kind, probability \
        in B[0][number_of_shots].items():
        
        if number_of_events_of_each_kind==(0,0):continue
            
        number_of_hits_after_streak_of_hits=number_of_events_of_each_kind[1]
        
        number_of_shots_after_streak_of_hits= \
            number_of_events_of_each_kind[0]+number_of_events_of_each_kind[1]
            
        proportion_hits_after_streak_of_hits= \
            number_of_hits_after_streak_of_hits/number_of_shots_after_streak_of_hits
            
        num += probability*proportion_hits_after_streak_of_hits
        den += probability
    return num/den

# Build a dictionary where the items defined unique proprotions
# key = proportion of hits after streak_length hits
# value = probability of occurance
def histogram_counts(number_of_shots,streak_length,B): 
    histogram = Counter({})           
    for number_of_events_of_each_kind, probability \
        in B[0][number_of_shots].items():
        
        if number_of_events_of_each_kind==(0,0):continue
            
        number_of_hits_after_streak_of_hits=number_of_events_of_each_kind[1]
        
        number_of_shots_after_streak_of_hits= \
            number_of_events_of_each_kind[0]+number_of_events_of_each_kind[1]
            
        proportion_hits_after_streak_of_hits= \
            number_of_hits_after_streak_of_hits/number_of_shots_after_streak_of_hits
            
        histogram = histogram + Counter({proportion_hits_after_streak_of_hits:probability})
    return histogram

# Calulate the weighted expected proportionw with the matrix of dictionaries B as input
def weighted_expected_proportion(number_of_shots,streak_length,sensitivity_to_sample_size,B):    
    num=0
    den=0
    for number_of_events_of_each_kind, probability \
        in B[0][number_of_shots].items():
        
        if number_of_events_of_each_kind==(0,0):continue
            
        number_of_hits_after_streak_of_hits=number_of_events_of_each_kind[1]
        
        number_of_shots_after_streak_of_hits= \
            number_of_events_of_each_kind[0]+number_of_events_of_each_kind[1]
            
        proportion_hits_after_streak_of_hits= \
            number_of_hits_after_streak_of_hits/number_of_shots_after_streak_of_hits
        weight = number_of_shots_after_streak_of_hits**sensitivity_to_sample_size   
        num += weight*probability*proportion_hits_after_streak_of_hits
        den += weight*probability
    return num/den


def probability_propotion_undefined(number_of_shots,streak_length,B):    
    probability_missing = 0

    for number_of_events_of_each_kind, probability \
        in B[0][number_of_shots].items():
        
        if number_of_events_of_each_kind==(0,0):continue
        
        probability_missing += probability

    return probability_missing