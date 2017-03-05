# graphics-engine
Extra functionality implemented:

Stochastic LSystems:
	if one replacement rule is given for a character, it's probability is 1
	
	if n replacement rules are given, but there aren't specific probabilities associated with any, 
	their probability is 1/n for each
	
	a fixed probability can be assigned by using 'char -> "replacement rule" : double probability'
	
	if n replacement rules are given and k of them have a fixed probability,
	the probability for the other n-k is (1-total probability of k rules)/(n-k)

[Extra] rainbow:
	true  ~>  draws lines2D in rainbow colors instead of the specified colors
	false ~>  use the normal color specified in [2DLSystem] color or [Figurex] color

