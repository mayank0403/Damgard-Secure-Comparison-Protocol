# Damgard-Secure-Comparison-Protocol
Proof-of-concept implementation of Damgard et al. secure comparison protocol as proposed in "Secure comparison for online auctions" by Ivan Damgard et al.
The code has been written in c++ with Pari c Library for number theoretic operations.<br/> <br/>

Example Parameter values for running the protocol correctly - <br/>
* First enter "1" to use the protocol in 32 bit mode.
* Then input `l = 33`
* Then input `t = 20`.<br/>
You should now get one of the values in the output to be as "1" if `x1>m1`.
