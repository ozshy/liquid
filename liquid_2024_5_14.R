#liquid_2024_1_18.R corresponding to liquid_v15.tex
# Numerical simulations for a paper titled: "Required Liquidity, Lobbying, and Bank Competition"

### R packages needed
#library(nleqslv)# package numerical solution of a system of equations
library(ggplot2); theme_set(theme_bw())
library(latex2exp)# LaTeX in ggplot

# parameter definitions
phi1=0.05# failure probability.
phi2=0.10# failure probability. 
tau1=5# transp cost or market power parameter
tau2=6# transp cost or market power parameter
q1=0.1# reserve ratio (later will be drawn as a vector)
(q1.vec = seq(0.1, 0.7, 0.05))
q2=0.3# reserve ratio (later will be drawn as a vector)
(q2.vec = seq(0.2, 0.8, 0.05))
(deltaq.vec = q2.vec-q1.vec)
#
deltav1 = 12 #v2-v1 in paper

n=120# total num of depositors
mu1=0.02# deposit insurance premium
mu2=0.05
R1 = 2# return on bank risky investment/lending
R2 = 3


# Verify Assumption 4
deltav1 - R1*deltaq.vec > 0
deltav1 - R2*deltaq.vec > 0
deltav1 - R1*deltaq.vec < 3*tau1
deltav1 - R2*deltaq.vec < 3*tau1
deltav1 - R1*deltaq.vec < 3*tau2
deltav1 - R2*deltaq.vec < 3*tau2

# Verify Assumption 2
(1-phi1)*R1
(1-phi1)*R1 > 1
(1-phi2)*R1
(1-phi2)*R1 > 1
(1-phi1)*R2
(1-phi1)*R2 > 1
(1-phi2)*R2
(1-phi2)*R2 > 1

# equilibrium deposit interest rates r1, r2, eq (7) in the paper. Also, verify Result 1: r1 > r2 and partly Result 2a, parly Result 3
(r1_eql1.vec = ((1-phi1)*(3-2*q1.vec-q2.vec)*R1 +deltav1*(1-phi1) - 3*(mu1 + tau1*(1-phi1)))/(3*(1-phi1)))
(r2_eql1.vec = ((1-phi1)*(3-q1.vec-2*q2.vec)*R1 -deltav1*(1-phi1) - 3*(mu1 + tau1*(1-phi1)))/(3*(1-phi1)))
# higher phi
(r1_eql2.vec = ((1-phi2)*(3-2*q1.vec-q2.vec)*R1 +deltav1*(1-phi2) - 3*(mu1 + tau1*(1-phi2)))/(3*(1-phi2)))
(r2_eql2.vec = ((1-phi2)*(3-q1.vec-2*q2.vec)*R1 -deltav1*(1-phi2) - 3*(mu1 + tau1*(1-phi2)))/(3*(1-phi2)))
# higher tau
(r1_eql3.vec = ((1-phi1)*(3-2*q1.vec-q2.vec)*R1 +deltav1*(1-phi1) - 3*(mu1 + tau2*(1-phi1)))/(3*(1-phi1)))
(r2_eql3.vec = ((1-phi1)*(3-q1.vec-2*q2.vec)*R1 -deltav1*(1-phi1) - 3*(mu1 + tau2*(1-phi1)))/(3*(1-phi1)))
# higher mu
(r1_eql4.vec = ((1-phi1)*(3-2*q1.vec-q2.vec)*R1 +deltav1*(1-phi1) - 3*(mu2 + tau1*(1-phi1)))/(3*(1-phi1)))
(r2_eql4.vec = ((1-phi1)*(3-q1.vec-2*q2.vec)*R1 -deltav1*(1-phi1) - 3*(mu2 + tau1*(1-phi1)))/(3*(1-phi1)))
# higher R
(r1_eql5.vec = ((1-phi1)*(3-2*q1.vec-q2.vec)*R2 +deltav1*(1-phi1) - 3*(mu1 + tau1*(1-phi1)))/(3*(1-phi1)))
(r2_eql5.vec = ((1-phi1)*(3-q1.vec-2*q2.vec)*R2 -deltav1*(1-phi1) - 3*(mu1 + tau1*(1-phi1)))/(3*(1-phi1)))

# equilibrium xhat eq (8)
(xhat.vec = 1/2 - ((deltav1 - R1*deltaq.vec)/(6*tau1)))
# verify with eq (2)
(xhat.vec = 1/2 + ((-deltav1 + r1_eql1.vec - r2_eql1.vec)/(2*tau1)))

# eq (9) (r difference)
((2*deltav1 + R1*deltaq.vec)/3)
# direct subtraction (verification)
r1_eql1.vec - r2_eql1.vec

# Equilibrium expected profit eq (10) in the paper
(eprofit1.vec = (n*(1-phi1)*(deltav1 - R1*deltaq.vec -3*tau1)^2)/(18*tau1) )
(eprofit2.vec = (n*(1-phi1)*(deltav1 - R1*deltaq.vec +3*tau1)^2)/(18*tau1) )
# verify direct from eq (6)
(d1.vec = n*xhat.vec)
(d2.vec = n*(1-xhat.vec))
#
(eprofit1.vec = (1-phi1)*((1-q1.vec)*d1.vec*R1 - d1.vec*r1_eql1.vec)-mu1*d1.vec)
(eprofit2.vec = (1-phi1)*((1-q2.vec)*d2.vec*R1 - d2.vec*r2_eql1.vec)-mu1*d2.vec)

# difference below eq (10)
(eprofit2.vec-eprofit1.vec)
# verify from text
(2/3)*(n*(1-phi1)*(deltav1 - R1*deltaq.vec))

### Section 5.1: Lobbying delta q
lambda1 = 5 # new lobbying parameter
# eq (15) in paper -delta q =
n*mu1*tau1/(2*lambda1)
# In the paper, I put a maximum q_0 on delta q, no need for simulations.

### Section 4: Lobbying by individual banks
# the following is lower bound on lam eq (12) in paper
(2*n*R1^2*(1-phi1))/(9*tau1)
# try 
lam=120

# eql lobbying eq (17) and (18)
(qq1.vec = (n*R1*(1-phi1)*(n*R1^2*(phi1-1) - 3*lam *(deltav1 - R1*deltaq.vec -3*tau1)))/(6*lam *(9*lam*tau1 -n*R1^2*(1-phi1))) )
(qq2.vec = (n*R1*(1-phi1)*(n*R1^2*(phi1-1) + 3*lam *(deltav1 - R1*deltaq.vec +3*tau1)))/(6*lam *(9*lam*tau1 -n*R1^2*(1-phi1))) )
# difference eq (19)
qq2.vec-qq1.vec
# from (19) verify
(n*R1*(1-phi1)*(deltav1-R1*deltaq.vec)/(9*lam*tau1 - n*R1^2*(1-phi1)))
# show cases where reduction does not yield negative liquidity requirement
q1.vec-qq1.vec
q2.vec-qq2.vec

### Section 5: Lobbying by association
# eq (21)
(q2bar.vec = (n*R1*(1-phi1)*(deltav1 -R1*deltaq.vec))/ (9*lam*tau1 - n*R1^2*(1-phi1)) )
# verify change in liquidity
q2.vec - q2bar.vec

### Section 6.1: extension: unequal lobbying costs
lam # from above where lam1 = lam2
# now define different lam2
#(lam2.vec = seq(80, 160, 1))
#length(lam2.vec)
#(lam1.vec = rep(120, length(lam2.vec)))

(lam2.vec = seq(500, 530, 0.2))
length(lam2.vec)
(lam1.vec = rep(500,length(lam2.vec)))

deltav1
(deltav3 = 4)
R1
(tau3 = 50)
q1
q2
(deltaq = q2-q1)

#verify Assumption 6
lam1.vec > (n*R1^2*(1-phi1))/(18*tau3)
lam2.vec > (n*R1^2*(1-phi1))/(18*tau3)
#
(lam1.vec*lam2.vec) > (n*R1^2*(1-phi1))/(18*tau3)

# extension eql lobbying eq (23) and (24) from paper
(qq1_ext.vec = (n*R1*(1-phi1)*(3*lam2.vec*(3*tau3 -deltav3 +R1*deltaq) -n*R1^2*(1-phi1) ))/(3*(18*lam1.vec*lam2.vec*tau3 -n*R1^2*(1-phi1)*(lam1.vec+lam2.vec))) )
#
(qq2_ext.vec = (n*R1*(1-phi1)*(3*lam1.vec*(3*tau3 +deltav3 -R1*deltaq) -n*R1^2*(1-phi1) ))/(3*(18*lam1.vec*lam2.vec*tau3 -n*R1^2*(1-phi1)*(lam1.vec+lam2.vec))) )
#
## extension eql lobbying from Derive directly to verify no typos in paper
(qq1_ext_verify.vec = n*R1*(1-phi1)*(n*R1^2*(phi1-1)-3*lam2.vec *(q1*R1 -q2*R1 +deltav3-3*tau3))/(3*(n*R1^2*(lam1.vec+lam2.vec) *(phi1-1) +18*lam1.vec*lam2.vec*tau3)))
# test
qq1_ext.vec - qq1_ext_verify.vec
#
(qq2_ext_verify.vec= n*R1*(1-phi1)*(n*R1^2*(phi1-1)+3*lam1.vec *(q1*R1 -q2*R1 +deltav3+3*tau3))/(3*(n*R1^2*(lam1.vec+lam2.vec) *(phi1-1) +18*lam1.vec*lam2.vec*tau3)))
# test
qq2_ext.vec - qq2_ext_verify.vec
#
# verify that qq1 < q1 and qq2 < q2
qq1_ext.vec < q1
qq2_ext.vec < q2

(delta_qq_ext = qq2_ext.vec - qq1_ext.vec)
# verify that typing directly from Derive file
(delta_qq_ext_verify = (n*R1*(1-phi1)*(q1*R1*(lam1.vec +lam2.vec) -q2*R1*(lam1.vec+lam2.vec) +deltav3*(lam1.vec+lam2.vec) +3*tau3*(lam1.vec-lam2.vec)))/(n*R1^2 *(lam1.vec+lam2.vec)*(phi1-1) +18*lam1.vec*lam2.vec*tau3))
# test [almost]
delta_qq_ext - delta_qq_ext_verify

## redo above for higher phi to show a shift in the graph
phi1
phi3 = 0.06

#verify Assumption 6
lam1.vec > (n*R1^2*(1-phi3))/(18*tau3)
lam2.vec > (n*R1^2*(1-phi3))/(18*tau3)
#
(lam1.vec*lam2.vec) > (n*R1^2*(1-phi3))/(18*tau3)

# extension eql lobbying eq (23) and (24) from paper
(qq1_ext2.vec = (n*R1*(1-phi3)*(3*lam2.vec*(3*tau3 -deltav3 +R1*deltaq) -n*R1^2*(1-phi3) ))/(3*(18*lam1.vec*lam2.vec*tau3 -n*R1^2*(1-phi3)*(lam1.vec+lam2.vec))) )
#
(qq2_ext2.vec = (n*R1*(1-phi3)*(3*lam1.vec*(3*tau3 +deltav3 -R1*deltaq) -n*R1^2*(1-phi3) ))/(3*(18*lam1.vec*lam2.vec*tau3 -n*R1^2*(1-phi3)*(lam1.vec+lam2.vec))) )
#
## extension eql lobbying from Derive directly to verify no typos in paper
(qq1_ext_verify2.vec = n*R1*(1-phi3)*(n*R1^2*(phi3-1)-3*lam2.vec *(q1*R1 -q2*R1 +deltav3-3*tau3))/(3*(n*R1^2*(lam1.vec+lam2.vec) *(phi3-1) +18*lam1.vec*lam2.vec*tau3)))
# test
qq1_ext2.vec - qq1_ext_verify2.vec
#
(qq2_ext_verify2.vec= n*R1*(1-phi3)*(n*R1^2*(phi3-1)+3*lam1.vec *(q1*R1 -q2*R1 +deltav3+3*tau3))/(3*(n*R1^2*(lam1.vec+lam2.vec) *(phi3-1) +18*lam1.vec*lam2.vec*tau3)))
# test
qq2_ext2.vec - qq2_ext_verify2.vec
#
# verify that qq1 < q1 and qq2 < q2
qq1_ext2.vec < q1
qq2_ext2.vec < q2

(delta_qq_ext2 = qq2_ext2.vec - qq1_ext2.vec)
# verify that typing directly from Derive file
(delta_qq_ext_verify2 = (n*R1*(1-phi3)*(q1*R1*(lam1.vec +lam2.vec) -q2*R1*(lam1.vec+lam2.vec) +deltav3*(lam1.vec+lam2.vec) +3*tau3*(lam1.vec-lam2.vec)))/(n*R1^2 *(lam1.vec+lam2.vec)*(phi3-1) +18*lam1.vec*lam2.vec*tau3))
# test [almost]
delta_qq_ext2 - delta_qq_ext_verify2

## Making a data frame
qq_ext.df = data.frame(lam1.vec, lam2.vec, qq1_ext.vec, qq2_ext.vec, qq1_ext2.vec, qq2_ext2.vec)
dim(qq_ext.df)
#
# Plot Figure 2 (section 6.1) in the paper
ggplot(qq_ext.df, aes(x=lam2.vec)) +geom_line(aes(y=qq1_ext.vec), linetype="solid", size=1.2) +geom_line(aes(y=qq2_ext.vec), linetype="solid", size=1.2) +geom_line(aes(y=qq1_ext2.vec), linetype="longdash", size=1.2, color="red") +geom_line(aes(y=qq2_ext2.vec), linetype="longdash", size=1.2, color={"red"}) + scale_x_continuous(breaks = seq(500,530,5)) + scale_y_continuous(breaks = seq(0.06, 0.08, 0.001))+labs(x=TeX("Bank 2's lobbying cost parameter: $\\lambda_2$"), y=TeX("Lobbying outcomes: $\\tilde{q}_1$ and $\\tilde{q}_2$"))  +theme(axis.text.x = element_text(size = 14, color = "black"),  axis.text.y = element_text(size = 16, color = "black"), text = element_text(size = 20)) +annotate("text", x = 510, y = 0.0767, label =TeX("$\\tilde{q}_2$"), size = 8, color="black") +annotate("text", x = 510, y = 0.0752, label =TeX("$\\tilde{q}_2$"), size = 8, color="red") +annotate("text", x = 505, y = 0.0744, label =TeX("$\\tilde{q}_1$"), size = 8, color="black") +annotate("text", x = 505, y = 0.0731, label =TeX("$\\tilde{q}_1$"), size = 8, color="red") + geom_segment(aes(x = 510, y = 0.0762, xend = 510, yend = 0.0756), arrow = arrow(length = unit(0.6, "cm")), color="red", size=1) + geom_segment(aes(x = 505, y = 0.0741, xend = 505, yend = 0.0735), arrow = arrow(length = unit(0.6, "cm")), color="red", size=1) + geom_segment(aes(x = 500, y = 0.0728, xend = 500, yend = 0.078), color="blue", size=1.5, linetype = "dotted") +annotate("text", x = 500, y = 0.0726, label =TeX("$\\lambda_1$"), size = 8, color="blue")

# Finding intersections of qq1_ext and qq2_ext for the disussion below Figure 

# verify Assumption 3 (from the main model)
(1-phi1)*R1 > phi1
(1-phi3)*R1 > phi3
# verify Assumption 4 (from the main model
deltav3-R1*deltaq > 0
deltav3-R1*deltaq < 3*tau3
