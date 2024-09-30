##############################
#
# David Tompkins - Addendum to the A Exam materials
# A simple simulation of data based on a causal graph
# with recovery of causal effects in a secondary population
#
##############################

#### 0 Loading Packages/setting up ####
if (!require(causaleffect)) install.packages("causaleffect")
if (!require(dosearch)) install.packages("dosearch")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(patchwork)) install.packages("patchwork")
library(dosearch)
library(causaleffect)
library(igraph)
library(ggplot2)
library(dplyr)
library(patchwork)


#### 1 Transporting causal effect  ####
# Will create graph for both original and target populations
original <- graph.formula(A -+ C, B -+ A, B -+ C)
target <- graph.formula(A -+ C, B -+ A, B -+ C, S -+ B)

# manipulations to identify selection as a square node in line with Pearl's work. 
# adapted this bit of code.
target <- set_vertex_attr(target, "description", 4, "S") # important step to identify as selection variable
V(target)$shape <- ifelse(V(target)$name =="S","square","circle") # unimportant step to make the S node square

layout1 <- matrix(c(0,1, 2,1, 1,1.5), ncol=2, byrow=TRUE) # strictly for visualization - these are positional pairs for vertices 
layout2 <- matrix(c(0,1, 2,1, 1,1.5, 1,1.8), ncol=2, byrow=TRUE)

# visualizing both graphs
plot(original,layout=layout1)
plot(target,layout=layout2)

# Calculating transport function
transport_formula = transport(y="C",x="A", D=target, steps =F) # Oddly y before x is their suggested order. Note that original graph is not needed here

print(transport_formula) # Sum_b(P*(C|B,A)P*(B)) 
## This formula is equivalent to that provided by Pearl & Bareinboim, if you apply Rule 2 of Do-Calculus.
## In the mutilated graph Gb_ that trims all connections leaving B, B and C are unrelated (as B only acts upon C directly)
## this means we could also write this as sum_b(P(C | do(A), B) x P*(B)) (as in P&B)

# As the causaleffect package was a bit tricky to navigate, I implemented it also in dosearch, which I believe has gotten a bit more traction (and shares authors)
# the real benefit (to me) in using dosearch is that it readily accepts do() nomenclature
# this meant that the outputted formula was more readily applicable.

# defining graph structure again (dosearch actually accepts igraph as well, but figured it best to show in their language)
target_graph <- "
  A -> C
  B -> A
  B -> C
  S -> B
"

# defining what we care about (what dosearch will attempt to solve for)
query <- "P(C|do(A),S)" #equivalent to P*(C|do(A))

# defining what data is available in this graph
data <- "
  P(B|S)
  P(C|B,do(A))
"
# this is saying we know B in the new domain and The effect of do(A) on C for the distribution of B in the source domain.

transport_formula2 <- dosearch(data,query,target_graph,transportability="S")
print(transport_formula2) #  Sum_b(P(B|S) x P(C |do(A),B)) - this is exactly as in P&B without manipulations. (P*(B) and P(B|S) are by definition equivalent)


##### 2 Simulating Data ####

# Structure of graph is A <- B -> C; A -> C (picked ABC over XYZ arbitrarily)
# Will generate B, then A as function of B, and C as function of A and B.
# For experimental group, will generate B, set A (atomic intervention), generate C
# We will say A is a treatment of 'puzzle exposure', B is of wealth, C is of spatial skill.
# will assume normal distributions 


set.seed(1234) # setting seed to keep results consistent between runs
n = 500 # number of participants per group
B_A = 0.75 # direct causal relation from B to A - note these are differed in the multi-iteration simulation.
B_C = 0.75 # direct causal relation from B to C
A_C = 0.75 # direct causal relation from A to C

experimental_data <- data.frame(
  condition = c(rep('experimental',n),rep('control',n)),
  puzzle = c(rep(1,n),rep(0,n)),
  wealth = rnorm(2*n, mean=0, sd=1),
  spatial = rnorm(2*n, mean=0, sd=1)
) %>% # being a bit cheeky here, but it's efficient!
  mutate(spatial = spatial + (A_C*puzzle)+(B_C*wealth))

# Visualing group differences:
ggplot(data=experimental_data, aes(x=condition,y=spatial))+
  geom_violin(aes(fill=condition),alpha=0.4)+
  geom_jitter(width=0.2,height=0)

ggplot(data=experimental_data, aes(x=wealth,y=spatial,color=condition))+
  geom_jitter(alpha=0.3)+
  geom_smooth(method="lm")

# simulating target domain (ordering here is for clearer narrative, would otherwise generate alongside experimental_data)
target_data <- data.frame(
  #no condition here
  puzzle = rnorm(2*n, mean = 0, sd=1), #arguably odd that here our values here include negatives, when our assignments previously were 0,1
  wealth = rnorm(2*n, mean = 1, sd=0.5), # different than in original domain. Richer and more homogenous 
  spatial = rnorm(2*n, mean = 0, sd=1)
) %>% # same trick again
  mutate(puzzle = puzzle + B_A*wealth) %>% # applying wealth effect on puzzle
  mutate(spatial =  (A_C*puzzle)+(B_C*wealth))

#### Estimating and Transporting Causal effect A->C ####
## Estimating P(C|do(A)|B) - this model is the estimate effectively
experimental_model <- lm(spatial~puzzle*wealth, data=experimental_data) # the relation we've built above is actually a simple linear relation without interaction, 
                                                                        # but it woudl be pretty reasonable to think there should be an interaction. Running without
                                                                        # the interaction will give a pretty similar result, but of course in that case you don't need 
                                                                        # to calculate individual effects below - each person's diff between expected and real would be identical

#checking recovered coefficients in experimental model
summary(experimental_model)
experimental_model$coefficients["puzzle"] # estimated at 0.78 - not too far off!

## estimating P*(B) or P(B|S) (again these are equivalent)
target_data$density_wealth <- dnorm(target_data$wealth, 
                               mean = mean(target_data$wealth), 
                               sd = sd(target_data$wealth))

# using model to predict new spatial for target, then weighting by B
target_data$predicted_exp <- predict(experimental_model, newdata = data.frame(puzzle=1, wealth=target_data$wealth))
target_data$predicted_con <- predict(experimental_model, newdata = data.frame(puzzle=0, wealth=target_data$wealth)) # predicting each condition


#calculating P*(B) x P(C|do(A))

target_data <- target_data %>% # normally would create a new df, but for cleanness of inspection, editing in place
  mutate(individual_effect = predicted_exp - predicted_con) %>% # no barriers to counterfactuals in simulation  
  mutate(weighted_effect= individual_effect * density_wealth)
  
target_summary <- target_data %>%
  ungroup() %>% #best practice imo
  summarise(sum_density = sum(density_wealth),
            sum_weighted_effect = sum(weighted_effect)) %>%
  mutate(average_causal_effect = sum_weighted_effect / sum_density)

target_summary$average_causal_effect
  # average estimated causal effect in target population  = 0.82
  # This is pretty reasonable in the grand scheme. If you crank n up, it will 
  # get closer to true value of 0.75 (In this simple example, the real relation is perfectly linear,
  # and the causal effect in target domain is orthogonal to B-wealth). In the below multi-sample 
  # simulation, I've added an interaction which leads to a different 'true' value that is dependent
  # on P*(b)

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#--- splitting without making a section
#### 3 Running again with interaction and for many seeds ####
# This has a lot of code reuse. A cleaner route would be using functions, etc.
# However, in my opinion, that makes it less intelligible for understanding the 
# mechanisms involved.

results <- data.frame(recovered_effect = numeric(), trasported_effect = numeric())
for (i in 1:1000){

set.seed(i) # setting seed to keep results consistent between runs
n = 500 # number of participants per group
B_A = 0.2 # direct causal relation from B to A 
B_C = 0.3 # direct causal relation from B to C
A_C = 0.4 # direct causal relation from A to C - These all varied to keep it interesting
interaction_modifier = 0.2 

experimental_data <- data.frame(
  condition = c(rep('experimental',n),rep('control',n)),
  puzzle = c(rep(1,n),rep(0,n)),
  wealth = rnorm(2*n, mean=0, sd=1),
  spatial = rnorm(2*n, mean=0, sd=1)
) %>% #
  mutate(spatial = spatial + (A_C*puzzle)+(B_C*wealth) +(interaction_modifier*puzzle*wealth) )


# simulating target domain (ordering here is for clearer narrative, would otherwise generate alongside experimental_data)
target_data <- data.frame(
  #no condition here
  puzzle = rnorm(2*n, mean = 0, sd=1), # arguably odd that here our values here include negatives, when our assignments previously were 0,1
  wealth = rnorm(2*n, mean = 1.5, sd=0.5), # different than in original domain. Richer and more homogenous 
  spatial = rnorm(2*n, mean = 0, sd=1)
) %>% 
  mutate(puzzle = puzzle + B_A*wealth) %>% # applying wealth effect on puzzle
  mutate(spatial =  (A_C*puzzle)+(B_C*wealth)+(interaction_modifier*puzzle*wealth))

#### Estimating and Transporting Causal effect A->C ####
## Estimating P(C|do(A),B) - this model is the estimate,   effectively
experimental_model <- lm(spatial~puzzle*wealth, data=experimental_data) 
 

## estimating P*(B) or P(B|S) (again these are equivalent)
target_data$density_wealth <- dnorm(target_data$wealth, 
                                    mean = mean(target_data$wealth), 
                                    sd = sd(target_data$wealth))

# using model to predict new spatial for target, then weighting by B
target_data$predicted_exp <- predict(experimental_model, newdata = data.frame(puzzle=1, wealth=target_data$wealth))
target_data$predicted_con <- predict(experimental_model, newdata = data.frame(puzzle=0, wealth=target_data$wealth)) # predicting each condition


#calculating P*(B) x P(C|do(A),B)
target_data <- target_data %>% # normally would create a new df, but for cleanness of inspection, editing in place
  mutate(individual_effect = predicted_exp - predicted_con,
         individual_effect2 = experimental_model$coefficients["puzzle"] + (experimental_model$coefficients["puzzle:wealth"] * wealth)) %>% # GPT o1 suggested this simpler route when I asked for graphing advice  
                                                                                                                                           # This is honestly a much cleaner and smarter way to do this. I've left 
                                                                                                                                           # my human logic in above, but they are effectively identical in output.
                                                                                                                                           # despite my politeness to the AI solution - it would be rather involved in more complicated models, and relies on how contrasts are set (i.e. the effect being from 0 to 1). For that reason, probably better to use my original method.

  mutate(weighted_effect= individual_effect * density_wealth)


target_summary <- target_data %>%
  ungroup() %>% #best practice imo
  summarise(sum_density = sum(density_wealth),
            sum_weighted_effect = sum(weighted_effect),
            transported_effect2 = mean(individual_effect2)) %>%
  mutate(average_causal_effect = sum_weighted_effect / sum_density) #GPT o1's recomendation (as in line 209)

avg_wealth=mean(experimental_data$wealth) # not needed in model without interaction, but if we're seeking the total effect, very important


results <- rbind(results,data.frame(recovered_effect = (experimental_model$coefficients["puzzle"]+(experimental_model$coefficients["puzzle:wealth"]*avg_wealth)), 
                                    transported_effect = target_summary$average_causal_effect,
                                    transported_effect2 =target_summary$transported_effect2,
                                    iter = i))
}

# Calculate the true transported causal effect - this is basically the function used to generate spatial score above, excluding the base distribution
calculate_true_effect <- function(puzzle, wealth) {
  return((A_C * puzzle) + (B_C * wealth) + (interaction_modifier * puzzle * wealth))
}

# Generate wealth values for the target population
n <- 5000000  # number of participants per group - large number to get 'real' value easily. 
set.seed(42)  # Set a seed for reproducibility
target_wealth <- rnorm(2*n, mean = 1.5, sd = 0.5)

# Calculate the true transported causal effect, using access to counterfactuals to calculate individual effects
true_effect_puzzle_1 <- mean(calculate_true_effect(1, target_wealth))
true_effect_puzzle_0 <- mean(calculate_true_effect(0, target_wealth))
true_transported_effect <- true_effect_puzzle_1 - true_effect_puzzle_0 #0.7


mean(results$transported_effect) # 0.695 - very good!
mean(results$transported_effect2) # the same

# In this case our fit is better than realistic, as we have full knowledge of the 
# causal structure and as the relations will extrapolate indefinitely. Neither of these
# are realistic, sadly, and applying in practice would likely yield a much worse estimate 
# due to additional causal relations, outside-of-distribution values, etc.



#### 4. Clean graphs for paper ####
# Nothing new here, just cleaning up graphs. Normally would have just included above, but 
# quick graphs I think are better for explaining in code than these multi-line graphs
# Claude Sonnet 3.5 helped with aesthetics. I refrained from using AI to write code above for understanding's sake
# (Though, I did bounce ideas and work through the logic with GPT 1o-preview and Sonnet 3.5)


# Create the plot
plot_data <- results %>%
  select(iter, recovered_effect, transported_effect) %>%
  tidyr::pivot_longer(cols = c(recovered_effect, transported_effect),
                      names_to = "effect_type",
                      values_to = "effect_value")

# Create a named vector for label mapping
effect_labels <- c(
  recovered_effect = "Recovered Effect",
  transported_effect = "Transported Effect"
)

ggplot(plot_data, aes(x = effect_type, y = effect_value, fill = effect_type)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.2, height =0, alpha = 0.4, size = 1) +
  geom_segment(aes(x = 1.7, xend = 2.3, y = true_transported_effect, yend = true_transported_effect),
               color = "red", alpha = 0.7, linewidth = 1) + # this segment makes the graphing slow! But it's cleaner than geom_hline 
  geom_segment(aes(x = 0.7, xend = 1.3, y = 0.4, yend = 0.4),
               color = "red", alpha = 0.7, linewidth = 1) + # this segment makes the graphing slow! But it's cleaner than geom_hline 
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "darkred") +
  labs(title = "Comparison of Recovered and Transported Effects",
       x = "Effect Type",
       y = "Effect Value") +
  scale_x_discrete(labels = effect_labels) +
  scale_fill_manual(values = c("recovered_effect" = "#1f77b4", "transported_effect" = "#ff7f0e"),
                    labels = effect_labels) +
  theme_minimal() +
  theme(legend.position = "bottom")


### Claude had suggested these latter plots, but I actually prefer ^ this one. 
### Left these for interest's sake.

# 1. Scatter plot with smoothed line
p1 <- ggplot(data = results, aes(x = recovered_effect, y = transported_effect)) +
  geom_point(alpha = 0.2, color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Recovered vs Transported Effect",
       x = "Recovered Effect",
       y = "Transported Effect") +
  theme_minimal()

# 2. Density plot of differences
results <- results %>%
  mutate(difference = transported_effect - recovered_effect)

p2 <- ggplot(data = results, aes(x = difference)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(difference)), color = "red", linetype = "dashed") +
  labs(title = "Density of Differences",
       x = "Transported Effect - Recovered Effect",
       y = "Density") +
  theme_minimal()

# 3. Boxplot of effects
p3 <- results %>%
  tidyr::pivot_longer(cols = c(recovered_effect, transported_effect), names_to = "effect_type", values_to = "effect") %>%
  ggplot(aes(x = effect_type, y = effect, fill = effect_type)) +
  geom_violin(alpha=0.3) + #geom_boxplot() + #claude suggested boxplot, but violins are better imo (esp. with such similarity)
  labs(title = "Distribution of Effects",
       x = "Effect Type",
       y = "Effect Size") +
  theme_minimal() +
  theme(legend.position = "none")


# Combine all plots
combined_plot <- p1 + p2 + p3
combined_plot <- combined_plot +
  plot_annotation(
    title = "Summary of Simulation Results: Recovered vs Transported Effects across 1000 iterations",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

# Display the combined plot
print(combined_plot) # these are the Claude suggested comparisons - but I don't think they're as interesting as a simple geom_violin with geom_points.
# left in for curiosity's sake, but not 'transported' over to the paper

# thanks for reading!
